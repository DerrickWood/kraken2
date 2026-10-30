#!/usr/bin/env bash
# End-to-end tests for an installed Kraken 2.
#
# Usage: tests/e2e.sh INSTALL_DIR [WORK_DIR]
#
# INSTALL_DIR is a directory produced by install_kraken2.sh.  WORK_DIR holds
# every database, output file and log the run produces; it defaults to a new
# temporary directory and is kept afterwards so CI can upload it.
#
# The suite builds a tiny database from data/ (38 viral sequences, 40 taxa),
# classifies frozen simulated reads from tests/data in several input modes and
# compares results with golden files in tests/golden.  It needs no network
# access and finishes in well under a minute.
#
# Every check prints one line: PASS, FAIL, XFAIL (a known bug reproduced, which
# does not fail the run) or XPASS (a known bug that no longer reproduces).  The
# script exits non-zero if anything FAILs or XPASSes.  An XPASS means a bug was
# fixed, so the corresponding xfail_check call should become a plain check.
#
# Environment:
#   E2E_UPDATE_GOLDEN=1  overwrite tests/golden instead of comparing to it
#   E2E_XPASS_OK=1       report XPASS without failing the run
#   E2E_SKIP_K2=1        skip the k2 front end
#   E2E_PYTHON_PRELOAD=LIB  run Python with LD_PRELOAD=LIB; the sanitizer build
#                        sets this to the ASan runtime so k2 can load an
#                        instrumented libtax.so
#   E2E_ONLY=REGEX       run only checks whose name matches REGEX

set -uo pipefail
export LC_ALL=C

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  echo "Usage: $0 INSTALL_DIR [WORK_DIR]" >&2
  exit 64
fi

TESTS=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$TESTS/.." && pwd)
INST=$(cd "$1" && pwd) || { echo "No such install dir: $1" >&2; exit 64; }
if [ $# -eq 2 ]; then
  mkdir -p "$2"
  WORK=$(cd "$2" && pwd)
else
  WORK=$(mktemp -d "${TMPDIR:-/tmp}/k2e2e.XXXXXX")
fi
GOLDEN="$TESTS/golden"
READS="$TESTS/data"
ACTUAL="$WORK/actual"
LOGS="$WORK/logs"
mkdir -p "$ACTUAL" "$LOGS"

KRAKEN2="$INST/kraken2"
KBUILD="$INST/kraken2-build"
KINSPECT="$INST/kraken2-inspect"
K2="$INST/k2"
PYTHON=${PYTHON:-python3}
if [ -n "${E2E_PYTHON_PRELOAD:-}" ]; then
  mkdir -p "$WORK/bin"
  printf '#!/bin/sh\nLD_PRELOAD="%s" exec %s "$@"\n' "$E2E_PYTHON_PRELOAD" "$(command -v "$PYTHON")" > "$WORK/bin/python"
  chmod +x "$WORK/bin/python"
  PYTHON="$WORK/bin/python"
fi

# Accuracy thresholds.  The fixture reads are error-light fragments of the
# reference genomes, so a working classifier gets essentially all of them.
MIN_EXACT_FRAC=0.97       # simulated reads called at exactly their true taxid
MIN_LINEAGE_FRAC=0.99     # simulated reads called at true taxid or an ancestor
MAX_RANDOM_CLASSIFIED=2   # of the 20 random-sequence pairs

n_pass=0 n_fail=0 n_xfail=0 n_xpass=0 n_skip=0
RESULTS="$WORK/results.tsv"
: > "$RESULTS"

say() { printf '%s\n' "$*"; }

record() {
  # record STATUS NAME DETAIL
  printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$RESULTS"
  if [ -n "$3" ]; then say "$1 $2: $3"; else say "$1 $2"; fi
}

wanted() {
  [ -z "${E2E_ONLY:-}" ] || [[ "$1" =~ $E2E_ONLY ]]
}

# check NAME CMD...  runs CMD; exit status 0 is a pass.  Output of CMD goes to
# logs/NAME.log and the tail is echoed on failure.
check() {
  local name=$1; shift
  wanted "$name" || return 0
  local log="$LOGS/$name.log"
  if "$@" > "$log" 2>&1; then
    n_pass=$((n_pass + 1)); record PASS "$name" ""
  else
    n_fail=$((n_fail + 1)); record FAIL "$name" "see logs/$name.log"
    show_log "$log"
  fi
}

# xfail_check NAME REASON CMD...  like check, but CMD is expected to fail
# because of a known bug described by REASON.
xfail_check() {
  local name=$1 reason=$2; shift 2
  wanted "$name" || return 0
  local log="$LOGS/$name.log"
  # Seen by golden(), which is called from inside CMD.
  local in_xfail=1
  if "$@" > "$log" 2>&1; then
    n_xpass=$((n_xpass + 1))
    record XPASS "$name" "expected to fail ($reason) but passed; if the bug is fixed, change xfail_check to check in tests/e2e.sh"
  else
    n_xfail=$((n_xfail + 1)); record XFAIL "$name" "known bug: $reason"
  fi
}

show_log() {
  # Echo a log indented, eliding the middle of long ones.
  local n
  n=$(wc -l < "$1")
  if [ "$n" -le 30 ]; then
    sed 's/^/    | /' "$1"
  else
    sed -n '1,15s/^/    | /p' "$1"
    echo "    | ... ($((n - 30)) lines omitted; full log in logs/$(basename "$1"))"
    tail -n 15 "$1" | sed 's/^/    | /'
  fi
}

skip() {
  wanted "$1" || return 0
  n_skip=$((n_skip + 1)); record SKIP "$1" "$2"
}

# golden NAME FILE  compares FILE with tests/golden/NAME, or refreshes the
# golden copy when E2E_UPDATE_GOLDEN=1.  Golden files used by an xfail_check
# hold the output expected once the known bug is fixed, so they are refreshed
# only when E2E_UPDATE_XFAIL_GOLDEN=1 as well, from a build that has the fix.
golden() {
  local name=$1 file=$2
  cp "$file" "$ACTUAL/$name"
  if [ "${E2E_UPDATE_GOLDEN:-0}" = 1 ]; then
    if [ "${in_xfail:-0}" = 1 ] && [ "${E2E_UPDATE_XFAIL_GOLDEN:-0}" != 1 ]; then
      echo "not updating tests/golden/$name: it belongs to an expected-failure check (set E2E_UPDATE_XFAIL_GOLDEN=1 to update it from a build with the fix)"
      return 0
    fi
    cp "$file" "$GOLDEN/$name"
    return 0
  fi
  if [ ! -f "$GOLDEN/$name" ]; then
    echo "missing golden file tests/golden/$name; run tests/update_golden.sh"
    return 1
  fi
  if ! diff -u --label "tests/golden/$name" --label "new/$name" "$GOLDEN/$name" "$file" > "$ACTUAL/$name.diff"; then
    echo "differs from tests/golden/$name ($(grep -c '^[-+][^-+]' "$ACTUAL/$name.diff") changed lines)"
    case $name in
      *.out)
        # Summarize per-read call changes before the raw diff.
        awk -F '\t' 'NR == FNR { c[$2] = $1 $3; next }
             ($2 in c) && c[$2] != $1 $3 { n++; if (n <= 10) print "  call changed: " $2 " golden=" c[$2] " now=" $1 $3 }
             END { print "  " n + 0 " reads changed classification status or taxid" }' \
          "$GOLDEN/$name" "$file" ;;
    esac
    head -n 40 "$ACTUAL/$name.diff"
    return 1
  fi
  rm -f "$ACTUAL/$name.diff"
}

same() {
  # same A B  succeeds when files are byte-identical, else shows a diff.
  if ! cmp -s "$1" "$2"; then
    echo "$1 and $2 differ:"
    diff -u "$1" "$2" | head -n 30
    return 1
  fi
}

must_fail() {
  # must_fail CMD...  succeeds when CMD exits with a clean non-zero status:
  # not success, not death by signal (exit status 129-159) and no sanitizer
  # report, since ASan exits with status 1 by default.
  local out rc
  out=$("$@" 2>&1)
  rc=$?
  printf '%s\n' "$out" | tail -n 20
  echo "exit status $rc"
  if printf '%s\n' "$out" | grep -q -e 'AddressSanitizer' -e 'runtime error:'; then
    echo "sanitizer report in output"
    return 1
  fi
  [ $rc -ne 0 ] && { [ $rc -lt 129 ] || [ $rc -gt 159 ]; }
}

with_timeout() {
  # with_timeout SECS CMD...  runs CMD in its own process group and kills the
  # whole group if it outlives SECS, returning 124 like timeout(1).
  perl -e '
    my $secs = shift;
    my $pid = fork() // die "fork: $!";
    if ($pid == 0) { setpgrp(0, 0); exec @ARGV or die "exec: $!" }
    local $SIG{ALRM} = sub { kill "KILL", -$pid; waitpid($pid, 0); print STDERR "timed out after ${secs}s\n"; exit 124 };
    alarm $secs;
    waitpid($pid, 0);
    exit($? & 127 ? 128 + ($? & 127) : $? >> 8);
  ' "$@"
}

# Truth table: read id -> true taxid, parsed from the simulator's headers
# ("@FRAG.1/1 kraken:taxid|10710|NC_001416.1:20345") and taxid 0 for the
# random-sequence pairs.
TRUTH="$WORK/truth.tsv"
awk 'NR % 4 == 1 {
       id = substr($1, 2); sub(/\/[12]$/, "", id)
       t = 0
       if (match($2, /taxid\|[0-9]+\|/)) t = substr($2, RSTART + 6, RLENGTH - 7)
       print id "\t" t
     }' "$READS/sim_1.fq" > "$TRUTH"

# Parent table from the fixture taxonomy, used for lineage checks.
PARENTS="$WORK/parents.tsv"
awk -F '\t[|]\t' '{ print $1 "\t" $2 }' "$REPO/data/nodes.dmp" > "$PARENTS"

# accuracy KRAKEN_OUTPUT  checks per-read calls against the truth table.
accuracy() {
  awk -v min_exact="$MIN_EXACT_FRAC" -v min_lineage="$MIN_LINEAGE_FRAC" \
      -v max_rand="$MAX_RANDOM_CLASSIFIED" -F '\t' '
    FILENAME == ARGV[1] { parent[$1] = $2; next }
    FILENAME == ARGV[2] { truth[$1] = $2; next }
    {
      id = $2; sub(/\/[12]$/, "", id); call = $3; t = truth[id]
      # A chimera is a tie only when both mates are classified together.
      if (id ~ /^CHIMERA/ && $4 !~ /[|]/) next
      if (!(id in truth)) { print "read not in truth table: " id; bad = 1; next }
      if (t == 0) { nr++; if (call != 0) rc++; next }
      ns++
      if (call == t) exact++
      ok = 0
      for (x = t; x != "" ; x = parent[x]) { if (x == call) { ok = 1; break } if (x == 1) break }
      if (call != 0 && ok) lineage++
      else if (wrong_shown++ < 5) print "wrong call: " id " true=" t " called=" call
    }
    END {
      ef = exact / ns; lf = lineage / ns
      printf "simulated=%d exact=%d (%.4f) lineage=%d (%.4f) random=%d random_classified=%d\n", ns, exact, ef, lineage, lf, nr, rc
      if (bad || ef < min_exact || lf < min_lineage || rc > max_rand) {
        print "thresholds: exact>=" min_exact " lineage>=" min_lineage " random_classified<=" max_rand
        exit 1
      }
    }' "$PARENTS" "$TRUTH" "$1"
}

# hit_fraction KRAKEN_OUTPUT  prints the fraction of k-mers with a non-zero
# LCA across all reads (ambiguous "A:" and mate separators excluded).
hit_fraction() {
  awk -F '\t' '{
      n = split($5, f, " ")
      for (i = 1; i <= n; i++) {
        if (f[i] == "|:|") continue
        split(f[i], p, ":"); if (p[1] == "A") continue
        tot += p[2]; if (p[1] != 0) hit += p[2]
      }
    } END { printf "%.4f\n", hit / tot }' "$1"
}

# build_db refuses more threads than omp_get_max_threads(), which is 1 in
# builds without OpenMP (they use stub omp_* functions) and the core count
# otherwise.  Probe it so builds use as many threads as allowed, up to 4.
build_threads_ok() {
  local out
  out=$("$INST/build_db" -p "$1" < /dev/null 2>&1)
  ! printf '%s' "$out" | grep -q 'only wants you to use'
}
BT=1 BT4=1
build_threads_ok 2 && BT=2
build_threads_ok 4 && BT4=4
[ $BT4 = 1 ] && [ $BT = 2 ] && BT4=2
OPENMP=0; [ $BT = 2 ] && OPENMP=1

# Timing helper: seconds since epoch with fractional part, portable to macOS.
now() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }
T0=$(now)

say "Kraken 2 end-to-end tests"
say "  install: $INST"
say "  work:    $WORK"
say ""

cd "$WORK" || exit 1

# ---------------------------------------------------------------------------
# Inputs derived from the frozen reads.
R1="$READS/sim_1.fq"
R2="$READS/sim_2.fq"
awk 'NR % 4 == 1 { print ">" substr($0, 2) } NR % 4 == 2 { print }' "$R1" > r1.fa
gzip -c -n "$R1" > r1.fq.gz
gzip -c -n "$R2" > r2.fq.gz
paste -d '\n' <(paste - - - - < "$R1") <(paste - - - - < "$R2") | tr '\t' '\n' > interleaved.fq

# The reference library, concatenated in sorted order so the built hash table
# is byte-for-byte reproducible.
LIB_FILES=("$REPO"/data/*.fa)
cat "${LIB_FILES[@]}" > lib.fa

new_db_dir() {
  rm -rf "$1"; mkdir -p "$1/taxonomy"
  cp "$REPO/data/nodes.dmp" "$REPO/data/names.dmp" "$1/taxonomy/"
}

# ---------------------------------------------------------------------------
# Database build with kraken2-build.

# Every build pins the sub-block size.  Its default is block size / threads,
# and sub-block boundaries set the insertion order into the compact hash
# table, so without this the database bytes (and a few calls) would depend on
# the thread count, which varies with the OpenMP support of the build.
SUBBLOCK=8192

build_db_kraken2_build() {
  # build_db_kraken2_build DIR THREADS [extra build args]
  local dir=$1 threads=$2; shift 2
  new_db_dir "$dir"
  "$KBUILD" --db "$dir" --add-to-library lib.fa &&
    "$KBUILD" --db "$dir" --build --threads "$threads" --subblock-size $SUBBLOCK "$@"
}

t=$(now)
check build.kraken2-build build_db_kraken2_build db "$BT"
BUILD_SECONDS=$(perl -e 'printf "%.2f", $ARGV[1] - $ARGV[0]' "$t" "$(now)")
if [ ! -s db/hash.k2d ]; then
  say "Database build failed; nothing else can run."
  tail -n 40 "$LOGS/build.kraken2-build.log"
  exit 1
fi

check build.db_check "$PYTHON" "$TESTS/k2db_check.py" db
opts_same() {
  # opts_same A B  compares two opts.k2d files field by field.  The files are
  # a raw IndexOptions struct whose padding bytes are uninitialized, so a
  # byte comparison can differ between otherwise identical builds.
  "$PYTHON" - "$1" "$2" <<'PY'
import struct, sys
fmt = "<QQQQ?7xQiii4x"
a, b = (struct.unpack(fmt, open(p, "rb").read()) for p in sys.argv[1:3])
print("opts fields:", a)
sys.exit(a != b)
PY
}

deterministic_build() {
  local n f
  for n in 1 "$BT4"; do
    build_db_kraken2_build db_t$n "$n" || return 1
    for f in hash.k2d taxo.k2d; do cmp "db/$f" "db_t$n/$f" || return 1; done
    opts_same db/opts.k2d "db_t$n/opts.k2d" || return 1
  done
  echo "hash.k2d and taxo.k2d identical, opts.k2d fields identical, at --threads 1, $BT and $BT4"
}
if [ $OPENMP = 1 ]; then
  check build.deterministic_threads deterministic_build
else
  skip build.deterministic_threads "built without OpenMP"
fi

# Piped through cat because of the redirect bug checked just below.
inspect_db() { "$KINSPECT" --db db | cat > inspect.txt && golden inspect.txt inspect.txt; }
check inspect.golden inspect_db

inspect_to_file() {
  # dump_table writes the header with std::cout, then reopens /dev/fd/1 with
  # an ofstream for the per-taxon lines.  On Linux that open truncates a
  # regular file, erasing the header.
  "$KINSPECT" --db db > inspect_file.txt || return 1
  grep '^#' inspect_file.txt
  grep -q '^# Database options' inspect_file.txt
}
if [ "$(uname -s)" = Linux ]; then
  xfail_check inspect.redirect_keeps_header "on Linux, kraken2-inspect > file loses the header lines (dump_table.cc opens /dev/fd/1 with truncation)" \
    inspect_to_file
else
  check inspect.redirect_keeps_header inspect_to_file
fi

# ---------------------------------------------------------------------------
# Classification with the kraken2 wrapper.

classify() {
  # classify NAME [kraken2 args]  writes NAME.out and NAME.report
  local name=$1; shift
  "$KRAKEN2" --db db --output "$name.out" --report "$name.report" "$@"
}

paired_basic() {
  classify paired --paired "$R1" "$R2" --threads 2 &&
    golden paired.out paired.out && golden paired.report paired.report
}
check classify.paired.golden paired_basic
check classify.paired.accuracy accuracy paired.out

single_basic() {
  classify single "$R1" && golden single.out single.out && golden single.report single.report
}
check classify.single_fastq.golden single_basic
check classify.single_fastq.accuracy accuracy single.out

# A read shorter than k has no minimizers and an empty hit list.  It goes
# first, so the only thread starts on it.
short_first_read() {
  { printf '@SHORT.1\nACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJ\n'; cat "$R1"; } > short_first.fq &&
    "$KRAKEN2" --db db --threads 1 --output short_first.out short_first.fq &&
    [ "$(head -n 1 short_first.out)" = "$(printf 'U\tSHORT.1\t0\t20\t0:0')" ] &&
    cmp <(tail -n +2 short_first.out) single.out
}
check classify.single_short_first_read short_first_read
check classify.single_fasta.matches_fastq bash -c \
  "'$KRAKEN2' --db db --output single_fa.out r1.fa && cmp single.out single_fa.out"
check classify.paired_gzip.matches_plain bash -c \
  "'$KRAKEN2' --db db --paired --output paired_gz.out r1.fq.gz r2.fq.gz && cmp paired.out paired_gz.out"
check classify.paired_gzip_flag.matches_plain bash -c \
  "'$KRAKEN2' --db db --paired --gzip-compressed --output paired_gzf.out r1.fq.gz r2.fq.gz && cmp paired.out paired_gzf.out"
check classify.interleaved.matches_paired bash -c \
  "'$INST/classify' -H db/hash.k2d -t db/taxo.k2d -o db/opts.k2d -S -O interleaved.out interleaved.fq && cmp paired.out interleaved.out"

report_variants() {
  classify mm --paired "$R1" "$R2" --report-minimizer-data &&
    golden paired.minimizer.report mm.report &&
  classify mpa --paired "$R1" "$R2" --use-mpa-style &&
    golden paired.mpa.report mpa.report &&
  classify zero --paired "$R1" "$R2" --report-zero-counts &&
    golden paired.zero.report zero.report &&
  classify names --paired "$R1" "$R2" --use-names &&
    golden paired.names.out names.out
}
check classify.report_variants.golden report_variants

variant() {
  # variant GOLDEN_PREFIX NAME [kraken2 args]  classifies the paired reads with
  # extra options and compares output and report with golden files.
  local prefix=$1 name=$2; shift 2
  classify "$name" --paired "$@" "$R1" "$R2" &&
    golden "$prefix.out" "$name.out" && golden "$prefix.report" "$name.report"
}
check classify.confidence_0.5.golden variant paired.conf0.5 c05 --confidence 0.5
check classify.min_hit_groups_5.golden variant paired.mhg5 mhg5 --minimum-hit-groups 5

# Every simulated pair clears --minimum-hit-groups 5, so a threshold that is
# not applied goes unnoticed above.  tests/data/mhg_*.fq hold shorter pairs,
# several of which fall short of 10 hit groups.
mhg_reads() {
  classify mhg10 --paired --minimum-hit-groups 10 "$READS/mhg_1.fq" "$READS/mhg_2.fq" &&
    golden mhg10.out mhg10.out && golden mhg10.report mhg10.report
}
check classify.min_hit_groups_10.golden mhg_reads
check classify.quick.golden variant paired.quick quick --quick
check classify.min_base_quality_20.golden variant paired.mbq20 mbq --minimum-base-quality 20

check classify.confidence_monotone bash -c "
  a=\$(grep -c '^C' paired.out); b=\$(grep -c '^C' c05.out)
  echo \"classified at confidence 0: \$a, at 0.5: \$b\"; [ \$b -le \$a ] && [ \$b -gt 0 ]"

classified_split() {
  "$KRAKEN2" --db db --paired --output /dev/null \
    --classified-out 'cl#.fq' --unclassified-out 'un#.fq' "$R1" "$R2" || return 1
  local nu nc
  nu=$(grep -c '^@RANDOM' un_1.fq); nc=$(awk 'NR % 4 == 1' cl_1.fq | wc -l | tr -d ' ')
  echo "unclassified-out has $nu RANDOM reads of $(awk 'NR % 4 == 1' un_1.fq | wc -l | tr -d ' '); classified-out has $nc reads"
  [ "$(awk 'NR % 4 == 1' un_1.fq | wc -l)" -eq "$(grep -c '^U' paired.out)" ] &&
    [ "$nc" -eq "$(grep -c '^C' paired.out)" ] && cmp <(awk 'NR % 4 == 1' un_1.fq | cut -c2- | cut -d/ -f1 | cut -d' ' -f1) <(awk -F '\t' '$1 == "U" { print $2 }' paired.out)
}
check classify.classified_unclassified_out classified_split

threads_det() {
  local n
  for n in 1 3 4 8; do
    "$KRAKEN2" --db db --paired --threads $n --output t$n.out --report t$n.report --report-minimizer-data "$R1" "$R2" || return 1
    same paired.out t$n.out || return 1
    same mm.report t$n.report || return 1
  done
  echo "--output and --report-minimizer-data identical at --threads 1, 2, 3, 4, 8"
}
check classify.thread_determinism threads_det
check classify.memory_mapping.matches bash -c \
  "'$KRAKEN2' --db db --paired --memory-mapping --output mmap.out '$R1' '$R2' && cmp paired.out mmap.out"

# With several threads, classify writes output blocks in input order, and a
# block that finishes early (or whose thread finds another thread writing) is
# spilled to a temporary file and appended later.  The fixture reads fit in
# one 8 MB input block, so the spill checks repeat them 400 times (176,000
# pairs, 5 blocks per mate) with a copy number in each read name.
SPILL_COPIES=400
for m in 1 2; do
  awk -v n=$SPILL_COPIES '{ l[NR] = $0 }
    END {
      for (c = 0; c < n; c++)
        for (i = 1; i < NR; i += 4) {
          split(l[i], h, " "); id = h[1]; mate = substr(id, length(id) - 1)
          print substr(id, 1, length(id) - 2) "#" c mate; print l[i + 1]; print "+"; print l[i + 3]
        }
    }' "$READS/sim_$m.fq" > spill_$m.fq
done

spill_classify() {
  # spill_classify DIR THREADS [kraken2 args]  classifies the spill reads in
  # DIR with sequence output.  A file size limit and a timeout stop runs that
  # write without bound.
  local dir=$1 threads=$2; shift 2
  rm -rf "$dir"; mkdir -p "$dir"
  ( cd "$dir" && ulimit -f 102400 &&
      with_timeout 120 "$KRAKEN2" --db ../db --paired --threads "$threads" \
        --report k.report --classified-out 'c#.fq' --unclassified-out 'u#.fq' \
        "$@" ../spill_1.fq ../spill_2.fq )
}

spill_same() {
  # spill_same DIR  compares the outputs in DIR with the single-threaded run
  # and fails if any spill file is left in DIR or in spill_tmp.
  local f left=()
  for f in k.out k.report c_1.fq c_2.fq u_1.fq u_2.fq; do
    same spill_p1/$f "$1/$f" || return 1
  done
  for f in "$1"/* "$1"/.[!.]* spill_tmp/* spill_tmp/.[!.]*; do
    [ -e "$f" ] || continue
    case ${f#"$1"/} in
      k.out|k.report|[cu]_[12].fq) ;;
      *) left+=("$f") ;;
    esac
  done
  [ ${#left[@]} -eq 0 ] || { echo "files left behind:"; printf '%s\n' "${left[@]}"; return 1; }
}

spill_baseline() {
  [ -s spill_p1/k.out ] || spill_classify spill_p1 1 --output k.out
}

# Which blocks spill depends on thread timing, and a run may spill none, so
# each check repeats the multi-threaded run.
SPILL_RUNS=3

spill_to_file() {
  spill_baseline || return 1
  for _ in $(seq $SPILL_RUNS); do
    spill_classify spill_p4 4 --output k.out || return 1
    spill_same spill_p4 || return 1
  done
  echo "Kraken output, report and sequence files identical at --threads 1 and 4 in $SPILL_RUNS runs"
}

spill_to_stdout() {
  # Kraken output on stdout puts spill files in TMPDIR.
  spill_baseline || return 1
  for _ in $(seq $SPILL_RUNS); do
    rm -rf spill_tmp; mkdir spill_tmp
    TMPDIR="$WORK/spill_tmp" spill_classify spill_p4_stdout 4 > spill_stdout.out || return 1
    mv spill_stdout.out spill_p4_stdout/k.out
    spill_same spill_p4_stdout || return 1
  done
  echo "Kraken output on stdout, report and sequence files identical at --threads 1 and 4 in $SPILL_RUNS runs"
}
mkdir -p spill_tmp
check classify.spill.output_file_matches_one_thread spill_to_file
check classify.spill.stdout_matches_one_thread spill_to_stdout
# The spill inputs and outputs take about 270 MB; the logs keep any differences.
rm -rf spill_1.fq spill_2.fq spill_p1 spill_p4 spill_p4_stdout spill_tmp

self_classify() {
  # Each reference sequence should come back at its own taxid or an ancestor.
  "$KRAKEN2" --db db --output self.out lib.fa || return 1
  awk -F '\t' '
    FILENAME == ARGV[1] { parent[$1] = $2; next }
    function on_path(a, b,   x) {  # is a equal to b or an ancestor of b?
      for (x = b; x != ""; x = parent[x]) { if (x == a) return 1; if (x == 1) break }
      return 0
    }
    { n++; split($2, f, "|"); t = f[2]; c = $3
      # A reference may be called at a descendant of its own taxid (SARS at
      # SARS-CoV-2) because Kraken 2 scores root-to-leaf paths, so accept
      # anything on the same lineage.
      if (c == 0 || !(on_path(c, t) || on_path(t, c))) { bad++; print "self-classification off: " $2 " called=" c } }
    END { print n " reference sequences, " bad + 0 " off-lineage"; exit (bad > 0 || n != 38) }' \
    "$PARENTS" self.out
}
check classify.self_classification self_classify

# ---------------------------------------------------------------------------
# Fast (non-deterministic) build: calls may shift slightly, accuracy may not.
fast_build() {
  build_db_kraken2_build dbfast "$BT4" --fast-build &&
    "$KRAKEN2" --db dbfast --paired --threads 2 --output fast.out "$R1" "$R2"
}
check build.fast_build fast_build
check build.fast_build.accuracy accuracy fast.out

# ---------------------------------------------------------------------------
# Negative cases: bad input must produce a clean error, not a crash or a
# silent success.

corrupt_db() {
  # corrupt_db NAME FILE BYTES  copies db to dbNAME with FILE cut to BYTES.
  rm -rf "db$1"; cp -R db "db$1"
  head -c "$3" "db/$2" > "db$1/$2"
}
corrupt_db _trunc_hash hash.k2d 4096
corrupt_db _trunc_taxo taxo.k2d 100
corrupt_db _short_opts opts.k2d 32
rm -rf db_long_opts; cp -R db db_long_opts; { cat db/opts.k2d; head -c 4096 /dev/zero; } > db_long_opts/opts.k2d

check negative.missing_db must_fail "$KRAKEN2" --db no_such_db "$R1"
check negative.missing_reads must_fail "$KRAKEN2" --db db --output /dev/null no_such_reads.fq
check negative.missing_mate must_fail "$KRAKEN2" --db db --paired --output /dev/null "$R1" no_such_reads.fq
check negative.k2_missing_reads must_fail "$PYTHON" "$K2" classify --db db --output /dev/null no_such_reads.fq
check negative.paired_odd_files must_fail "$KRAKEN2" --db db --paired "$R1"
check negative.confidence_out_of_range must_fail "$KRAKEN2" --db db --confidence 1.5 "$R1"
check negative.truncated_hash must_fail "$KRAKEN2" --db db_trunc_hash --output /dev/null "$R1"
check negative.truncated_taxo must_fail "$KRAKEN2" --db db_trunc_taxo --output /dev/null "$R1"
check negative.short_opts must_fail "$KRAKEN2" --db db_short_opts --output /dev/null "$R1"
rejects_opts() {
  # The error must blame opts.k2d, not whatever the overflow clobbered.
  local out
  out=$(must_fail "$@" 2>&1); local rc=$?
  printf '%s\n' "$out"
  [ $rc -eq 0 ] && printf '%s' "$out" | grep -q 'opts'
}
xfail_check negative.long_opts "an oversized opts.k2d overflows a stack struct in classify.cc; crashes on macOS, misleading error on Linux" \
  rejects_opts "$KRAKEN2" --db db_long_opts --output /dev/null "$R1"

names_no_trailer() {
  # names.dmp lines must end in a tab-pipe because the parser strips the last
  # two characters of every line.  Without it "scientific name" no longer
  # matches and every taxon silently loses its name.  The build should be
  # rejected or the names kept.
  new_db_dir dbnt
  perl -pi -e 's/\t\|$//' dbnt/taxonomy/names.dmp
  "$KBUILD" --db dbnt --add-to-library lib.fa > /dev/null 2>&1 || return 1
  if ! "$KBUILD" --db dbnt --build --threads "$BT" --subblock-size $SUBBLOCK; then
    echo "build rejected the file"
    return 0
  fi
  "$KINSPECT" --db dbnt | grep -v '^#' | head -n 3
  "$KINSPECT" --db dbnt | grep -q 'Viruses$'
}
xfail_check negative.names_dmp_without_trailer "names.dmp lines without a trailing tab-pipe give every taxon an empty name, with no error (taxonomy.cc)" \
  names_no_trailer

orphan_build() {
  # A library sequence whose taxid is missing from nodes.dmp.  The build
  # should fail or at least name the taxid; today it is dropped silently.
  rm -rf dborph; mkdir -p dborph/taxonomy; cp "$REPO"/data/*.dmp dborph/taxonomy/
  { cat lib.fa; printf '>kraken:taxid|999999|ORPHAN.1 orphan\n'; awk 'NR > 1 && NR < 40' "$REPO/data/Lambda.fa" | tr ACGT TGCA; } > orph.fa
  "$KBUILD" --db dborph --add-to-library orph.fa || return 1
  local rc=0
  "$KBUILD" --db dborph --build --threads "$BT" --subblock-size $SUBBLOCK > orph_build.log 2>&1 || rc=$?
  cat orph_build.log
  [ $rc -ne 0 ] || grep -q 999999 orph_build.log
}
xfail_check build.orphan_taxid_reported "taxids absent from nodes.dmp are dropped without a warning (taxonomy.cc MarkNode)" \
  orphan_build

# ---------------------------------------------------------------------------
# The k2 front end.

k2() { "$PYTHON" "$K2" "$@"; }

if [ "${E2E_SKIP_K2:-0}" = 1 ]; then
  skip k2 "E2E_SKIP_K2=1"
else
  k2_build() {
    # k2_build DIR [k2 build args] [-- FASTA...]  builds with the k2 front end
    local dir=$1; shift
    local args=() files=()
    while [ $# -gt 0 ] && [ "$1" != -- ]; do args+=("$1"); shift; done
    [ $# -gt 0 ] && shift
    files=("$@"); [ ${#files[@]} -eq 0 ] && files=(lib.fa)
    new_db_dir "$dir"
    local f
    for f in "${files[@]}"; do
      k2 add-to-library --db "$dir" --files "$f" || return 1
    done
    k2 build --db "$dir" --sub-block-size $SUBBLOCK "${args[@]}"
  }

  check k2.build k2_build k2db --threads "$BT"
  check k2.build.db_check "$PYTHON" "$TESTS/k2db_check.py" k2db
  check k2.classify.matches_kraken2 bash -c "
    '$PYTHON' '$K2' classify --db k2db --paired --output k2.out --report k2.report '$R1' '$R2' &&
    cmp paired.out k2.out && cmp paired.report k2.report"
  check k2.classify.interleaved.matches_paired bash -c "
    '$PYTHON' '$K2' classify --db k2db --interleaved --output k2i.out interleaved.fq && cmp paired.out k2i.out"
  check k2.classify.gzip.matches_plain bash -c "
    '$PYTHON' '$K2' classify --db k2db --paired --output k2gz.out r1.fq.gz r2.fq.gz && cmp paired.out k2gz.out"
  check k2.inspect.matches_kraken2_inspect bash -c "
    '$PYTHON' '$K2' inspect --db k2db | grep -v '^# Database: ' > k2inspect.txt && diff inspect.txt k2inspect.txt"

  # Two databases that split the references between them.  Multi-database
  # classification runs classify once per database and combines the results
  # with the merge binary.
  split_dbs() {
    cat "$REPO"/data/Flu*.fa "$REPO"/data/HIV*.fa > libA.fa
    cat "$REPO"/data/COVID_19.fa "$REPO"/data/Lambda.fa "$REPO"/data/MERS.fa "$REPO"/data/SARS.fa > libB.fa
    k2_build k2dbA --threads "$BT" -- libA.fa && k2_build k2dbB --threads "$BT" -- libB.fa
  }
  check k2.multi_db.build split_dbs

  multi_classify() {
    # multi_classify NAME [k2 classify args]
    local name=$1; shift
    k2 classify --db k2dbA,k2dbB --paired --output "$name.out" --report "$name.report" "$@" "$R1" "$R2"
  }
  multi_lca() {
    multi_classify multi --threads "$BT" &&
      golden multi.out multi.out && golden multi.report multi.report
  }
  check k2.multi_db.resolve_lca.golden multi_lca
  check k2.multi_db.resolve_lca.accuracy accuracy multi.out

  # merge's parallel path needs at least two threads, because one of them
  # only writes the blocks the others produce.  With a single thread, from a
  # build without OpenMP or from OMP_THREAD_LIMIT=1, merge must fall back to
  # its serial path rather than wait forever.
  threads_one_team() (
    # threads_one_team [NAME=VALUE]  runs in a subshell, with NAME exported
    if [ $# -gt 0 ]; then export "${1?}"; fi
    with_timeout 15 "$PYTHON" "$K2" classify --db k2dbA,k2dbB --paired --threads 2 \
      --output multi_t2.out --report multi_t2.report "$R1" "$R2" &&
      cmp multi.out multi_t2.out && cmp multi.report multi_t2.report
  )
  check k2.multi_db.order_independent bash -c "
    '$PYTHON' '$K2' classify --db k2dbB,k2dbA --paired --output multi_ba.out --report multi_ba.report '$R1' '$R2' &&
    cmp multi.out multi_ba.out && cmp multi.report multi_ba.report"

  # Quick mode prints no hit counts, so per-database quick calls cannot be
  # merged, and k2 must refuse them before classifying.
  multi_quick() {
    rm -f multi_quick.out
    local out
    out=$(must_fail k2 classify --db k2dbA,k2dbB --paired --quick --output multi_quick.out "$R1" "$R2")
    local rc=$?
    printf '%s\n' "$out"
    [ $rc -eq 0 ] && ! printf '%s\n' "$out" | grep -q Traceback &&
      printf '%s\n' "$out" | grep -q -- --quick && [ ! -s multi_quick.out ]
  }
  check k2.multi_db.quick_rejected multi_quick
  # merge labels reads itself with --use-names, so its labels must match
  # classify's, including "unclassified" for taxid 0.
  multi_names() {
    k2 classify --db k2db --paired --use-names --output names_single.out "$R1" "$R2" &&
      multi_classify multi_names --use-names &&
      cmp <(sed 's/[[:space:]]*$//' names_single.out) <(sed 's/[[:space:]]*$//' multi_names.out)
  }
  check k2.multi_db.use_names_matches_single_db multi_names
  # With only the first 20 pairs, some taxa get no reads, so the report
  # differs from the default one only by its zero-count rows.
  multi_zero() {
    head -n 80 "$R1" > zero_1.fq; head -n 80 "$R2" > zero_2.fq
    k2 classify --db k2db --paired --report-zero-counts --output zero_single.out \
      --report zero_single.report zero_1.fq zero_2.fq &&
    k2 classify --db k2dbA,k2dbB --paired --report-zero-counts --output multi_zero.out \
      --report multi_zero.report zero_1.fq zero_2.fq &&
    k2 classify --db k2db --paired --output /dev/null --report nozero.report zero_1.fq zero_2.fq || return 1
    echo "rows: $(wc -l < nozero.report) without --report-zero-counts, $(wc -l < zero_single.report) with"
    ! cmp -s nozero.report zero_single.report && same zero_single.report multi_zero.report
  }
  check k2.multi_db.report_zero_counts_matches_single_db multi_zero
  if [ $OPENMP = 0 ]; then
    check k2.multi_db.threads_without_openmp threads_one_team
  else
    check k2.multi_db.threads_with_thread_limit_1 threads_one_team OMP_THREAD_LIMIT=1
  fi
  # classify ends each hit list with a space and merge does not, so trailing
  # whitespace is not compared.
  check k2.multi_db.matches_single_db bash -c "
    cmp <(sed 's/[[:space:]]*\$//' paired.out) <(sed 's/[[:space:]]*\$//' multi.out) && cmp paired.report multi.report"
  multi_policies() {
    multi_classify first --merge-policy first-classified &&
      golden multi.first.out first.out &&
    multi_classify last --merge-policy last-classified &&
      golden multi.last.out last.out
  }
  check k2.multi_db.policies.golden multi_policies
  multi_mm() {
    multi_classify multimm --report-minimizer-data && golden multi.minimizer.report multimm.report
  }
  check k2.multi_db.minimizer_report.golden multi_mm
  check k2.multi_db.accuracy_first accuracy first.out

  multi_missing_map() {
    rm -rf k2dbB_nomap; cp -R k2dbB k2dbB_nomap; rm k2dbB_nomap/seqid2taxid.map
    must_fail k2 classify --db k2dbA,k2dbB_nomap --paired --output nomap.out "$R1" "$R2"
  }
  check k2.multi_db.missing_seqid2taxid_fails multi_missing_map

  # A classify or merge that fails must make k2 exit non-zero with a message
  # naming the binary: not exit 0, raise a Python traceback, or wait forever
  # on the FIFO that feeds compressed input to classify.  k2 must also remove
  # its temporary files.  A database with a truncated hash.k2d makes classify
  # fail, and a merge earlier in PATH that exits 3 makes merge fail.
  k2_fails_cleanly() {
    # k2_fails_cleanly [k2 classify args]
    rm -rf k2fail_tmp; mkdir k2fail_tmp
    local out rc left
    out=$(TMPDIR="$WORK/k2fail_tmp" with_timeout 30 "$PYTHON" "$K2" classify "$@" 2>&1)
    rc=$?
    printf '%s\n' "$out" | tail -n 20
    echo "exit status $rc"
    left=$(ls -A k2fail_tmp)
    [ -z "$left" ] || echo "temporary files left: $left"
    [ $rc -ne 0 ] && [ $rc -ne 124 ] && { [ $rc -lt 129 ] || [ $rc -gt 159 ]; } &&
      ! printf '%s\n' "$out" | grep -q Traceback &&
      printf '%s\n' "$out" | grep -q -E '(classify|merge) (exited with status|was terminated by)' &&
      [ -z "$left" ]
  }
  env_fails_cleanly() (
    # env_fails_cleanly NAME=VALUE [k2 classify args]  runs k2_fails_cleanly
    # in a subshell with one variable exported.
    export "${1?}"; shift
    k2_fails_cleanly "$@"
  )
  rm -rf k2db_trunc k2dbB_trunc k2fail_bin
  cp -R k2db k2db_trunc; cp -R k2dbB k2dbB_trunc; mkdir k2fail_bin
  head -c 4096 k2db/hash.k2d > k2db_trunc/hash.k2d
  head -c 4096 k2dbB/hash.k2d > k2dbB_trunc/hash.k2d
  printf '#!/bin/sh\nexit 3\n' > k2fail_bin/merge; chmod +x k2fail_bin/merge
  check k2.classify.failed_classify_exits_nonzero \
    k2_fails_cleanly --db k2db_trunc --paired --output fail.out "$R1" "$R2"
  check k2.classify.gzip.failed_classify_exits_nonzero \
    k2_fails_cleanly --db k2db_trunc --paired --output fail.out r1.fq.gz r2.fq.gz
  check k2.multi_db.failed_classify_exits_nonzero \
    k2_fails_cleanly --db k2dbA,k2dbB_trunc --paired --output fail.out "$R1" "$R2"
  check k2.multi_db.failed_merge_exits_nonzero \
    env_fails_cleanly PATH="$WORK/k2fail_bin:$PATH" --db k2dbA,k2dbB --paired --output fail.out "$R1" "$R2"

  # Three- and four-way splits.  The influenza A subtypes, the HIV types and
  # the betacoronaviruses are spread across databases, so many calls depend
  # on combining hits from several of them.  k2 merges the per-database
  # outputs pairwise, so these exercise merges of already merged output,
  # which two databases never reach.
  split_more() {
    local d="$REPO/data"
    cat "$d/FluA_H1N1.fa" "$d/FluA_H2N2.fa" "$d/HIV_1.fa" "$d/COVID_19.fa" > lib3a.fa
    cat "$d/FluA_H3N2.fa" "$d/HIV_2.fa" "$d/SARS.fa" "$d/Lambda.fa" > lib3b.fa
    cat "$d/FluB.fa" "$d/MERS.fa" > lib3c.fa
    cat "$d/FluA_H1N1.fa" "$d/HIV_1.fa" "$d/COVID_19.fa" > lib4a.fa
    cat "$d/FluA_H2N2.fa" "$d/HIV_2.fa" "$d/SARS.fa" > lib4b.fa
    cat "$d/FluA_H3N2.fa" "$d/Lambda.fa" > lib4c.fa
    cat "$d/FluB.fa" "$d/MERS.fa" > lib4d.fa
    local x
    for x in 3a 3b 3c 4a 4b 4c 4d; do
      k2_build k2db$x --threads "$BT" -- lib$x.fa || return 1
    done
  }
  check k2.multi_db.split_3_and_4.build split_more
  DB3=k2db3a,k2db3b,k2db3c
  DB4=k2db4a,k2db4b,k2db4c,k2db4d

  same_ignoring_trailing_space() {
    cmp <(sed 's/[[:space:]]*$//' "$1") <(sed 's/[[:space:]]*$//' "$2")
  }
  # matches_single DBS NAME SINGLE_PREFIX [k2 classify args]  classifies the
  # fixture reads with DBS and compares output and report with the
  # single-database files SINGLE_PREFIX.out and SINGLE_PREFIX.report.
  matches_single() {
    local dbs=$1 name=$2 single=$3; shift 3
    k2 classify --db "$dbs" --paired --output "$name.out" --report "$name.report" "$@" "$R1" "$R2" &&
      same_ignoring_trailing_space "$single.out" "$name.out" &&
      cmp "$single.report" "$name.report"
  }
  check k2.multi_db.three.matches_single_db matches_single "$DB3" multi3 paired
  check k2.multi_db.four.matches_single_db matches_single "$DB4" multi4 paired

  # k2 merge-classifications combines outputs that k2 classify wrote for each
  # database separately.  With three inputs it also merges merged output.
  merge_files_matches_single() {
    local d files=()
    for d in k2db3a k2db3b k2db3c; do
      k2 classify --db "$d" --paired --output "mc_$d.out" "$R1" "$R2" || return 1
      files+=("mc_$d.out")
    done
    k2 merge-classifications --db "$DB3" --output mc3.out --report mc3.report "${files[@]}" &&
      same_ignoring_trailing_space paired.out mc3.out && cmp paired.report mc3.report
  }
  check k2.merge_classifications.three.matches_single_db merge_files_matches_single
  check k2.multi_db.four.order_independent bash -c "
    '$PYTHON' '$K2' classify --db k2db4d,k2db4b,k2db4a,k2db4c --paired --output multi4_r.out --report multi4_r.report '$R1' '$R2' &&
    cmp multi4.out multi4_r.out && cmp multi4.report multi4_r.report"
  check k2.multi_db.three.confidence_matches_single_db matches_single "$DB3" multi3_c05 c05 --confidence 0.5

  # tests/data/mhg_*.fq hold read pairs whose calls at --minimum-hit-groups 5
  # and 10 differ between a single database and these splits when hit groups
  # are not combined correctly across more than two databases.
  mhg_matches_single() {
    # mhg_matches_single DBS NAME
    local dbs=$1 name=$2 m
    for m in 5 10; do
      "$KRAKEN2" --db db --paired --minimum-hit-groups $m --output "mhg_single.$m.out" \
        "$MR1" "$MR2" 2> /dev/null &&
      k2 classify --db "$dbs" --paired --minimum-hit-groups $m --output "$name.$m.out" "$MR1" "$MR2" &&
      same_ignoring_trailing_space "mhg_single.$m.out" "$name.$m.out" || return 1
    done
  }
  MR1="$READS/mhg_1.fq"
  MR2="$READS/mhg_2.fq"
  check k2.multi_db.three.min_hit_groups_matches_single_db \
    mhg_matches_single "$DB3" mhg3
  check k2.multi_db.four.min_hit_groups_matches_single_db \
    mhg_matches_single "$DB4" mhg4

  # With several databases, k2 writes --classified-out and --unclassified-out
  # itself from the read names merge reports as classified, so the sequence
  # files must equal those classify writes for the single database.
  seqout_matches_single() {
    # seqout_matches_single DBS NAME [k2 classify args and reads]
    local dbs=$1 name=$2; shift 2
    local a pat="" suffixes=("") s k
    for a in "$@"; do
      case $a in --paired|--interleaved) pat="#"; suffixes=(_1 _2) ;; esac
    done
    rm -f "$name".*
    k2 classify --db k2db --output "$name.single.out" \
      --classified-out "$name.single.cl$pat.seq" --unclassified-out "$name.single.un$pat.seq" "$@" &&
    k2 classify --db "$dbs" --output "$name.out" \
      --classified-out "$name.cl$pat.seq" --unclassified-out "$name.un$pat.seq" "$@" || return 1
    for s in "${suffixes[@]}"; do
      for k in cl un; do
        [ -s "$name.single.$k$s.seq" ] || { echo "$name.single.$k$s.seq is empty"; return 1; }
        same "$name.single.$k$s.seq" "$name.$k$s.seq" || return 1
      done
    done
    echo "$(grep -c '^C' "$name.out") classified and $(grep -c '^U' "$name.out") unclassified reads match the single database"
  }
  check k2.multi_db.seqout.paired_matches_single_db \
    seqout_matches_single k2dbA,k2dbB seqout_p --paired "$R1" "$R2"
  check k2.multi_db.seqout.gzip_matches_single_db \
    seqout_matches_single "$DB3" seqout_gz --paired r1.fq.gz r2.fq.gz
  check k2.multi_db.seqout.interleaved_matches_single_db \
    seqout_matches_single k2dbA,k2dbB seqout_i --interleaved interleaved.fq
  check k2.multi_db.seqout.single_fasta_matches_single_db \
    seqout_matches_single "$DB3" seqout_fa r1.fa
  check k2.multi_db.seqout.min_base_quality_matches_single_db \
    seqout_matches_single k2dbA,k2dbB seqout_mbq --paired --minimum-base-quality 20 "$R1" "$R2"

  k2_orphan() {
    rm -rf k2orph; new_db_dir k2orph
    k2 add-to-library --db k2orph --files orph.fa || return 1
    k2 build --db k2orph --threads "$BT" --sub-block-size $SUBBLOCK > k2orph.log 2>&1; local rc=$?
    cat k2orph.log
    [ $rc -ne 0 ] || grep -q 999999 k2orph.log
  }
  check k2.build.orphan_taxid_reported k2_orphan

  cell40() {
    # 40-bit cells should classify as well as 32-bit cells.  With fewer than
    # 129 taxa value_bits is below 8 and hashed_key drops key bits.
    k2_build k2db40 --threads "$BT" --cht-cell-size 40 || return 1
    k2 classify --db k2db40 --paired --output c40.out "$R1" "$R2" || return 1
    local h32 h40
    h32=$(hit_fraction paired.out); h40=$(hit_fraction c40.out)
    echo "k-mer hit fraction: 32-bit cells $h32, 40-bit cells $h40"
    "$PYTHON" "$TESTS/k2db_check.py" k2db40
    local rc_check=$?
    accuracy c40.out && [ $rc_check -eq 0 ] &&
      perl -e 'exit !($ARGV[1] >= 0.98 * $ARGV[0])' "$h32" "$h40"
  }
  xfail_check k2.build.cht_cell_40 "40-bit cells with value_bits < 8 lose key bits in 32-bit int shifts (compact_hash.h CompactHashCell40::hashed_key)" \
    cell40
fi

# ---------------------------------------------------------------------------
T1=$(now)
TOTAL_SECONDS=$(perl -e 'printf "%.1f", $ARGV[1] - $ARGV[0]' "$T0" "$T1")

say ""
say "Summary: $n_pass passed, $n_fail failed, $n_xfail expected failures, $n_xpass unexpected passes, $n_skip skipped (${TOTAL_SECONDS}s; kraken2-build ${BUILD_SECONDS}s)"

if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then
  {
    echo "### e2e: ${E2E_LABEL:-$(uname -s)}"
    echo ""
    echo "$n_pass passed, $n_fail failed, $n_xfail expected failures (known bugs), $n_xpass unexpected passes, $n_skip skipped in ${TOTAL_SECONDS}s"
    echo ""
    if [ $((n_fail + n_xpass + n_xfail)) -gt 0 ]; then
      echo "| Status | Check | Detail |"
      echo "|---|---|---|"
      awk -F '\t' '$1 != "PASS" { printf "| %s | `%s` | %s |\n", $1, $2, $3 }' "$RESULTS" |
        sort -t '|' -k2,2 -s | sed 's/| FAIL |/| **FAIL** |/; s/| XPASS |/| **XPASS** |/'
      echo ""
    fi
  } >> "$GITHUB_STEP_SUMMARY"
fi

if [ $n_fail -gt 0 ]; then
  exit 1
fi
if [ $n_xpass -gt 0 ] && [ "${E2E_XPASS_OK:-0}" != 1 ]; then
  say "Unexpected passes: a known bug appears fixed. Turn those xfail_check calls into check calls."
  exit 1
fi
exit 0
