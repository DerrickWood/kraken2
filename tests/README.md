# Kraken 2 tests

Everything here runs offline in under a minute. It uses the small viral
reference set in `data/`, which has 38 sequences, 40 taxa and a 420 KB hash table.

| File | Purpose |
|---|---|
| `e2e.sh INSTALL_DIR [WORK_DIR]` | End-to-end suite against an install made by `install_kraken2.sh` |
| `update_golden.sh INSTALL_DIR` | Regenerate `golden/` after an intended output change |
| `pr_compare.py BASE HEAD WORK [BIG]` | Compare two installs read by read (used by the PR job) |
| `k2db_check.py DB` | Header-only sanity check of a built database (load, bit budgets, orphans) |
| `data/sim_1.fq`, `data/sim_2.fq` | Frozen reads: 400 simulated pairs, 20 chimeric tie pairs, 20 random pairs |
| `data/make_reads.sh` | How those reads were made (rerun only on purpose) |
| `data/mhg_1.fq`, `data/mhg_2.fq` | 18 read pairs whose calls at `--minimum-hit-groups` 5 and 10 depend on combining hit groups across three or more databases |
| `golden/` | Expected `--output`, `--report` and `kraken2-inspect` files |

## Running locally

```sh
./install_kraken2.sh /tmp/k2            # on macOS, see "macOS" below
tests/e2e.sh /tmp/k2                    # prints PASS/FAIL/XFAIL per check
```

The last line summarizes the run. The script exits non-zero if any check FAILs or XPASSes. The
work directory, which is printed at the top and kept after the run, has `logs/<check>.log`
for every check, plus `actual/` with each golden-compared output and a `.diff` next to any
that differ.

`E2E_ONLY=regex` runs a subset, for example `E2E_ONLY='^k2\.' tests/e2e.sh /tmp/k2`.

macOS with Homebrew libomp:

```sh
P=$(brew --prefix libomp)
CXX=clang++ KRAKEN2_SKIP_FOPENMP="-Xpreprocessor -fopenmp -I$P/include" \
  LDFLAGS="-L$P/lib -lomp" ./install_kraken2.sh /tmp/k2
```

To build without OpenMP, set `KRAKEN2_SKIP_FOPENMP=` (empty). The suite detects this and builds with one thread.

## What the suite checks

- **Build.** `kraken2-build` builds the fixture database, which must pass `k2db_check.py`. The
  database files must be byte-identical at 1, 2 and 4 threads. A `--fast-build` database must
  still meet the accuracy thresholds.
- **Classification.** Paired FASTQ, single FASTQ, FASTA, gzip (auto-detected and
  `--gzip-compressed`) and interleaved input. Each mode must produce the same calls as the
  others. `--classified-out` and `--unclassified-out` must match the calls.
  `--memory-mapping` must not change output, and neither may 1, 2, 3, 4 or 8 threads.
  On a larger input built from the fixture reads (176,000 pairs, 5 input blocks), 4 threads
  with `--classified-out` and `--unclassified-out` must give the same Kraken output, report
  and sequence files as 1 thread and leave no spill files, both with `--output` and with
  Kraken output on stdout (spill files in `TMPDIR`).
- **Accuracy.** Each read's call is checked against the true taxid in its name. At least 97%
  of simulated pairs must get the exact taxid and 99% must land on the true lineage. At
  most 2 of the 20 random pairs may be classified.
- **Tie-breaking.** In each chimeric pair, one mate comes from each of two taxa, and every
  k-mer maps to that mate's own taxon. The pair must be called at the LCA of the two taxa.
  These pairs exercise the tie branch of `ResolveTree`, which simulated reads almost never
  reach.
- **Golden files.** These cover `--output`, `--report`, `--report-minimizer-data`,
  `--use-mpa-style`, `--report-zero-counts`, `--use-names`, `--confidence 0.5`,
  `--minimum-hit-groups 5`, `--quick`, `--minimum-base-quality 20` and `kraken2-inspect`.
  Each reference sequence must also classify to its own lineage.
- **k2.** `k2 add-to-library`, `k2 build`, `k2 classify` and `k2 inspect` must agree with the
  Perl front end. Multi-database classification uses two databases that split the
  references. It is checked against golden files for each `--merge-policy`, against the
  single-database result, and for order independence. Three- and four-way splits, which
  spread the influenza A subtypes, the HIV types and the betacoronaviruses across databases,
  must also match the single-database output and report, including with `--confidence 0.5`
  and, on `data/mhg_*.fq`, with `--minimum-hit-groups` 5 and 10.
- **Negative cases.** Each of these must fail with a clean non-zero exit, not a crash:
  - a missing database;
  - a missing read file or mate file;
  - an odd number of `--paired` files;
  - a bad `--confidence`;
  - a truncated or oversized database file;
  - a `names.dmp` without the trailing delimiter;
  - a taxid that is absent from `nodes.dmp`.

## Reading CI results

The workflow is `.github/workflows/ci.yml`. On a pull request, open the run's **Summary** page.

1. **Build matrix.** The same build-and-test runs in five jobs:
   - ubuntu gcc;
   - ubuntu clang with libomp;
   - ubuntu without OpenMP;
   - macOS clang with libomp;
   - ubuntu clang with ASan and UBSan.

   Each job writes a short section with the compiler warning count and a table of every check
   that did not pass. If a job is red, look in that table for **FAIL** rows. The job log shows
   the first and last lines of the failing check's log. For a golden mismatch, the log lists
   the reads whose call changed. Download the `e2e-<job>` artifact for the full logs and diffs.
2. **XFAIL rows are expected.** Each one reproduces a known bug and names the source location
   (see the table below). They do not turn the job red.
3. **XPASS turns a job red on purpose.** It means a known bug no longer reproduces, which is
   usually because the PR fixed it. Change that `xfail_check` in `e2e.sh` to `check`. The fix
   then stays covered.
4. **PR vs base comparison.** This job builds both the base commit and the PR (merged into the
   base). Each builds its own fixture database and classifies the same reads. The summary shows:
   - whether the database files are identical;
   - every read whose call changed, with its true taxid;
   - report rows whose counts changed;
   - accuracy for both builds;
   - wall time and peak memory for build and classify, including a 500,000-pair timing set.

   This job never fails on a difference. It shows what the PR does to results, so a reviewer
   can decide whether the change was intended.
5. **Script checks.** These cover `perl -c`, `bash -n`, `py_compile` and a `k2 --help` smoke
   test. ShellCheck errors fail the job. ShellCheck warnings in the shipped scripts only appear
   in the summary.
6. **CMake build.** This job fails if CMake silently drops OpenMP or a core target is missing.
   It warns, without failing, that there is no `merge` target yet.

A PR that intentionally changes output should include regenerated golden files
(`tests/update_golden.sh /path/to/install`). Its PR comparison should show the same changes. Golden files that belong to expected-failure checks hold the output expected once the bug is fixed, so `update_golden.sh` leaves them alone unless `E2E_UPDATE_XFAIL_GOLDEN=1` is set; set it only when running against a build in which that bug is fixed.

## Known bugs marked as expected failures

| Check | Bug |
|---|---|
| `inspect.redirect_keeps_header` | On Linux, `kraken2-inspect --db DB > file` loses the `#` header lines, because `dump_table` reopens `/dev/fd/1` with truncation (`dump_table.cc`, `reports.cc`). Checked as a normal test on macOS. |
| `negative.long_opts` | An `opts.k2d` larger than 64 bytes is read into a stack struct and crashes `classify` (`classify.cc`) |
| `negative.names_dmp_without_trailer` | A `names.dmp` without a trailing tab-pipe gives every taxon an empty name, with no error (`taxonomy.cc`) |
| `build.orphan_taxid_reported` | `kraken2-build` silently drops sequences whose taxid is absent from `nodes.dmp` (`taxonomy.cc` `MarkNode`). `k2 build` reports them. |
| `k2.build.cht_cell_40` | 40-bit cells with fewer than 129 taxa lose key bits in `int` shifts (`compact_hash.h` `CompactHashCell40::hashed_key`) |
| `k2.multi_db.threads_without_openmp` | In a build without OpenMP, `merge -T 2` busy-waits forever (`merge.cc` `merge_classification_output_parallel`). Checked only in that build. |

## How `data/mhg_*.fq` was made

The pairs were selected from reads simulated with `data/simulator.pl` (`PERL_HASH_SEED=0
PERL_PERTURB_KEYS=0`, 60 bp reads: `--random-seed 7 --num-frags 4000`, names prefixed `a.`,
and `--random-seed 11 --num-frags 20000`, prefixed `b.`). Each kept pair gets a different call
at `--minimum-hit-groups` 5 or 10 from `k2 classify` with the three- or four-way split than
from the single database, on master at e7733dd. With hit groups combined correctly, every
pair gets the same call as the single database. The two-way split shows no difference on any
of these reads.
