#!/usr/bin/env python3
"""Compare two Kraken 2 installs (a PR and its base) on the test fixture.

Usage: tests/pr_compare.py BASE_INSTALL HEAD_INSTALL WORK_DIR [BIG_READS_PREFIX]

Each install builds its own database from data/ with kraken2-build and then
classifies the frozen fixture reads in tests/data (and, when given, a larger
simulated read set BIG_READS_PREFIX_1.fq / _2.fq used for timing).  The script
prints a Markdown summary: whether the database files are identical, reads
whose call changed, report rows whose counts changed, accuracy against the
truth in read names, and wall time and peak memory of build and classify.

It never fails on a difference, since a difference may be the point of the
PR; it exits non-zero only if a command cannot run at all.
"""
import os
import resource
import shutil
import struct
import subprocess
import sys
import time

TESTS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(TESTS)
REPS = 3  # timing repetitions; the minimum wall time is reported


def run(argv, log, cwd=None):
    """Run argv, append its output to log, return (wall seconds, peak RSS MiB).

    Peak RSS comes from getrusage(RUSAGE_CHILDREN), which reports the largest
    child seen so far, so each measurement runs in its own Python child
    process to keep the numbers independent.
    """
    code = (
        "import resource, subprocess, sys, time\n"
        "t = time.perf_counter()\n"
        "rc = subprocess.call(sys.argv[2:], stdout=open(sys.argv[1], 'ab'), stderr=subprocess.STDOUT)\n"
        "w = time.perf_counter() - t\n"
        "m = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss\n"
        "print(rc, w, m)\n"
    )
    out = subprocess.run([sys.executable, "-c", code, log] + argv, cwd=cwd,
                         check=True, capture_output=True, text=True).stdout.split()
    rc, wall, maxrss = int(out[0]), float(out[1]), int(out[2])
    if rc != 0:
        raise RuntimeError(f"{' '.join(argv)} exited with {rc}; see {log}")
    # ru_maxrss is KiB on Linux and bytes on macOS.
    mib = maxrss / (1024 * 1024) if sys.platform == "darwin" else maxrss / 1024
    return wall, mib


def best_of(argv, log, cwd=None, reps=REPS):
    results = [run(argv, log, cwd) for _ in range(reps)]
    return min(r[0] for r in results), max(r[1] for r in results)


def build(inst, db, log):
    if os.path.exists(db):
        shutil.rmtree(db)
    os.makedirs(os.path.join(db, "taxonomy"))
    for f in ("nodes.dmp", "names.dmp"):
        shutil.copy(os.path.join(REPO, "data", f), os.path.join(db, "taxonomy"))
    lib = os.path.join(os.path.dirname(os.path.dirname(db)), "lib.fa")
    run([os.path.join(inst, "kraken2-build"), "--db", db, "--add-to-library", lib], log)
    return run([os.path.join(inst, "kraken2-build"), "--db", db, "--build", "--threads", "2"], log)


def read_output(path):
    calls = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            calls[fields[1]] = fields[2]
    return calls


def read_report(path):
    rows = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            # pct, clade reads, taxon reads, [minimizers, distinct,] rank, taxid, name
            taxid = fields[-2]
            rows[taxid] = (int(fields[1]), int(fields[2]), fields[-1].strip())
    return rows


def truth_table(fastq):
    truth = {}
    with open(fastq) as f:
        for i, line in enumerate(f):
            if i % 4:
                continue
            name, _, rest = line[1:].rstrip("\n").partition(" ")
            name = name.rsplit("/", 1)[0]
            taxid = "0"
            if rest.startswith("kraken:taxid|"):
                taxid = rest.split("|")[1]
            truth[name] = taxid
    return truth


def db_check(db):
    p = subprocess.run([sys.executable, os.path.join(TESTS, "k2db_check.py"), db],
                       capture_output=True, text=True)
    return p.stdout.strip()


def pct(a, b):
    if a == 0:
        return ""
    return f"{(b - a) / a * 100:+.0f}%"


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    base, head, work = (os.path.abspath(p) for p in sys.argv[1:4])
    big = os.path.abspath(sys.argv[4]) if len(sys.argv) == 5 else None
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(work, "lib.fa"), "wb") as out:
        for name in sorted(os.listdir(os.path.join(REPO, "data"))):
            if name.endswith(".fa"):
                with open(os.path.join(REPO, "data", name), "rb") as f:
                    shutil.copyfileobj(f, out)
    r1 = os.path.join(TESTS, "data", "sim_1.fq")
    r2 = os.path.join(TESTS, "data", "sim_2.fq")
    truth = truth_table(r1)

    res = {}
    for label, inst in (("base", base), ("head", head)):
        d = os.path.join(work, label)
        os.makedirs(d, exist_ok=True)
        log = os.path.join(d, "log.txt")
        db = os.path.join(d, "db")
        r = res[label] = {}
        r["build"] = build(inst, db, log)
        r["check"] = db_check(db)
        k2 = os.path.join(inst, "kraken2")
        out = os.path.join(d, "fixture.out")
        rep = os.path.join(d, "fixture.report")
        r["classify"] = best_of([k2, "--db", db, "--paired", "--output", out,
                                 "--report", rep, "--report-minimizer-data", r1, r2], log)
        r["calls"] = read_output(out)
        r["report"] = read_report(rep)
        if big:
            bout = os.path.join(d, "big.out")
            r["big"] = best_of([k2, "--db", db, "--paired", "--threads", "2", "--output",
                                bout, big + "_1.fq", big + "_2.fq"], log)
            r["big_calls"] = read_output(bout)
        r["db_files"] = {f: open(os.path.join(db, f), "rb").read()
                         for f in ("hash.k2d", "taxo.k2d")}
        # opts.k2d is a raw struct with uninitialized padding, so compare fields.
        r["db_files"]["opts.k2d"] = struct.unpack(
            "<QQQQ?7xQiii4x", open(os.path.join(db, "opts.k2d"), "rb").read())

    b, h = res["base"], res["head"]
    md = []
    md.append("## PR vs base on the test fixture")
    md.append("")
    md.append("Both commits built a database from `data/` with `kraken2-build` and classified "
              f"the {len(truth)} fixture read pairs from `tests/data`. "
              "A difference is not necessarily a problem; it shows what the PR changes.")
    md.append("")

    same_db = [f for f in b["db_files"] if b["db_files"][f] == h["db_files"][f]]
    diff_db = [f for f in b["db_files"] if f not in same_db]
    if diff_db:
        md.append(f"**Database files differ:** {', '.join(diff_db)} "
                  f"(identical: {', '.join(same_db) or 'none'})")
    else:
        md.append("**Database files:** hash.k2d and taxo.k2d are byte-identical and opts.k2d fields match.")
    if b["check"] != h["check"]:
        md.append("")
        md.append("k2db_check.py output:")
        md.append("```")
        md.append("base: " + b["check"].replace("\n", "\n      "))
        md.append("head: " + h["check"].replace("\n", "\n      "))
        md.append("```")
    md.append("")

    def accuracy(calls):
        sim = [r for r in truth if truth[r] != "0"]
        exact = sum(calls.get(r) == truth[r] for r in sim)
        rand_cl = sum(calls.get(r, "0") != "0" for r in truth if truth[r] == "0")
        return exact, len(sim), rand_cl

    changed = sorted((r for r in truth if b["calls"].get(r) != h["calls"].get(r)),
                     key=lambda r: (r.split(".")[0], int(r.split(".")[1])))
    be, n, brc = accuracy(b["calls"])
    he, _, hrc = accuracy(h["calls"])
    md.append(f"**Changed calls (fixture):** {len(changed)} of {len(truth)} reads. "
              f"Exact-taxid accuracy {be}/{n} (base) vs {he}/{n} (head); "
              f"random reads classified {brc} vs {hrc}.")
    if changed:
        md.append("")
        md.append("| Read | True taxid | Base call | Head call |")
        md.append("|---|---|---|---|")
        for r in changed[:25]:
            md.append(f"| {r} | {truth[r]} | {b['calls'].get(r)} | {h['calls'].get(r)} |")
        if len(changed) > 25:
            md.append(f"| ... | {len(changed) - 25} more | | |")
    md.append("")

    if big:
        nb = len(b["big_calls"])
        bchanged = sum(b["big_calls"].get(r) != c for r, c in h["big_calls"].items())
        md.append(f"**Changed calls (timing set, {nb} pairs):** {bchanged}")
        md.append("")

    taxids = sorted(set(b["report"]) | set(h["report"]), key=lambda t: int(t))
    rows = []
    for t in taxids:
        bo, ho = b["report"].get(t), h["report"].get(t)
        if bo is None or ho is None or bo[:2] != ho[:2]:
            name = (bo or ho)[2]
            rows.append(f"| {t} | {name} | {bo[0] if bo else '-'} | {ho[0] if ho else '-'} "
                        f"| {bo[1] if bo else '-'} | {ho[1] if ho else '-'} |")
    if rows:
        md.append("**Report rows that changed:**")
        md.append("")
        md.append("| Taxid | Name | Clade reads (base) | Clade reads (head) | Taxon reads (base) | Taxon reads (head) |")
        md.append("|---|---|---|---|---|---|")
        md.extend(rows)
    else:
        md.append("**Report:** clade and taxon read counts identical for every taxon.")
    md.append("")

    md.append("**Resources** (wall time is the best of "
              f"{REPS} runs for classify; peak RSS of the largest child process)")
    md.append("")
    md.append("| Step | Base wall (s) | Head wall (s) | Change | Base peak RSS (MiB) | Head peak RSS (MiB) | Change |")
    md.append("|---|---|---|---|---|---|---|")
    steps = [("kraken2-build --build", "build"), ("classify fixture (1 thread)", "classify")]
    if big:
        steps.append((f"classify {len(b['big_calls'])} pairs (2 threads)", "big"))
    for title, key in steps:
        (bw, bm), (hw, hm) = b[key], h[key]
        md.append(f"| {title} | {bw:.2f} | {hw:.2f} | {pct(bw, hw)} | {bm:.0f} | {hm:.0f} | {pct(bm, hm)} |")
    md.append("")
    md.append("Timings on shared CI runners vary by 10-20% between runs; treat small changes as noise.")
    print("\n".join(md))


if __name__ == "__main__":
    main()
