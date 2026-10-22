# ResolveTree comparison

Run from the repository root:

```sh
python3 tests/benchmarks/run_resolve_tree.py --rounds 7 --output /tmp/resolve-tree.csv
python3 tests/benchmarks/run_resolve_tree.py --check-only
python3 tests/benchmarks/run_resolve_tree.py --sanitize
```

`CXX` selects the compiler (for example `CXX=/usr/bin/g++`). The runner compiles
with C++17 and `-O3` into a temporary directory and removes that directory on
exit. The sanitizer option requires installed ASan/UBSan runtimes, uses `-O1`,
and only checks correctness. No database downloads or extra Python packages
are required. The benchmark uses the serial OpenMP compatibility functions;
OpenMP is not involved in either tree resolver.

`resolve_tree_benchmark.cc` includes the actual `src/classify.cc` with its
program entry point renamed. Dead-section elimination removes unrelated
classifier functions at link time. `resolve_tree_reference.h` is a frozen copy
of the original resolver: do not update it to match production optimizations.
Timing and reference code stay outside `src/`.

Timing uses the standard `std::chrono::steady_clock` and
`std::chrono::duration<double, std::nano>`. The clock surrounds a batch of
function calls, not each individual call. Input generation, tree loading,
copying maps, correctness checks and printing are outside the measured region.
Both implementations receive identical inputs, are warmed up, and independently
calibrate repetitions to at least about 5 ms per sample (up to a repetition
cap). Samples alternate reference/current execution order. CSV columns report
the median nanoseconds per invocation and reference/current speedup.
Checksums and indirect calls prevent elimination of the measured computation.
Warmup also stabilizes the original function's insertion of a missing candidate
into `hit_counts`; performance measurements describe repeated, warmed inputs.

Fixtures cover balanced, chain, star and irregular parent trees. Each timed
workload cycles through 16 generated hit maps, with 1/2/4/8/16/64/256 distinct taxa
and confidence 0/0.5/0.9. The chain has 256 real nodes; the other trees have
4096. A fixed random seed makes generation reproducible. Correctness checks
compare both the returned taxid and the final map, including:

- empty input, zero sentinel, zero votes and a single hit;
- ancestor votes, ties and confidence fallback or failure;
- randomized insertion orders and 64-bit counts exercising 32-bit score wraparound;
- all 2187 assignments of absent/one/two votes over seven nodes, across five
  confidence thresholds, with and without an unmatched position.

The complete timing run checks 32,934 cases. `--check-only` checks 31,590 cases
because the 1344 timed-workload input checks are omitted.

## Production algorithm

Let D be the number of distinct hit taxa and H the tree height. The original
candidate scoring compares every pair of hits using ancestor walks,
O(D²H); confidence fallback can add O(DH²).

The new scoring walks each candidate's ancestors and looks up direct votes,
stopping below the smallest hit ID: expected O(DH), with no scratch storage.
It keeps the original candidate iteration order, tie LCA and uint32_t score
arithmetic. A singleton is resolved directly in constant time. With exactly two
hits, scoring uses an ancestor test rather than hashing every intervening
parent, and confidence uses at most two intersections in a stack array. This
avoids allocation overhead for small inputs without changing the O(DH) bound.

When confidence requires clade support, only the candidate's ancestor chain is
stored. Each hit contributes at LCA(candidate, hit), its first intersection with
that chain. Prefix sums then give clade support for successive ancestors. This
takes O(DH + D log H + H) time and O(H) temporary space; the D log H term comes
from binary search along decreasing internal IDs. For nonempty input, the
overall expected bound is O(DH). No confidence scratch allocation occurs when
direct support already suffices, including the default confidence of zero.

Including the existing hit map, space is O(D + H) during fallback and O(D)
otherwise. The original fallback used O(1) scratch space, so the speedup has a
small per-call memory tradeoff. There is no global tree index or persistent
thread-local cache.

These are function microbenchmarks, not end-to-end classification speedups.
Tree shape, hit distribution, CPU load and clock frequency affect timings;
the long-chain case is intentionally a stress case.

## Recorded run (2026-10-22)

[resolve_tree_results.csv](resolve_tree_results.csv) contains a seven-sample run
on an AMD Ryzen 7 8845H, Linux x86-64, GCC 16.2.1, C++17, `-O3`, without OpenMP.
Both functions run in the same process. These are synthetic warmed workloads,
with no CPU affinity or frequency pinning; small differences may be noise.

Representative balanced-tree results (microseconds per call):

| Distinct taxa | Confidence | Original (µs) | Optimized (µs) | Speedup |
|---:|---:|---:|---:|---:|
| 1 | 0 | 0.0062 | 0.0030 | 2.08x |
| 2 | 0.9 | 0.0456 | 0.0241 | 1.89x |
| 16 | 0 | 0.6914 | 0.1413 | 4.89x |
| 64 | 0.9 | 12.1593 | 1.0989 | 11.06x |
| 256 | 0 | 198.9476 | 3.7544 | 52.99x |
| 256 | 0.9 | 210.9462 | 6.0041 | 35.13x |

Not every workload improves materially: the chain with two hits and confidence
zero was 116.686 ns versus 116.553 ns (1.001x), effectively
unchanged. Report all rows rather than only the largest speedups. Before
optimizing, the same harness measured the original production function against
the frozen reference, giving approximately 1x.

Validation of this version:

- All 32,934 differential cases passed in the timing run.
- All 31,590 check-only cases passed with `_GLIBCXX_ASSERTIONS` and UBSan trap
  mode (`CXX='/usr/bin/g++ -D_GLIBCXX_ASSERTIONS -fsanitize=undefined
  -fsanitize-undefined-trap-on-error'`). ASan could not link because its runtime
  is absent on this host; no ASan result is claimed.
- Full `tests/e2e.sh` after rebuilding: 72 passed, 0 failed, 5 existing expected
  failures, 0 unexpected passes. The suite needs local multiprocessing sockets,
  so it was run outside the restricted sandbox.
