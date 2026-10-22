#!/usr/bin/env python3
"""Build and run the isolated ResolveTree correctness/performance benchmark."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rounds", type=int, default=5, help="timed samples; report median")
    parser.add_argument("--check-only", action="store_true", help="skip timing, check equivalence")
    parser.add_argument("--sanitize", action="store_true", help="enable ASan/UBSan; skip timing")
    parser.add_argument("--output", type=Path, help="also save benchmark CSV to this path")
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error("--rounds must be positive")
    root = Path(__file__).resolve().parents[2]
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    with tempfile.TemporaryDirectory(prefix="k2-resolve-tree-") as temporary:
        binary = Path(temporary) / "benchmark"
        flags = ["-std=c++17", "-O3", "-Wall", "-Wno-unknown-pragmas", "-DLINEAR_PROBING",
                 "-ffunction-sections", "-fdata-sections", "-I", str(root / "src")]
        if args.sanitize:
            flags += ["-O1", "-g", "-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
        # omp_hack supplies serial fallbacks; benchmark only measures tree decisions.
        sources = [root / "tests/benchmarks/resolve_tree_benchmark.cc",
                   root / "src/taxonomy.cc", root / "src/mmap_file.cc",
                   root / "src/omp_hack.cc"]
        link_gc = "-Wl,-dead_strip" if os.uname().sysname == "Darwin" else "-Wl,--gc-sections"
        subprocess.run(compiler + flags + [str(p) for p in sources]
                       + [link_gc, "-o", str(binary)], check=True)
        result = subprocess.run([str(binary), temporary, str(args.rounds),
                                 str(int(args.check_only or args.sanitize))],
                                check=True, text=True, stdout=subprocess.PIPE)
        print(result.stdout, end="")
        if args.output:
            args.output.write_text(result.stdout)


if __name__ == "__main__":
    main()
