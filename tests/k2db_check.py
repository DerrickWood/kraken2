#!/usr/bin/env python3
"""Headroom checker for a built Kraken 2 database directory.

Reads only fixed-size headers plus the taxonomy node array, so it runs in
seconds on multi-hundred-GB databases (it never loads hash.k2d's table).
Layouts follow src/compact_hash.h (WriteTable), src/taxonomy.cc (WriteToDisk)
and src/kraken2_data.h (IndexOptions), assuming an LP64 little-endian build.
"""
import math, os, struct, sys

WARN, FAIL = [], []


def check(cond, level, msg):
    if not cond:
        (FAIL if level == "fail" else WARN).append(msg)


def read_hash_header(path):
    with open(path, "rb") as f:
        capacity, size, key_bits, value_bits = struct.unpack("<4Q", f.read(32))
    cell_bits = key_bits + value_bits
    check(cell_bits in (32, 40), "fail", f"unknown cell width {cell_bits}")
    expected = 32 + capacity * (cell_bits // 8)
    check(os.path.getsize(path) == expected, "fail",
          f"hash.k2d size {os.path.getsize(path)} != 32 + capacity*{cell_bits // 8} = {expected}")
    return dict(capacity=capacity, size=size, key_bits=key_bits,
                value_bits=value_bits, cell_bits=cell_bits)


def read_taxonomy(path):
    with open(path, "rb") as f:
        magic = f.read(8)
        check(magic == b"K2TAXDAT", "fail", f"bad taxo.k2d magic {magic!r}")
        node_count, name_len, rank_len = struct.unpack("<3Q", f.read(24))
        raw = f.read(56 * node_count)
    check(os.path.getsize(path) == 32 + 56 * node_count + name_len + rank_len,
          "fail", "taxo.k2d size inconsistent with header")
    max_ext, orphans, order_violations = 0, 0, 0
    for i in range(1, node_count):
        parent, first_child, child_ct, name_off, rank_off, ext, _ = \
            struct.unpack_from("<7Q", raw, 56 * i)
        max_ext = max(max_ext, ext)
        if ext == 0:
            orphans += 1          # slot allocated but never reached by BFS
        if i > 1 and parent >= i:
            order_violations += 1  # LCA assumes parent id < child id
    return dict(node_count=node_count, max_external=max_ext,
                orphans=orphans, order_violations=order_violations)


def read_opts(path):
    b = open(path, "rb").read()
    check(len(b) == 64, "fail",
          f"opts.k2d is {len(b)} bytes; classify reads st_size bytes into a 64-byte struct")
    b = b.ljust(64, b"\0")[:64]
    k, l, seed, toggle, dna, min_hash, revcom, db_ver, db_type = \
        struct.unpack("<QQQQ?7xQiii4x", b)
    return dict(k=k, l=l, dna=dna, min_hash=min_hash, revcom=revcom)


def unsuccessful_probe_fp(load, key_bits):
    # Knuth: expected cells examined in an unsuccessful linear-probing search
    # is 0.5*(1 + 1/(1-a)^2); all but the final empty one are key compares.
    compares = 0.5 * (1 + 1 / (1 - load) ** 2) - 1
    return compares * 2.0 ** -key_bits


def main(db):
    h = read_hash_header(os.path.join(db, "hash.k2d"))
    t = read_taxonomy(os.path.join(db, "taxo.k2d"))
    o = read_opts(os.path.join(db, "opts.k2d"))
    load = h["size"] / h["capacity"]
    need_vb = max(1, math.ceil(math.log2(t["node_count"])))
    vb_headroom = (1 << h["value_bits"]) - t["node_count"]
    fp = unsuccessful_probe_fp(min(load, 0.999), h["key_bits"])

    check(load < 0.85, "warn", f"occupancy {load:.3f}; probe chains grow as 1/(1-a)^2")
    check(h["value_bits"] >= need_vb, "fail", "value_bits too small for node_count")
    check(t["node_count"] <= 0.8 * (1 << h["value_bits"]), "warn",
          f"node_count {t['node_count']} within 20% of 2^{h['value_bits']}; "
          f"next growth costs a key bit")
    check(h["key_bits"] >= 12, "warn", f"only {h['key_bits']} key bits")
    check(fp < 1e-3, "warn", f"est. false key match per absent minimizer {fp:.2e}")
    check(not (h["cell_bits"] == 40 and h["value_bits"] < 8), "fail",
          "40-bit cells with value_bits < 8 lose key bits (compact_hash.h:62)")
    check(t["max_external"] < 2 ** 31, "warn",
          "external taxid >= 2^31 (build_db.h:155 int, reports.cc:176 uint32_t)")
    check(t["orphans"] == 0, "warn",
          f"{t['orphans']} taxonomy slots unreachable: seqid2taxid.map names taxids absent from nodes.dmp")
    check(t["order_violations"] == 0, "fail", "parent id >= child id; LCA would loop")
    check(o["l"] <= (31 if o["dna"] else 15), "fail", "l exceeds scanner limit")
    check(o["revcom"] == 1, "warn", "built with outdated revcom version")

    print(f"capacity={h['capacity']} size={h['size']} load={load:.3f} "
          f"cell={h['cell_bits']} key_bits={h['key_bits']} value_bits={h['value_bits']}")
    print(f"nodes={t['node_count']} (needs {need_vb} bits, headroom {vb_headroom}) "
          f"max_external_taxid={t['max_external']} k={o['k']} l={o['l']} "
          f"subsample_keep={1 - o['min_hash'] / 2**64:.3f} est_fp={fp:.2e}")
    for m in FAIL:
        print("FAIL:", m)
    for m in WARN:
        print("WARN:", m)
    return 2 if FAIL else (1 if WARN else 0)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
