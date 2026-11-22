#!/usr/bin/env python3
"""Draw query states for NLS7D_stream_kde from an existing NLS7D01 sample file.

Usage: python3 make_queries.py SOURCE.bin OUT.bin COUNT [SEED]
"""
import struct, sys
import numpy as np

def read_nls7d(path):
    with open(path, "rb") as f:
        if f.read(8) != b"NLS7D01\0":
            raise ValueError(f"bad header: {path}")
        n, d, b = struct.unpack("<QII", f.read(16))
        if d != 7 or b != 4:
            raise ValueError(f"unsupported layout: {path}")
        x = np.fromfile(f, "<f4")
    if x.size != n * 7:
        raise ValueError(f"truncated: {path}")
    return x.reshape(n, 7)

def write_nls7d(path, x):
    with open(path, "wb") as f:
        f.write(b"NLS7D01\0")
        f.write(struct.pack("<QII", len(x), 7, 4))
        f.write(np.asarray(x, "<f4").tobytes())

if __name__ == "__main__":
    src, out, count = sys.argv[1], sys.argv[2], int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    x = read_nls7d(src)
    idx = np.random.default_rng(seed).choice(len(x), min(count, len(x)), replace=False)
    write_nls7d(out, x[idx])
    print(f"wrote {out}: {len(idx)} queries from {src} (seed {seed})")
