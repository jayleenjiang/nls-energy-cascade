#!/usr/bin/env python3
"""Draw query states for nls5d_stream from an existing NLS5D01 sample file (port of make_queries.py).

Usage: python3 make_queries_n3.py SOURCE.bin OUT.bin COUNT [SEED]
"""
import struct
import sys
from pathlib import Path

import numpy as np


def read_nls5d(path):
    with open(path, "rb") as f:
        if f.read(8) != b"NLS5D01\0":
            raise ValueError(f"bad header: {path}")
        n, d, b = struct.unpack("<QII", f.read(16))
        if d != 5 or b != 4:
            raise ValueError(f"unsupported layout: {path}")
        x = np.fromfile(f, "<f4")
    if x.size != n * 5:
        raise ValueError(f"truncated: {path}")
    return x.reshape(n, 5)


def write_nls5d(path, x):
    with open(path, "wb") as f:
        f.write(b"NLS5D01\0")
        f.write(struct.pack("<QII", len(x), 5, 4))
        f.write(np.asarray(x, "<f4").tobytes())


if __name__ == "__main__":
    src, out, count = sys.argv[1], sys.argv[2], int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    if Path(out).exists():
        raise SystemExit(f"refusing to overwrite {out}")
    x = read_nls5d(src)
    idx = np.random.default_rng(seed).choice(len(x), min(count, len(x)), replace=False)
    write_nls5d(out, x[idx])
    print(f"wrote {out}: {len(idx)} queries from {src} (seed {seed})")
