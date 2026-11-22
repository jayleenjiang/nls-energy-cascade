#!/usr/bin/env python3
"""Equilibrium known-answer check for a streamed KDE run at T1=T4=T.

Compares log KDE with the exact Gibbs log density -E/T (E = H/2) at the query
states, after removing the mean offset (centred log RMSE; smaller is better).
Usage: python3 eval_eq.py QUERIES.bin RUN_kde.bin [T=5]
"""
import struct, sys
import numpy as np
from make_queries import read_nls7d

def read_kde(path):
    with open(path, "rb") as f:
        if f.read(8) != b"NLS7KDE1":
            raise ValueError(f"bad header: {path}")
        q, n = struct.unpack("<QQ", f.read(16))
        return np.fromfile(f, "<f8"), int(n)

def energy(x):
    I, th = x[:, :4].astype(np.float64), x[:, 4:].astype(np.float64)
    m = I.sum(1)
    return (0.5 * m * m - 0.25 * (I * I).sum(1) + I[:, 0] * I[:, 1] * np.cos(th[:, 0])
            + I[:, 1] * I[:, 2] * np.cos(th[:, 1]) + I[:, 2] * I[:, 3] * np.cos(th[:, 2]))

if __name__ == "__main__":
    q = read_nls7d(sys.argv[1]); dens, n = read_kde(sys.argv[2])
    T = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
    good = dens > 0
    e = np.log(dens[good]) + energy(q[good]) / T
    e -= e.mean()
    print(f"samples={n} queries={len(q)} zero_density={int((~good).sum())} "
          f"centred_log_RMSE={np.sqrt((e * e).mean()):.4f}")
