#!/usr/bin/env python3
"""Generator Galerkin projection (gEDMD) for the n=3 chain, paper conventions.

State x = (I1, I2, I3, theta1, theta3), theta1 = 2(phi1-phi2), theta3 = 2(phi3-phi2).
Generator (2.18):  L g = sum_i f_i d_i g + sum_i a_i d_i^2 g,  a_i = D_ii / 2:
  a_I1 = 4 gamma T1 I1, a_th1 = 4 gamma T1 / I1, a_I3 = 4 gamma T3 I3, a_th3 = 4 gamma T3 / I3.
Drift (from (1.10)-(1.14) with n=3):
  dH1 = 2M - I1 + 2 I2 cos th1, dH2 = 2M - I2 + 2 I1 cos th1 + 2 I3 cos th3, dH3 = 2M - I3 + 2 I2 cos th3
  f_I1 = 4 I1 I2 sin th1 + 2 gamma (2 T1 - I1 dH1)
  f_I2 = -4 I2 (I1 sin th1 + I3 sin th3)
  f_I3 = 4 I2 I3 sin th3 + 2 gamma (2 T3 - I3 dH3)
  f_th1 = 2 (dH1 - dH2) + 4 gamma I2 sin th1
  f_th3 = 2 (dH3 - dH2) + 4 gamma I2 sin th3
Dictionary: phi = prod_j z_j^{p_j} h_j^{q_j} * T(m1 th1 + m3 th3),  z = I/s,  h = I/(s+I),
with gate q = (|m1|, |m1+m3|, |m3|) so that every function is smooth in the Cartesian amplitudes.
Galerkin: A = E[phi phi^T], B = E[phi (L phi)^T];  L psi = lambda psi  ->  B c = lambda A c.
"""
from __future__ import annotations
import itertools, json, struct, sys, time
import numpy as np

GAMMA = 0.1


def load(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"NLS3X01\0"
        ns, nt = struct.unpack("<QQ", f.read(16))
        x = np.fromfile(f, "<f8")
    return x.reshape(ns, nt, 5)


def make_dictionary(d, M):
    acts = [p for p in itertools.product(range(d + 1), repeat=3) if sum(p) <= d]
    modes = [(0, 0, "c")]
    for m1 in range(0, M + 1):
        for m3 in range(-M, M + 1):
            if m1 == 0 and m3 <= 0:
                continue
            modes += [(m1, m3, "c"), (m1, m3, "s")]
    terms = []
    for (m1, m3, kind) in modes:
        q = (abs(m1), abs(m1 + m3), abs(m3))
        for p in acts:
            terms.append((p, q, m1, m3, kind))
    return terms


def factor(I, s, p, q):
    """u = (I/s)^p (I/(s+I))^q and its first two derivatives."""
    z = I / s
    h = I / (s + I)
    h1 = s / (s + I) ** 2
    h2 = -2.0 * s / (s + I) ** 3
    zp = z ** p
    zp1 = p * z ** (p - 1) / s if p >= 1 else 0.0 * I
    zp2 = p * (p - 1) * z ** (p - 2) / s ** 2 if p >= 2 else 0.0 * I
    hq = h ** q
    hq1 = q * h ** (q - 1) * h1 if q >= 1 else 0.0 * I
    hq2 = (q * (q - 1) * h ** (q - 2) * h1 ** 2 if q >= 2 else 0.0 * I) + (q * h ** (q - 1) * h2 if q >= 1 else 0.0 * I)
    u = zp * hq
    du = zp1 * hq + zp * hq1
    d2u = zp2 * hq + 2.0 * zp1 * hq1 + zp * hq2
    return u, du, d2u


def drift(x, T1, T3):
    I1, I2, I3, t1, t3 = x.T
    c1, s1, c3, s3 = np.cos(t1), np.sin(t1), np.cos(t3), np.sin(t3)
    M = I1 + I2 + I3
    dH1 = 2 * M - I1 + 2 * I2 * c1
    dH2 = 2 * M - I2 + 2 * I1 * c1 + 2 * I3 * c3
    dH3 = 2 * M - I3 + 2 * I2 * c3
    f = np.empty_like(x)
    f[:, 0] = 4 * I1 * I2 * s1 + 2 * GAMMA * (2 * T1 - I1 * dH1)
    f[:, 1] = -4 * I2 * (I1 * s1 + I3 * s3)
    f[:, 2] = 4 * I2 * I3 * s3 + 2 * GAMMA * (2 * T3 - I3 * dH3)
    f[:, 3] = 2 * (dH1 - dH2) + 4 * GAMMA * I2 * s1
    f[:, 4] = 2 * (dH3 - dH2) + 4 * GAMMA * I2 * s3
    a = np.stack([4 * GAMMA * T1 * I1, 4 * GAMMA * T3 * I3, 4 * GAMMA * T1 / I1, 4 * GAMMA * T3 / I3], 1)
    return f, a  # a columns: I1, I3, th1, th3


def evaluate(x, terms, s, T1, T3, want_grad=False):
    """Phi, L Phi (n x K); optionally the noisy-coordinate gradient (n x K x 4)."""
    n, K = len(x), len(terms)
    f, a = drift(x, T1, T3)
    cache = {}
    def fac(j, p, q):
        key = (j, p, q)
        if key not in cache:
            cache[key] = factor(x[:, j], s[j], p, q)
        return cache[key]
    Phi = np.empty((n, K)); LPhi = np.empty((n, K))
    G = np.empty((n, K, 4)) if want_grad else None
    for k, (p, q, m1, m3, kind) in enumerate(terms):
        (u1, d1, e1), (u2, d2, e2), (u3, d3, e3) = fac(0, p[0], q[0]), fac(1, p[1], q[1]), fac(2, p[2], q[2])
        arg = m1 * x[:, 3] + m3 * x[:, 4]
        if kind == "c":
            T, dT = np.cos(arg), -np.sin(arg)
        else:
            T, dT = np.sin(arg), np.cos(arg)
        U = u1 * u2 * u3
        phi = U * T
        g_I1 = d1 * u2 * u3 * T
        g_I2 = u1 * d2 * u3 * T
        g_I3 = u1 * u2 * d3 * T
        g_t1 = U * m1 * dT
        g_t3 = U * m3 * dT
        Lphi = (f[:, 0] * g_I1 + f[:, 1] * g_I2 + f[:, 2] * g_I3 + f[:, 3] * g_t1 + f[:, 4] * g_t3
                + a[:, 0] * e1 * u2 * u3 * T + a[:, 1] * u1 * u2 * e3 * T
                - (a[:, 2] * m1 * m1 + a[:, 3] * m3 * m3) * phi)
        Phi[:, k] = phi; LPhi[:, k] = Lphi
        if want_grad:
            G[:, k, 0] = g_I1; G[:, k, 1] = g_I3; G[:, k, 2] = g_t1; G[:, k, 3] = g_t3
    return Phi, LPhi, f, a, G


def accumulate(X, terms, s, T1, T3, chunk=20000):
    """X: (streams, nt, 5). Returns per-half A, B, mean(L phi) and counts."""
    ns = X.shape[0]
    halves = [X[: ns // 2].reshape(-1, 5), X[ns // 2:].reshape(-1, 5)]
    K = len(terms)
    out = []
    for H in halves:
        A = np.zeros((K, K)); B = np.zeros((K, K)); m = np.zeros(K); m2 = np.zeros(K)
        for i in range(0, len(H), chunk):
            Phi, LPhi, *_ = evaluate(H[i:i + chunk], terms, s, T1, T3)
            A += Phi.T @ Phi; B += Phi.T @ LPhi; m += LPhi.sum(0); m2 += (LPhi ** 2).sum(0)
        n = len(H)
        out.append(dict(A=A / n, B=B / n, mL=m / n, sdL=np.sqrt(m2 / n - (m / n) ** 2), n=n))
    return out


def solve(A, B, cut=1e-10):
    w, V = np.linalg.eigh(A)
    keep = w > cut * w.max()
    W = V[:, keep] / np.sqrt(w[keep])
    Mw = W.T @ B @ W
    lam, Vw = np.linalg.eig(Mw)
    order = np.argsort(-lam.real)
    return lam[order], W @ Vw[:, order], int(keep.sum())
