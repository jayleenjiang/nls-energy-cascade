#!/usr/bin/env python3
"""Gibbs expectation of the site actions for the n-mode chain, exp(-E/T) with E = H/2 (the invariant measure of the canonical bath), by Metropolis in Cartesian
coordinates (vectorised over independent chains).  H = M^2 - 1/2 sum I_j^2 + 2 sum Re(c_j^2 conj(c_{j-1})^2).
usage: gibbs_chain.py n T chains sweeps seed OUT.csv   (T is the bath temperature; the Metropolis uses exp(-H/(2T)))
"""
import sys
import numpy as np

n, T, nch, nsweep, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(seed)
amp = np.sqrt(np.sqrt(2.0 * T / n))
c = amp * (rng.standard_normal((nch, n)) + 1j * rng.standard_normal((nch, n))) / np.sqrt(2) * np.sqrt(2)
I = np.abs(c) ** 2
M = I.sum(1)
c2 = c * c


def local_energy(j, cj):
    """part of H depending on c_j (M^2 handled separately): -I_j^2/2 + 2 Re(c_j^2 (conj(c_{j-1})^2 + conj(c_{j+1})^2))"""
    Ij = np.abs(cj) ** 2
    nb = np.zeros(nch, dtype=complex)
    if j > 0: nb += np.conj(c2[:, j - 1])
    if j + 1 < n: nb += np.conj(c2[:, j + 1])
    return -0.5 * Ij * Ij + 2.0 * np.real(cj * cj * nb), Ij


step = 0.35 * amp
burn = nsweep // 5
acc = 0; tot = 0
Isum = np.zeros(n); Isum2 = np.zeros(n); cnt = 0
chain_I = np.zeros((nch, n)); chain_cnt = 0
for sweep in range(nsweep):
    for j in rng.permutation(n):
        cj = c[:, j]
        prop = cj + step * (rng.standard_normal(nch) + 1j * rng.standard_normal(nch))
        e_old, Ij = local_energy(j, cj)
        e_new, Ip = local_energy(j, prop)
        Mn = M - Ij + Ip
        dH = (e_new - e_old) + (Mn * Mn - M * M)
        ok = np.log(rng.random(nch)) < -dH / (2.0 * T)
        c[ok, j] = prop[ok]; c2[ok, j] = prop[ok] ** 2; M[ok] = Mn[ok]
        acc += ok.sum(); tot += nch
    if sweep >= burn:
        I = np.abs(c) ** 2
        Isum += I.sum(0); cnt += nch
        chain_I += I; chain_cnt += 1
mean = Isum / cnt
per_chain = chain_I / chain_cnt
se = per_chain.std(0, ddof=1) / np.sqrt(nch)
Ibar = mean.mean(); Ibar_se = per_chain.mean(1).std(ddof=1) / np.sqrt(nch)
print(f"n={n} T={T} chains={nch} sweeps={nsweep} acceptance={acc/tot:.3f}  Ibar={Ibar:.5f} +- {Ibar_se:.5f}  sqrt(n)*Ibar={np.sqrt(n)*Ibar:.4f}  <M>={np.sqrt(n)*0+ (per_chain.sum(1)).mean():.3f}")
with open(sys.argv[6], "w") as f:
    f.write("j,mean_I_gibbs,se\n")
    for j in range(n): f.write(f"{j+1},{mean[j]:.6f},{se[j]:.6f}\n")
