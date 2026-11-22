#!/usr/bin/env python3
"""Independent random-walk Metropolis reference for the n-mode Gibbs density at T=5
(port of the n=4 gibbs_reference.py: same chains, burn, steps, thinning and step-size adaptation).

State x = (I1, I2, I3, theta2, theta3), theta_j = 2(phi_j - phi_{j-1}); density exp(-H/(2T)).
Usage: python3 gibbs_reference_n.py SEED OUT.bin [n=3]   (also writes OUT_diagnostics.json)
"""
from __future__ import annotations

import json
import struct
import sys
from pathlib import Path

import numpy as np

SEED = int(sys.argv[1])
OUTPUT = Path(sys.argv[2])
REPORT = OUTPUT.with_name(OUTPUT.stem + "_diagnostics.json")
TEMPERATURE = 5.0
CHAINS = 64
BURN_STEPS = 100_000
SAMPLE_STEPS = 200_000
THIN = 20
NA = int(sys.argv[3]) if len(sys.argv) > 3 else 3
NT = NA - 1


def wrap_pi(theta):
    return (theta + np.pi) % (2.0 * np.pi) - np.pi


def hamiltonian(y, theta):
    actions = np.exp(y)
    mass = actions.sum(axis=1)
    return (mass * mass - 0.5 * (actions * actions).sum(axis=1)
            + 2.0 * (actions[:, :-1] * actions[:, 1:] * np.cos(theta)).sum(axis=1))


def log_target(y, theta):
    # The sum(y) term is the Jacobian for I_j = exp(y_j).
    return -hamiltonian(y, theta) / (2.0 * TEMPERATURE) + y.sum(axis=1)


def main():
    if OUTPUT.exists():
        raise SystemExit(f"refusing to overwrite {OUTPUT}")
    rng = np.random.default_rng(SEED)
    y = np.log(rng.gamma(shape=2.0, scale=0.6, size=(CHAINS, NA)))
    theta = rng.uniform(-np.pi, np.pi, size=(CHAINS, NT))
    current = log_target(y, theta)
    scale_y, scale_theta = 0.20, 0.35
    accepted_burn = 0
    for step in range(BURN_STEPS):
        py = y + scale_y * rng.normal(size=y.shape)
        pt = wrap_pi(theta + scale_theta * rng.normal(size=theta.shape))
        prop = log_target(py, pt)
        acc = np.log(rng.random(CHAINS)) < prop - current
        y[acc], theta[acc], current[acc] = py[acc], pt[acc], prop[acc]
        accepted_burn += int(acc.sum())
        if (step + 1) % 1000 == 0:
            rate = accepted_burn / ((step + 1) * CHAINS)
            adj = np.exp(0.05 * (rate - 0.28))
            scale_y = float(np.clip(scale_y * adj, 0.02, 1.0))
            scale_theta = float(np.clip(scale_theta * adj, 0.04, 1.5))

    draws = SAMPLE_STEPS // THIN
    samples = np.empty((draws * CHAINS, NA + NT), dtype=np.float32)
    accepted, row = 0, 0
    chain_sum = np.zeros((CHAINS, NA))
    saved = 0
    for step in range(SAMPLE_STEPS):
        py = y + scale_y * rng.normal(size=y.shape)
        pt = wrap_pi(theta + scale_theta * rng.normal(size=theta.shape))
        prop = log_target(py, pt)
        acc = np.log(rng.random(CHAINS)) < prop - current
        y[acc], theta[acc], current[acc] = py[acc], pt[acc], prop[acc]
        accepted += int(acc.sum())
        if (step + 1) % THIN == 0:
            a = np.exp(y)
            samples[row:row + CHAINS, :NA] = a
            samples[row:row + CHAINS, NA:] = theta
            chain_sum += a
            saved += 1
            row += CHAINS

    with OUTPUT.open("wb") as h:
        h.write(("NLS%dD01" % (NA + NT)).encode() + b"\0")
        h.write(struct.pack("<QII", samples.shape[0], NA + NT, 4))
        samples.tofile(h)
    diag = {"seed": SEED, "temperature": TEMPERATURE, "chains": CHAINS, "burn_steps": BURN_STEPS,
            "sample_steps": SAMPLE_STEPS, "thin": THIN, "saved_samples": int(samples.shape[0]),
            "burn_acceptance": accepted_burn / (BURN_STEPS * CHAINS),
            "sample_acceptance": accepted / (SAMPLE_STEPS * CHAINS),
            "final_scale_y": scale_y, "final_scale_theta": scale_theta,
            "mean_actions": samples[:, :NA].mean(axis=0).astype(float).tolist(),
            "max_between_chain_action_mean_range": float(np.ptp(chain_sum / saved, axis=0).max())}
    REPORT.write_text(json.dumps(diag, indent=2) + "\n")
    print(json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()
