#!/usr/bin/env python3
"""Independent Metropolis reference for the equilibrium Gibbs density.

The target measure is

    rho(I, phi) dI dphi proportional to exp[-H(I,phi)/(2T)] dI dphi,

with I_j >= 0 and one global phase fixed because H depends only on phase
differences.  Actions are sampled in logarithmic coordinates, so the log-target
includes the Jacobian sum(log I_j).

The output is intended as an independent equilibrium reference for validating
the SDE integrator; it does not reuse any stochastic-dynamics code.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--temperature", type=float, default=2.0)
    parser.add_argument("--chains", type=int, default=256)
    parser.add_argument("--burn-sweeps", type=int, default=10_000)
    parser.add_argument("--samples-per-chain", type=int, default=5_000)
    parser.add_argument("--thin", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260619)
    parser.add_argument("--action-proposal", type=float, default=0.45)
    parser.add_argument("--phase-proposal", type=float, default=0.9)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def wrap_pi(value: np.ndarray) -> np.ndarray:
    return value - 2.0 * np.pi * np.rint(value / (2.0 * np.pi))


def hamiltonian(action: np.ndarray, phase: np.ndarray) -> np.ndarray:
    total = np.sum(action, axis=1)
    onsite = total**2 - 0.5 * np.sum(action**2, axis=1)
    angle = 2.0 * (phase[:, 1:] - phase[:, :-1])
    coupling = 2.0 * np.sum(
        action[:, :-1] * action[:, 1:] * np.cos(angle), axis=1
    )
    return onsite + coupling


def action_gradient(action: np.ndarray, phase: np.ndarray) -> np.ndarray:
    total = np.sum(action, axis=1)
    gradient = 2.0 * total[:, None] - action
    angle = 2.0 * (phase[:, 1:] - phase[:, :-1])
    cosine = np.cos(angle)
    gradient[:, :-1] += 2.0 * action[:, 1:] * cosine
    gradient[:, 1:] += 2.0 * action[:, :-1] * cosine
    return gradient


def log_target(
    log_action: np.ndarray, phase: np.ndarray, temperature: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    action = np.exp(log_action)
    energy = hamiltonian(action, phase)
    value = -energy / (2.0 * temperature) + np.sum(log_action, axis=1)
    return value, action, energy


def gelman_rubin(chain_sum: np.ndarray, chain_sum_sq: np.ndarray, n: int) -> np.ndarray:
    means = chain_sum / n
    variances = (chain_sum_sq - chain_sum**2 / n) / (n - 1)
    within = np.mean(variances, axis=0)
    between = n * np.var(means, axis=0, ddof=1)
    variance_hat = ((n - 1) / n) * within + between / n
    return np.sqrt(variance_hat / within)


def main() -> None:
    args = parse_args()
    if args.n < 2:
        raise ValueError("n must be at least 2")
    if args.temperature <= 0:
        raise ValueError("temperature must be positive")
    if args.chains < 4 or args.samples_per_chain < 2:
        raise ValueError("need at least four chains and two samples per chain")

    rng = np.random.default_rng(args.seed)
    initial_action = math.sqrt(args.temperature / args.n)
    log_action = (
        math.log(initial_action)
        + 0.25 * rng.standard_normal((args.chains, args.n))
    )
    phase = rng.uniform(-np.pi, np.pi, size=(args.chains, args.n))
    phase[:, 0] = 0.0
    current_log_target, action, energy = log_target(
        log_action, phase, args.temperature
    )

    action_accept = np.zeros(args.n, dtype=np.int64)
    action_attempt = np.zeros(args.n, dtype=np.int64)
    phase_accept = np.zeros(args.n, dtype=np.int64)
    phase_attempt = np.zeros(args.n, dtype=np.int64)

    chain_sum = np.zeros((args.chains, args.n), dtype=float)
    chain_sum_sq = np.zeros((args.chains, args.n), dtype=float)
    chain_energy_sum = np.zeros(args.chains, dtype=float)
    chain_mass_sum = np.zeros(args.chains, dtype=float)
    chain_identity_sum = np.zeros((args.chains, args.n), dtype=float)

    total_sweeps = (
        args.burn_sweeps + args.samples_per_chain * args.thin
    )
    stored = 0

    for sweep in range(total_sweeps):
        for j in range(args.n):
            proposal_log_action = log_action.copy()
            proposal_log_action[:, j] += (
                args.action_proposal * rng.standard_normal(args.chains)
            )
            proposal_target, proposal_action, proposal_energy = log_target(
                proposal_log_action, phase, args.temperature
            )
            accept = np.log(rng.random(args.chains)) < (
                proposal_target - current_log_target
            )
            log_action[accept, j] = proposal_log_action[accept, j]
            action[accept] = proposal_action[accept]
            energy[accept] = proposal_energy[accept]
            current_log_target[accept] = proposal_target[accept]
            action_accept[j] += int(np.sum(accept))
            action_attempt[j] += args.chains

        # phi_0 is fixed to remove the irrelevant global phase.
        for j in range(1, args.n):
            proposal_phase = phase.copy()
            proposal_phase[:, j] = wrap_pi(
                proposal_phase[:, j]
                + args.phase_proposal * rng.standard_normal(args.chains)
            )
            proposal_target, proposal_action, proposal_energy = log_target(
                log_action, proposal_phase, args.temperature
            )
            accept = np.log(rng.random(args.chains)) < (
                proposal_target - current_log_target
            )
            phase[accept, j] = proposal_phase[accept, j]
            action[accept] = proposal_action[accept]
            energy[accept] = proposal_energy[accept]
            current_log_target[accept] = proposal_target[accept]
            phase_accept[j] += int(np.sum(accept))
            phase_attempt[j] += args.chains

        if sweep >= args.burn_sweeps and (
            sweep - args.burn_sweeps
        ) % args.thin == 0:
            gradient = action_gradient(action, phase)
            chain_sum += action
            chain_sum_sq += action**2
            chain_energy_sum += energy
            chain_mass_sum += np.sum(action, axis=1)
            chain_identity_sum += action * gradient
            stored += 1

        if (sweep + 1) % max(1, total_sweeps // 10) == 0:
            print(f"completed {sweep + 1}/{total_sweeps} sweeps")

    if stored != args.samples_per_chain:
        raise RuntimeError(
            f"stored {stored} samples per chain, expected "
            f"{args.samples_per_chain}"
        )

    chain_means = chain_sum / stored
    profile_mean = np.mean(chain_means, axis=0)
    profile_se = np.std(chain_means, axis=0, ddof=1) / math.sqrt(args.chains)
    profile_rhat = gelman_rubin(chain_sum, chain_sum_sq, stored)
    identity_chain_means = chain_identity_sum / stored
    identity_mean = np.mean(identity_chain_means, axis=0)
    identity_se = (
        np.std(identity_chain_means, axis=0, ddof=1)
        / math.sqrt(args.chains)
    )

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    profile_path = output_prefix.with_name(
        output_prefix.name + "_profile.csv"
    )
    with profile_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "mean_action",
                "between_chain_se",
                "rhat",
                "mean_I_dH_dI",
                "I_dH_dI_between_chain_se",
            ]
        )
        for j in range(args.n):
            writer.writerow(
                [
                    j,
                    profile_mean[j],
                    profile_se[j],
                    profile_rhat[j],
                    identity_mean[j],
                    identity_se[j],
                ]
            )

    action_acceptance = action_accept / action_attempt
    phase_acceptance = np.divide(
        phase_accept[1:],
        phase_attempt[1:],
        out=np.zeros(args.n - 1, dtype=float),
        where=phase_attempt[1:] > 0,
    )
    summary = {
        "target": "exp(-H/(2T)) dI dphi",
        "n": args.n,
        "temperature": args.temperature,
        "chains": args.chains,
        "burn_sweeps": args.burn_sweeps,
        "samples_per_chain": args.samples_per_chain,
        "thin": args.thin,
        "seed": args.seed,
        "action_proposal": args.action_proposal,
        "phase_proposal": args.phase_proposal,
        "action_acceptance_min": float(np.min(action_acceptance)),
        "action_acceptance_max": float(np.max(action_acceptance)),
        "phase_acceptance_min": float(np.min(phase_acceptance)),
        "phase_acceptance_max": float(np.max(phase_acceptance)),
        "max_rhat": float(np.max(profile_rhat)),
        "mean_energy": float(np.mean(chain_energy_sum / stored)),
        "energy_between_chain_se": float(
            np.std(chain_energy_sum / stored, ddof=1)
            / math.sqrt(args.chains)
        ),
        "mean_total_action": float(np.mean(chain_mass_sum / stored)),
        "total_action_between_chain_se": float(
            np.std(chain_mass_sum / stored, ddof=1)
            / math.sqrt(args.chains)
        ),
        "integration_by_parts_target": 2.0 * args.temperature,
        "max_abs_integration_by_parts_z": float(
            np.max(
                np.abs(
                    (identity_mean - 2.0 * args.temperature) / identity_se
                )
            )
        ),
        "reflection_max_abs_z": float(
            np.max(
                np.abs(
                    profile_mean
                    - profile_mean[::-1]
                )
                / np.sqrt(profile_se**2 + profile_se[::-1] ** 2)
            )
        ),
    }
    summary_path = output_prefix.with_name(
        output_prefix.name + "_summary.json"
    )
    with summary_path.open("w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    print(json.dumps(summary, indent=2))
    print(f"wrote {profile_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
