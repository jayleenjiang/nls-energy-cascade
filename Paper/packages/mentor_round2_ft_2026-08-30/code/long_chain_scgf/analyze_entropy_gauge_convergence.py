#!/usr/bin/env python3
"""Joint-bootstrap the medium-entropy/right-current gauge difference.

For the recorded heat convention,

    Sigma_R - Sigma_m = beta_L Delta E + O(first-law residual).

The two finite-time SCGFs need not be equal, but their difference should decay
at long times wherever direct exponential averages retain support.  Both SCGFs
are evaluated from the same block groups and the difference is bootstrapped
with the same stream multiplicities, so the endpoint correlation is preserved.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_entropy_current_scgf import load_matrices
from analyze_entropy_scgf import (
    aggregate_nonoverlapping,
    estimate_tau,
    infer_n,
    write_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--temperature-left", type=float, default=10.0)
    parser.add_argument("--temperature-right", type=float, default=2.0)
    parser.add_argument("--base-time", type=float, default=20.0)
    parser.add_argument("--max-multiple", type=int, default=10)
    parser.add_argument(
        "--k-values", nargs="+", type=float, default=[0.05, 0.1, 0.2, 0.3]
    )
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260829)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.temperature_left <= 0.0 or args.temperature_right <= 0.0:
        raise ValueError("temperatures must be positive")
    if args.max_multiple < 1 or args.bootstraps < 1:
        raise ValueError("invalid analysis parameters")
    k_values = np.asarray(args.k_values, dtype=float)
    if np.any(k_values < 0.0) or np.any(k_values > 1.0):
        raise ValueError("k values must lie in [0,1]")

    beta_left = 1.0 / args.temperature_left
    delta_beta = 1.0 / args.temperature_right - beta_left
    rows: list[dict[str, object]] = []
    identity_rows: list[dict[str, object]] = []

    for path_index, path in enumerate(args.inputs):
        n = infer_n(path)
        matrices = load_matrices(path)
        right = -delta_beta * matrices["q_right"]
        medium = matrices["medium"]
        identity_error = (
            right - medium - beta_left * matrices["delta_energy"]
        )
        identity_rows.append(
            {
                "n": n,
                "samples": int(identity_error.size),
                "gauge_identity_rms": float(
                    np.sqrt(np.mean(np.square(identity_error)))
                ),
                "gauge_identity_max_abs": float(
                    np.max(np.abs(identity_error))
                ),
            }
        )

        for multiple in range(1, args.max_multiple + 1):
            tau = args.base_time * multiple
            medium_agg = aggregate_nonoverlapping(medium, multiple)
            right_agg = aggregate_nonoverlapping(right, multiple)
            shared_seed = args.seed + 100003 * path_index + 1009 * multiple
            medium_est = estimate_tau(
                medium_agg, tau, k_values, args.bootstraps, shared_seed
            )
            right_est = estimate_tau(
                right_agg, tau, k_values, args.bootstraps, shared_seed
            )
            for index, k in enumerate(k_values):
                difference_bootstrap = (
                    right_est.psi_bootstrap[:, index]
                    - medium_est.psi_bootstrap[:, index]
                )
                low, high = np.quantile(difference_bootstrap, [0.025, 0.975])
                medium_psi = float(medium_est.rows[index]["psi"])
                right_psi = float(right_est.rows[index]["psi"])
                difference = right_psi - medium_psi
                rows.append(
                    {
                        "n": n,
                        "tau": tau,
                        "inverse_tau": 1.0 / tau,
                        "k": float(k),
                        "psi_medium": medium_psi,
                        "psi_right": right_psi,
                        "right_minus_medium": difference,
                        "difference_ci_low": float(low),
                        "difference_ci_high": float(high),
                        "tau_times_difference": tau * difference,
                        "medium_reliability_gate": int(
                            medium_est.reliable[index]
                        ),
                        "right_reliability_gate": int(
                            right_est.reliable[index]
                        ),
                        "paired_reliability_gate": int(
                            medium_est.reliable[index]
                            and right_est.reliable[index]
                        ),
                    }
                )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "gauge_convergence.csv", rows)
    write_csv(args.output_dir / "gauge_identity.csv", identity_rows)

    for n in sorted({int(row["n"]) for row in rows}):
        figure, axis = plt.subplots(figsize=(6.8, 4.6))
        for k in k_values:
            selected = sorted(
                [
                    row for row in rows
                    if int(row["n"]) == n
                    and math.isclose(float(row["k"]), float(k))
                    and int(row["paired_reliability_gate"]) == 1
                ],
                key=lambda row: float(row["inverse_tau"]),
            )
            if not selected:
                continue
            x = np.asarray([float(row["inverse_tau"]) for row in selected])
            y = np.asarray([float(row["right_minus_medium"]) for row in selected])
            low = np.asarray([float(row["difference_ci_low"]) for row in selected])
            high = np.asarray([float(row["difference_ci_high"]) for row in selected])
            axis.errorbar(
                x, y, yerr=np.vstack([y - low, high - y]), marker="o",
                capsize=2.5, label=rf"$k={k:g}$",
            )
        axis.axhline(0.0, color="black", linestyle="--", linewidth=1.0)
        axis.set_xlabel(r"inverse observation time $1/t$")
        axis.set_ylabel(r"$\psi_t^{R}(k)-\psi_t^{m}(k)$")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
        figure.tight_layout()
        figure.savefig(args.output_dir / f"gauge_convergence_n{n}.png", dpi=220)
        plt.close(figure)

    metadata = {
        "temperature_left": args.temperature_left,
        "temperature_right": args.temperature_right,
        "base_time": args.base_time,
        "max_multiple": args.max_multiple,
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "k_values": " ".join(f"{value:g}" for value in k_values),
    }
    with (args.output_dir / "metadata.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metadata))
        writer.writeheader()
        writer.writerow(metadata)
    print(f"Wrote {len(rows)} paired gauge comparisons")


if __name__ == "__main__":
    main()
