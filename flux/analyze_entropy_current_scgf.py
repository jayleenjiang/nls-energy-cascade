#!/usr/bin/env python3
"""SCGF diagnostics for boundary-heat entropy-current gauges.

For two reservoirs with inverse temperatures beta_L and beta_R, the medium
entropy and the right-bath entropy current are

    Sigma_m = -beta_L Q_L - beta_R Q_R,
    Sigma_R = -(beta_R-beta_L) Q_R.

They differ by the boundary term beta_L Delta E (up to the recorded first-law
residual) and therefore target the same long-time SCGF when endpoint
exponential moments are controlled.  Sigma_R is analyzed because it removes
most of the endpoint-energy variance from the cold-bath representation.
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

from analyze_entropy_scgf import (
    aggregate_nonoverlapping,
    estimate_tau,
    infer_n,
    summarize,
    symmetry_rows,
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
    parser.add_argument("--k-step", type=float, default=0.025)
    parser.add_argument("--bootstraps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260828)
    return parser.parse_args()


def load_matrices(path: Path) -> dict[str, np.ndarray]:
    raw = np.loadtxt(
        path, delimiter=",", skiprows=1, usecols=(0, 1, 2, 3, 4, 5)
    )
    if raw.ndim != 2 or raw.shape[1] != 6 or raw.shape[0] == 0:
        raise ValueError(f"unexpected raw shape in {path}: {raw.shape}")
    if not np.all(np.isfinite(raw)):
        raise ValueError(f"nonfinite values in {path}")

    stream_ids = raw[:, 0].astype(np.int64)
    block_ids = raw[:, 1].astype(np.int64)
    streams, counts = np.unique(stream_ids, return_counts=True)
    if not np.array_equal(streams, np.arange(streams.size)):
        raise ValueError(f"stream IDs are not contiguous in {path}")
    if np.unique(counts).size != 1:
        raise ValueError(f"streams have unequal block counts in {path}")
    blocks_per_stream = int(counts[0])
    if not np.array_equal(stream_ids, np.repeat(streams, blocks_per_stream)):
        raise ValueError(f"rows are not grouped by stream in {path}")
    if not np.array_equal(
        block_ids, np.tile(np.arange(blocks_per_stream), streams.size)
    ):
        raise ValueError(f"block IDs are not ordered within streams in {path}")

    shape = (streams.size, blocks_per_stream)
    return {
        "q_left": raw[:, 2].reshape(shape),
        "q_right": raw[:, 3].reshape(shape),
        "delta_energy": raw[:, 4].reshape(shape),
        "medium": raw[:, 5].reshape(shape),
    }


def add_label(rows: list[dict[str, object]], label: str) -> None:
    for row in rows:
        row["observable"] = label


def main() -> None:
    args = parse_args()
    if not (args.temperature_left > 0 and args.temperature_right > 0):
        raise ValueError("temperatures must be positive")
    if args.max_multiple < 1 or args.bootstraps < 1:
        raise ValueError("invalid analysis parameters")

    beta_left = 1.0 / args.temperature_left
    beta_right = 1.0 / args.temperature_right
    delta_beta = beta_right - beta_left
    k_values = np.arange(0.0, 1.0 + 0.5 * args.k_step, args.k_step)

    all_points: list[dict[str, object]] = []
    all_symmetry: list[dict[str, object]] = []
    all_summary: list[dict[str, object]] = []
    identities: list[dict[str, object]] = []

    for path_index, path in enumerate(args.inputs):
        n = infer_n(path)
        matrices = load_matrices(path)
        balance = matrices["q_left"] + matrices["q_right"] - matrices["delta_energy"]
        right_current = -delta_beta * matrices["q_right"]
        identity_error = (
            right_current - matrices["medium"] - beta_left * matrices["delta_energy"]
        )
        identities.append(
            {
                "n": n,
                "samples": int(balance.size),
                "first_law_rms": float(np.sqrt(np.mean(balance**2))),
                "gauge_identity_rms": float(np.sqrt(np.mean(identity_error**2))),
                "medium_mean_rate": float(matrices["medium"].mean() / args.base_time),
                "right_current_mean_rate": float(right_current.mean() / args.base_time),
            }
        )

        observables = {
            "medium_entropy": matrices["medium"],
            "right_bath_entropy_current": right_current,
        }
        for label, base_matrix in observables.items():
            for multiple in range(1, args.max_multiple + 1):
                tau = args.base_time * multiple
                aggregated = aggregate_nonoverlapping(base_matrix, multiple)
                estimate = estimate_tau(
                    aggregated,
                    tau,
                    k_values,
                    args.bootstraps,
                    args.seed + 100003 * path_index + 1009 * multiple,
                )
                point_rows = [dict(row, n=n) for row in estimate.rows]
                add_label(point_rows, label)
                sym_rows = symmetry_rows(n, tau, k_values, estimate)
                add_label(sym_rows, label)
                summary = summarize(n, tau, k_values, estimate, sym_rows)
                summary["observable"] = label
                all_points.extend(point_rows)
                all_symmetry.extend(sym_rows)
                all_summary.append(summary)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "scgf_points.csv", all_points)
    write_csv(args.output_dir / "gc_symmetry_residuals.csv", all_symmetry)
    write_csv(args.output_dir / "scgf_reliability_summary.csv", all_summary)
    write_csv(args.output_dir / "gauge_identity_summary.csv", identities)

    for n in sorted({int(row["n"]) for row in all_points}):
        figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
        for tau, marker in [(20.0, "o"), (40.0, "s"), (80.0, "^")]:
            rows = [
                row
                for row in all_points
                if int(row["n"]) == n
                and row["observable"] == "right_bath_entropy_current"
                and math.isclose(float(row["tau"]), tau)
            ]
            if not rows:
                continue
            k = np.asarray([float(row["k"]) for row in rows])
            psi = np.asarray([float(row["psi"]) for row in rows])
            gate = np.asarray([bool(int(row["direct_reliability_gate"])) for row in rows])
            axes[0].plot(k, psi, marker=marker, ms=3, lw=1, label=rf"$t={tau:g}$")
            axes[0].scatter(k[gate], psi[gate], s=18)
            residual_rows = [
                row
                for row in all_symmetry
                if int(row["n"]) == n
                and row["observable"] == "right_bath_entropy_current"
                and math.isclose(float(row["tau"]), tau)
                and float(row["k"]) <= 0.5
            ]
            axes[1].plot(
                [float(row["k"]) for row in residual_rows],
                [float(row["symmetry_residual"]) for row in residual_rows],
                marker=marker,
                ms=3,
                lw=1,
                label=rf"$t={tau:g}$",
            )
        axes[0].set_xlabel(r"tilt $k$")
        axes[0].set_ylabel(r"$t^{-1}\log\mathbb{E}e^{-k\Sigma_R}$")
        axes[1].axhline(0.0, color="black", ls="--", lw=1)
        axes[1].set_xlabel(r"tilt $k\leq1/2$")
        axes[1].set_ylabel(r"$\psi_t(k)-\psi_t(1-k)$")
        for axis in axes:
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        figure.tight_layout()
        figure.savefig(args.output_dir / f"right_current_gc_n{n}.png", dpi=220)
        plt.close(figure)

    metadata = {
        "temperature_left": args.temperature_left,
        "temperature_right": args.temperature_right,
        "beta_left": beta_left,
        "beta_right": beta_right,
        "delta_beta": delta_beta,
        "base_time": args.base_time,
        "max_multiple": args.max_multiple,
        "k_step": args.k_step,
        "bootstraps": args.bootstraps,
    }
    with (args.output_dir / "analysis_metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata), lineterminator="\n")
        writer.writeheader()
        writer.writerow(metadata)

    print(f"Wrote {len(all_points)} SCGF points for {len(args.inputs)} chain lengths")


if __name__ == "__main__":
    main()
