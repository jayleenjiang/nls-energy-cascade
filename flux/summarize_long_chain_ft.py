#!/usr/bin/env python3
"""Compile the pre-specified long-chain GC evidence into one audited table."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment_root", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def close(left: object, right: object) -> bool:
    return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=1.0e-12)


def one(rows: list[dict[str, str]], **criteria: object) -> dict[str, str]:
    matches = [
        row for row in rows
        if all(close(row[key], value) for key, value in criteria.items())
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one row for {criteria}, found {len(matches)}")
    return matches[0]


def main() -> None:
    args = parse_args()
    specs = read_csv(args.experiment_root / "FINAL_RESULT_SPEC.csv")
    summary: list[dict[str, object]] = []

    for spec in specs:
        n = int(spec["n"])
        k = float(spec["k"])
        partner_k = 1.0 - k
        clone_count = int(spec["clone_count"])
        horizon = float(spec["horizon"])
        analysis = args.experiment_root / spec["analysis_dir"]

        aggregates = read_csv(analysis / "controlled_aggregate.csv")
        pairs = read_csv(analysis / "paired_gc_residuals.csv")
        time_rows = read_csv(analysis / "time_convergence_pairs.csv")
        population_pairs = read_csv(analysis / "population_convergence_pairs.csv")
        population_members = read_csv(
            analysis / "population_convergence_members.csv"
        )

        low = one(
            aggregates, n=n, clone_count=clone_count, horizon=horizon, k=k
        )
        high = one(
            aggregates, n=n, clone_count=clone_count, horizon=horizon,
            k=partner_k,
        )
        pair = one(
            pairs, n=n, clone_count=clone_count, horizon=horizon, k=k
        )
        time = one(
            time_rows, n=n, clone_count=clone_count, final_horizon=horizon, k=k
        )
        pop_pair = one(
            population_pairs, n=n, upper_clone_count=clone_count,
            horizon=horizon, k=k,
        )
        pop_low = one(
            population_members, n=n, upper_clone_count=clone_count,
            horizon=horizon, k=k,
        )
        pop_high = one(
            population_members, n=n, upper_clone_count=clone_count,
            horizon=horizon, k=partner_k,
        )

        gate = all(
            [
                int(low["support_gate"]),
                int(high["support_gate"]),
                int(pair["paired_gc_gate"]),
                int(time["time_convergence_gate"]),
                int(pop_pair["pair_population_gate"]),
                int(pop_low["population_gate"]),
                int(pop_high["population_gate"]),
            ]
        )
        summary.append(
            {
                "n": n,
                "k": k,
                "one_minus_k": partner_k,
                "clone_count": clone_count,
                "horizon": horizon,
                "psi_k": float(low["mean_scgf"]),
                "psi_k_ci_low": float(low["scgf_ci_low"]),
                "psi_k_ci_high": float(low["scgf_ci_high"]),
                "psi_one_minus_k": float(high["mean_scgf"]),
                "psi_one_minus_k_ci_low": float(high["scgf_ci_low"]),
                "psi_one_minus_k_ci_high": float(high["scgf_ci_high"]),
                "gc_residual": float(pair["residual"]),
                "gc_residual_ci_low": float(pair["residual_ci_low"]),
                "gc_residual_ci_high": float(pair["residual_ci_high"]),
                "late_half_residual": float(pair["late_half_residual"]),
                "late_half_ci_low": float(pair["late_half_ci_low"]),
                "late_half_ci_high": float(pair["late_half_ci_high"]),
                "minimum_weight_ess_k": float(low["minimum_weight_ess"]),
                "minimum_weight_ess_one_minus_k": float(
                    high["minimum_weight_ess"]
                ),
                "minimum_final_roots_k": int(low["minimum_final_unique_roots"]),
                "minimum_final_roots_one_minus_k": int(
                    high["minimum_final_unique_roots"]
                ),
                "support_gate": int(
                    int(low["support_gate"]) and int(high["support_gate"])
                ),
                "paired_gc_gate": int(pair["paired_gc_gate"]),
                "time_convergence_gate": int(time["time_convergence_gate"]),
                "population_member_gates": int(
                    int(pop_low["population_gate"])
                    and int(pop_high["population_gate"])
                ),
                "population_pair_gate": int(pop_pair["pair_population_gate"]),
                "primary_gate": int(gate),
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "final_gc_summary.csv", summary)
    audit = {
        "pre_specified_rows": len(specs),
        "compiled_rows": len(summary),
        "passing_primary_rows": sum(int(row["primary_gate"]) for row in summary),
        "all_primary_gates_pass": all(int(row["primary_gate"]) for row in summary),
    }
    (args.output_dir / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    tex_rows = []
    for row in summary:
        tex_rows.append(
            f"{int(row['n'])} & "
            f"$({float(row['k']):.1f},{float(row['one_minus_k']):.1f})$ & "
            f"{int(row['clone_count'])} & {float(row['horizon']):.0f} & "
            f"{float(row['psi_k']):.5f} & "
            f"{float(row['psi_one_minus_k']):.5f} & "
            f"{float(row['gc_residual']):+.5f} & "
            f"$[{float(row['gc_residual_ci_low']):+.5f},"
            f"{float(row['gc_residual_ci_high']):+.5f}]$ & "
            f"{'pass' if int(row['primary_gate']) else 'unresolved'} \\\\"  # noqa: W605
        )
    (args.output_dir / "final_gc_rows.tex").write_text(
        "\n".join(tex_rows) + "\n\\bottomrule\n"
    )

    ns = np.asarray([int(row["n"]) for row in summary])
    x = np.arange(len(summary), dtype=float)
    residual = np.asarray([float(row["gc_residual"]) for row in summary])
    residual_low = np.asarray(
        [float(row["gc_residual_ci_low"]) for row in summary]
    )
    residual_high = np.asarray(
        [float(row["gc_residual_ci_high"]) for row in summary]
    )
    psi_low = np.asarray([float(row["psi_k"]) for row in summary])
    psi_high = np.asarray([float(row["psi_one_minus_k"]) for row in summary])
    psi_low_lo = np.asarray([float(row["psi_k_ci_low"]) for row in summary])
    psi_low_hi = np.asarray([float(row["psi_k_ci_high"]) for row in summary])
    psi_high_lo = np.asarray(
        [float(row["psi_one_minus_k_ci_low"]) for row in summary]
    )
    psi_high_hi = np.asarray(
        [float(row["psi_one_minus_k_ci_high"]) for row in summary]
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3))
    axes[0].errorbar(
        x - 0.06, psi_low,
        yerr=np.vstack((psi_low - psi_low_lo, psi_low_hi - psi_low)),
        fmt="o", capsize=4, label=r"$\psi_t(k)$",
    )
    axes[0].errorbar(
        x + 0.06, psi_high,
        yerr=np.vstack((psi_high - psi_high_lo, psi_high_hi - psi_high)),
        fmt="s", capsize=4, label=r"$\psi_t(1-k)$",
    )
    axes[0].set_ylabel("finite-time SCGF")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)

    axes[1].errorbar(
        x, residual,
        yerr=np.vstack((residual - residual_low, residual_high - residual)),
        fmt="o", capsize=4, color="tab:purple",
    )
    axes[1].axhline(0.0, color="black", linestyle="--", linewidth=1.2)
    axes[1].set_ylabel(r"$\psi_t(k)-\psi_t(1-k)$")
    axes[1].grid(alpha=0.25)

    tick_labels = [
        f"n={n}{'' if int(row['primary_gate']) else '*'}\n"
        f"({float(row['k']):.1f},{float(row['one_minus_k']):.1f})"
        for n, row in zip(ns, summary)
    ]
    for axis in axes:
        axis.set_xticks(x, tick_labels)
        axis.set_xlabel("chain length and tested tilt pair")
    fig.tight_layout()
    fig.savefig(args.output_dir / "long_chain_gc_summary.png", dpi=220)
    plt.close(fig)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
