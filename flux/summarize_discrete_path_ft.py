#!/usr/bin/env python3
"""Summarize timestep convergence of discrete path entropy toward heat entropy."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("audits", type=Path, nargs="+")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_single_row(path: Path) -> dict[str, str]:
    with path.open() as stream:
        return next(csv.DictReader(stream))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for audit_path in args.audits:
        audit = json.loads(audit_path.read_text())
        forward_path = Path(audit["forward_file"])
        stem = str(forward_path)
        if not stem.endswith("_forward.csv"):
            raise ValueError(f"unexpected forward path: {forward_path}")
        summary_path = Path(stem.removesuffix("_forward.csv") +
                            "_summary.csv")
        metadata = read_single_row(summary_path)
        for direction in ("forward", "reverse"):
            numerical = audit["numerical"][direction]
            ift = audit[f"{direction}_ift"]
            rows.append({
                "label": summary_path.stem.removesuffix("_summary"),
                "direction": direction,
                "dt": float(metadata["dt"]),
                "horizon": float(metadata["horizon"]),
                "trajectories": int(metadata["trajectories_per_direction"]),
                "mean_kernel_minus_heat":
                    numerical["mean_kernel_minus_heat"],
                "rms_kernel_minus_heat":
                    numerical["rms_kernel_minus_heat"],
                "rms_energy_balance_error":
                    numerical["rms_energy_balance_error"],
                "log_ift": ift["log_mean_exp_minus_sigma"],
                "log_ift_ci_low": ift["log_mean_exp_ci95"][0],
                "log_ift_ci_high": ift["log_mean_exp_ci95"][1],
                "ift_gate": audit["gates"][f"{direction}_ift"],
                "overall_gate": audit["gates"]["overall"],
            })
    rows.sort(key=lambda row: (float(row["dt"]), str(row["direction"])))
    with (args.output_dir / "timestep_convergence.csv").open(
            "w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    slopes: dict[str, dict[str, float]] = {}
    figure, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
    for direction, marker in (("forward", "o"), ("reverse", "s")):
        selected = [row for row in rows if row["direction"] == direction]
        dt = np.asarray([row["dt"] for row in selected], dtype=float)
        mean = np.abs(np.asarray(
            [row["mean_kernel_minus_heat"] for row in selected],
            dtype=float))
        rms = np.asarray([row["rms_kernel_minus_heat"] for row in selected],
                         dtype=float)
        order = np.argsort(dt)
        dt, mean, rms = dt[order], mean[order], rms[order]
        mean_slope, mean_intercept = np.polyfit(np.log(dt), np.log(mean), 1)
        rms_slope, rms_intercept = np.polyfit(np.log(dt), np.log(rms), 1)
        slopes[direction] = {
            "absolute_mean_order": float(mean_slope),
            "rms_order": float(rms_slope),
        }
        axes[0].loglog(dt, mean, marker=marker, label=direction)
        axes[0].loglog(dt, np.exp(mean_intercept) * dt ** mean_slope,
                       linestyle="--", alpha=0.7,
                       label=f"{direction} fit: dt^{mean_slope:.2f}")
        axes[1].loglog(dt, rms, marker=marker, label=direction)
        axes[1].loglog(dt, np.exp(rms_intercept) * dt ** rms_slope,
                       linestyle="--", alpha=0.7,
                       label=f"{direction} fit: dt^{rms_slope:.2f}")
    axes[0].set_xlabel(r"timestep $\Delta t$")
    axes[0].set_ylabel(r"$|\mathbb{E}[\Sigma_{K}-\Sigma_Q]|$")
    axes[1].set_xlabel(r"timestep $\Delta t$")
    axes[1].set_ylabel(r"RMS$(\Sigma_{K}-\Sigma_Q)$")
    for axis in axes:
        axis.legend(frameon=False, fontsize=8)
        axis.grid(alpha=0.2, which="both")
    figure.tight_layout()
    figure.savefig(args.output_dir / "kernel_heat_timestep_convergence.png",
                   dpi=220)
    figure.savefig(args.output_dir / "kernel_heat_timestep_convergence.pdf")
    plt.close(figure)

    result = {
        "points": len(rows),
        "orders": slopes,
        "interpretation": (
            "The exact discrete transition-kernel entropy approaches the "
            "bath-energy heat entropy with the fitted timestep orders."
        ),
    }
    (args.output_dir / "timestep_convergence.json").write_text(
        json.dumps(result, indent=2) + "\n")
    with (args.output_dir / "timestep_convergence.md").open("w") as stream:
        stream.write("# Kernel-entropy / heat-entropy timestep convergence\n\n")
        stream.write("| direction | absolute-mean order | RMS order |\n")
        stream.write("|---|---:|---:|\n")
        for direction in ("forward", "reverse"):
            stream.write(
                f"| {direction} | "
                f"{slopes[direction]['absolute_mean_order']:.4f} | "
                f"{slopes[direction]['rms_order']:.4f} |\n"
            )
        stream.write(
            "\nOrders near one support first-order convergence of the "
            "finite-step kernel entropy to the physical bath-heat entropy.\n"
        )
    print(f"Wrote timestep summary to {args.output_dir}")


if __name__ == "__main__":
    main()
