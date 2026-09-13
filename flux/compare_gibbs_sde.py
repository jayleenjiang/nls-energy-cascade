#!/usr/bin/env python3
"""Compare equal-temperature SDE action profiles with Gibbs MCMC."""

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcmc-profile", type=Path, required=True)
    parser.add_argument(
        "--sde-profiles", type=Path, nargs="+", required=True
    )
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def read_profile(path: Path) -> dict[str, np.ndarray]:
    raw = np.genfromtxt(path, delimiter=",", names=True)
    return {
        name: np.atleast_1d(raw[name]).astype(float)
        for name in raw.dtype.names or ()
    }


def main() -> None:
    args = parse_args()
    mcmc = read_profile(args.mcmc_profile)
    sde = [read_profile(path) for path in args.sde_profiles]
    modes = mcmc["mode"].astype(int)
    if not sde:
        raise ValueError("at least one SDE profile is required")
    for path, profile in zip(args.sde_profiles, sde):
        if not np.array_equal(profile["mode"].astype(int), modes):
            raise ValueError(f"{path}: mode grid does not match MCMC")

    sde_matrix = np.stack([profile["mean_action"] for profile in sde])
    sde_mean = np.mean(sde_matrix, axis=0)
    if len(sde) < 2:
        raise ValueError("at least two independent SDE runs are required")
    sde_se = np.std(sde_matrix, axis=0, ddof=1) / math.sqrt(len(sde))
    mcmc_mean = mcmc["mean_action"]
    mcmc_se = mcmc["between_chain_se"]
    difference = sde_mean - mcmc_mean
    combined_se = np.sqrt(sde_se**2 + mcmc_se**2)
    z_score = difference / combined_se

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    table_path = output_prefix.with_name(
        output_prefix.name + "_comparison.csv"
    )
    with table_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "gibbs_mcmc_mean",
                "gibbs_mcmc_se",
                "sde_mean",
                "sde_between_run_se",
                "difference",
                "combined_z",
            ]
        )
        for values in zip(
            modes,
            mcmc_mean,
            mcmc_se,
            sde_mean,
            sde_se,
            difference,
            z_score,
        ):
            writer.writerow(values)

    sde_reflection_difference = sde_mean - sde_mean[::-1]
    sde_reflection_se = np.sqrt(sde_se**2 + sde_se[::-1] ** 2)
    summary = {
        "mcmc_profile": str(args.mcmc_profile),
        "sde_profiles": [str(path) for path in args.sde_profiles],
        "independent_sde_runs": len(sde),
        "max_abs_combined_z": float(np.max(np.abs(z_score))),
        "rms_relative_profile_difference": float(
            np.sqrt(np.mean((difference / mcmc_mean) ** 2))
        ),
        "max_abs_relative_profile_difference": float(
            np.max(np.abs(difference / mcmc_mean))
        ),
        "sde_reflection_max_abs_z": float(
            np.max(np.abs(sde_reflection_difference / sde_reflection_se))
        ),
        "interpretation": (
            "Combined z uses the between-chain MCMC SE and the between-run "
            "SDE SE. Agreement of profile means is a moment-level check of the "
            "Gibbs invariant measure, not a proof of full distributional "
            "equality."
        ),
    }
    summary_path = output_prefix.with_name(
        output_prefix.name + "_summary.json"
    )
    with summary_path.open("w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    fig, axis = plt.subplots(figsize=(6.0, 4.2))
    axis.errorbar(
        modes,
        mcmc_mean,
        yerr=1.959963984540054 * mcmc_se,
        fmt="o-",
        capsize=3,
        label="Gibbs MCMC (95% CI)",
        color="#1f4e79",
    )
    axis.errorbar(
        modes,
        sde_mean,
        yerr=1.959963984540054 * sde_se,
        fmt="s--",
        capsize=3,
        label="equal-temperature SDE (95% CI)",
        color="#b03a2e",
    )
    axis.set_xlabel("mode")
    axis.set_ylabel("mean action")
    axis.grid(True, alpha=0.2)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(
        output_prefix.with_name(output_prefix.name + "_profile.pdf")
    )
    fig.savefig(
        output_prefix.with_name(output_prefix.name + "_profile.png"), dpi=240
    )

    print(json.dumps(summary, indent=2))
    print(f"wrote {table_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
