#!/usr/bin/env python3
"""Frozen analysis for the BC1 canonical equal-temperature diagnostic."""

import csv
import hashlib
import math
from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
OUT = ROOT / "analysis"
FIG = ROOT / "figures"
SIMD_ROOT = ROOT.parent / "boundary_profiles_2026-09-05" / "curated_results" / "profiles"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)


def mean_se(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(values.size)), int(values.size)


def load_trajectory_profile(n, rep):
    path = RAW / f"canonical_T6_T6_n{n}_rep{rep}_trajectory_profiles.csv"
    action = np.full((256, n), np.nan)
    sine = np.full((256, n - 1), np.nan)
    with path.open() as handle:
        for row in csv.DictReader(handle):
            trajectory = int(row["trajectory_id"])
            j = int(row["j"]) - 1
            action[trajectory, j] = float(row["mean_I"])
            if j < n - 1:
                sine[trajectory, j] = float(row["mean_sin_theta"])
    if not np.isfinite(action).all() or not np.isfinite(sine).all():
        raise RuntimeError(f"non-finite or missing trajectory profile in {path}")
    return action, sine


def load_simd(n):
    path = SIMD_ROOT / f"BC1_T6_T6_n{n}_profile.csv"
    with path.open() as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    return {
        "mean_I": np.array([float(row["mean_I"]) for row in rows]),
        "se_I": np.array([float(row["se_mean_I"]) for row in rows]),
        "mean_sin": np.array([float(row["mean_sin_theta"]) for row in rows[:-1]]),
        "se_sin": np.array([float(row["se_mean_sin_theta"]) for row in rows[:-1]]),
    }


def load_run_summary(n, rep):
    path = RAW / f"canonical_T6_T6_n{n}_rep{rep}_summary.csv"
    with path.open() as handle:
        return next(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    profile_rows = []
    comparison_rows = []
    seed_rows = []
    plot_data = {}

    for n in (25, 100):
        actions = []
        sines = []
        summaries = []
        for rep in (0, 1):
            action, sine = load_trajectory_profile(n, rep)
            actions.append(action)
            sines.append(sine)
            summary = load_run_summary(n, rep)
            summaries.append(summary)
            sine_mean = sine.mean(axis=0)
            sine_se = sine.std(axis=0, ddof=1) / math.sqrt(sine.shape[0])
            z = np.abs(sine_mean / sine_se)
            seed_rows.append({
                "n": n,
                "replicate": rep,
                "seed": summary["seed"],
                "n_trajectories": sine.shape[0],
                "spatial_mean_I": action.mean(),
                "max_abs_z_sin": z.max(),
                "pointwise_95_failures": int((z > NormalDist().inv_cdf(0.975)).sum()),
                "projection_count": summary["projection_count"],
                "near_floor_count": summary["near_floor_count"],
                "minimum_proposed_action": summary["minimum_proposed_action"],
                "elapsed_seconds": summary["elapsed_seconds"],
            })

        action = np.concatenate(actions, axis=0)
        sine = np.concatenate(sines, axis=0)
        simd = load_simd(n)
        mean_I = action.mean(axis=0)
        se_I = action.std(axis=0, ddof=1) / math.sqrt(action.shape[0])
        mean_sin = sine.mean(axis=0)
        se_sin = sine.std(axis=0, ddof=1) / math.sqrt(sine.shape[0])
        z_sin = mean_sin / se_sin
        point_critical = NormalDist().inv_cdf(0.975)
        bonferroni_critical = NormalDist().inv_cdf(1.0 - 0.05 / (2.0 * (n - 1)))

        for j in range(n):
            profile_rows.append({
                "n": n,
                "j": j + 1,
                "canonical_mean_I": mean_I[j],
                "canonical_se_mean_I": se_I[j],
                "simd_mean_I": simd["mean_I"][j],
                "simd_se_mean_I": simd["se_I"][j],
                "canonical_mean_sin_theta": mean_sin[j] if j < n - 1 else "",
                "canonical_se_mean_sin_theta": se_sin[j] if j < n - 1 else "",
                "canonical_z_sin": z_sin[j] if j < n - 1 else "",
                "canonical_pointwise_95_fail": bool(abs(z_sin[j]) > point_critical) if j < n - 1 else "",
                "simd_mean_sin_theta": simd["mean_sin"][j] if j < n - 1 else "",
                "simd_se_mean_sin_theta": simd["se_sin"][j] if j < n - 1 else "",
            })

        simd_z = np.abs(simd["mean_sin"] / simd["se_sin"])
        canonical_spatial_by_trajectory = action.mean(axis=1)
        canonical_spatial_mean, canonical_spatial_se, _ = mean_se(canonical_spatial_by_trajectory)
        simd_spatial_mean = float(simd["mean_I"].mean())
        relative_mean_difference = abs(canonical_spatial_mean - simd_spatial_mean) / simd_spatial_mean
        rms_profile_relative_difference = float(
            np.sqrt(np.mean((mean_I - simd["mean_I"]) ** 2)) / simd_spatial_mean
        )
        comparison_rows.append({
            "n": n,
            "canonical_max_abs_z_sin": float(np.abs(z_sin).max()),
            "canonical_pointwise_95_failures": int((np.abs(z_sin) > point_critical).sum()),
            "number_of_bonds": n - 1,
            "bonferroni_critical": bonferroni_critical,
            "canonical_bonferroni_pass": bool(np.abs(z_sin).max() <= bonferroni_critical),
            "simd_max_abs_z_sin": float(simd_z.max()),
            "simd_pointwise_95_failures": int((simd_z > point_critical).sum()),
            "simd_bonferroni_pass": bool(simd_z.max() <= bonferroni_critical),
            "canonical_spatial_mean_I": canonical_spatial_mean,
            "canonical_se_spatial_mean_I": canonical_spatial_se,
            "simd_spatial_mean_I": simd_spatial_mean,
            "relative_spatial_mean_difference": relative_mean_difference,
            "action_mean_within_five_percent": bool(relative_mean_difference <= 0.05),
            "rms_action_profile_relative_difference": rms_profile_relative_difference,
            "projection_count": sum(int(row["projection_count"]) for row in summaries),
            "near_floor_count": sum(int(row["near_floor_count"]) for row in summaries),
            "minimum_proposed_action": min(float(row["minimum_proposed_action"]) for row in summaries),
        })
        plot_data[n] = (mean_I, se_I, mean_sin, se_sin, simd)

    write_csv(
        OUT / "canonical_vs_simd_profiles.csv",
        profile_rows,
        list(profile_rows[0]),
    )
    write_csv(OUT / "comparison_summary.csv", comparison_rows, list(comparison_rows[0]))
    write_csv(OUT / "per_seed_summary.csv", seed_rows, list(seed_rows[0]))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for axis, n in zip(axes, (25, 100)):
        _, _, mean_sin, se_sin, simd = plot_data[n]
        x = np.arange(1, n)
        axis.axhline(0.0, color="0.4", lw=0.8)
        axis.plot(x, simd["mean_sin"], color="#d95f02", lw=1.4, label="SIMD")
        axis.fill_between(x, simd["mean_sin"] - 1.96 * simd["se_sin"],
                          simd["mean_sin"] + 1.96 * simd["se_sin"],
                          color="#d95f02", alpha=0.18)
        axis.plot(x, mean_sin, color="#1b9e77", lw=1.4, label="canonical")
        axis.fill_between(x, mean_sin - 1.96 * se_sin, mean_sin + 1.96 * se_sin,
                          color="#1b9e77", alpha=0.18)
        axis.set_title(f"n={n}")
        axis.set_xlabel("bond j")
        axis.set_ylabel(r"$\langle\sin\theta_j\rangle$")
        axis.legend(frameon=False)
    fig.savefig(FIG / "bond_sine_canonical_vs_simd.png", dpi=220)
    fig.savefig(FIG / "bond_sine_canonical_vs_simd.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for axis, n in zip(axes, (25, 100)):
        mean_I, se_I, _, _, simd = plot_data[n]
        x = np.arange(1, n + 1)
        axis.plot(x, simd["mean_I"], color="#d95f02", lw=1.4, label="SIMD")
        axis.fill_between(x, simd["mean_I"] - 1.96 * simd["se_I"],
                          simd["mean_I"] + 1.96 * simd["se_I"],
                          color="#d95f02", alpha=0.18)
        axis.plot(x, mean_I, color="#1b9e77", lw=1.4, label="canonical")
        axis.fill_between(x, mean_I - 1.96 * se_I, mean_I + 1.96 * se_I,
                          color="#1b9e77", alpha=0.18)
        axis.set_title(f"n={n}")
        axis.set_xlabel("site j")
        axis.set_ylabel(r"$\langle I_j\rangle$")
        axis.legend(frameon=False)
    fig.savefig(FIG / "action_profile_canonical_vs_simd.png", dpi=220)
    fig.savefig(FIG / "action_profile_canonical_vs_simd.pdf")
    plt.close(fig)

    provenance = []
    for path in (
        ROOT / "source_archive" / "NLS_flux_canonical.original.cpp",
        ROOT / "src" / "NLS_flux_canonical_profile.cpp",
        ROOT / "bin" / "NLS_flux_canonical_profile",
    ):
        provenance.append({
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
        })
    write_csv(OUT / "provenance_hashes.csv", provenance, list(provenance[0]))

    for row in comparison_rows:
        print(row)


if __name__ == "__main__":
    main()

