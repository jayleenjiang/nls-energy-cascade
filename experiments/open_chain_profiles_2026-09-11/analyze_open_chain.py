#!/usr/bin/env python3
"""Frozen analysis for the one-bath open-chain profile diagnostic."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


Z95 = 1.959963984540054
T975_DF1 = 12.706204736432095
NS = (25, 50, 100)
TEMPS = ("T10", "T6")
INTERIOR_LO = 0.10
INTERIOR_HI = 0.90
PENETRATION_FRACTION = 0.50
EXPECTED_SOURCE_SHA = "dfee7ed0b416a30f2b27d919c524dfaf37e91ca44312ffbb9d68b7fc1cd36110"
EXPECTED_BINARY_SHA = "09d801e71b7f3bd13cea56aba3d432b471dbf16a726a9866894148a91397575d"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(lines))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"empty rows: {path}")
        fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def f(value) -> float:
    if value in (None, ""):
        return math.nan
    return float(value)


def combine(means: list[float], ses: list[float], counts: list[int]) -> tuple[float, float, int]:
    keep = [
        i for i, (mean, se, count) in enumerate(zip(means, ses, counts))
        if count > 0 and math.isfinite(mean) and math.isfinite(se)
    ]
    if not keep:
        return math.nan, math.nan, 0
    total = sum(counts[i] for i in keep)
    mean = sum(counts[i] * means[i] for i in keep) / total
    if total <= 1:
        return mean, math.nan, total
    ss = 0.0
    for i in keep:
        variance = ses[i] ** 2 * counts[i]
        ss += (counts[i] - 1) * variance + counts[i] * (means[i] - mean) ** 2
    variance = ss / (total - 1)
    return mean, math.sqrt(variance / total), total


def merge_rows(rows_by_rep: list[list[dict[str, str]]], keys: tuple[str, ...],
               pairs: tuple[tuple[str, str, str], ...], passthrough: tuple[str, ...]) -> list[dict]:
    indexed = []
    for rows in rows_by_rep:
        indexed.append({tuple(row[key] for key in keys): row for row in rows})
    if set(indexed[0]) != set(indexed[1]):
        raise RuntimeError(f"replicate keys differ for {keys}")
    merged = []
    for key in sorted(indexed[0], key=lambda item: tuple(float(x) if x.replace('.', '', 1).isdigit() else x for x in item)):
        reps = [table[key] for table in indexed]
        out = {name: value for name, value in zip(keys, key)}
        for name in passthrough:
            out[name] = reps[0][name]
        for mean_field, se_field, count_field in pairs:
            counts = [int(rep[count_field]) for rep in reps]
            mean, se, count = combine(
                [f(rep[mean_field]) for rep in reps],
                [f(rep[se_field]) for rep in reps], counts,
            )
            out[mean_field] = mean
            out[se_field] = se
            if count > 0 or count_field not in out:
                out[count_field] = count
        merged.append(out)
    return merged


def merge_profile(paths: list[Path]) -> list[dict]:
    return merge_rows(
        [read_csv(path) for path in paths], ("j",),
        (("mean_I", "se_mean_I", "sample_count"),
         ("mean_sin_theta", "se_mean_sin_theta", "sample_count")), (),
    )


def merge_quarters(paths: list[Path]) -> list[dict]:
    return merge_rows(
        [read_csv(path) for path in paths], ("quarter", "j"),
        (("mean_I", "se_mean_I", "sample_count"),
         ("mean_sin_theta", "se_mean_sin_theta", "sample_count")),
        ("interval_start", "interval_end"),
    )


def merge_state(paths: list[Path]) -> list[dict]:
    return merge_rows(
        [read_csv(path) for path in paths], ("phase", "checkpoint", "j"),
        (("mean_I", "se_mean_I", "sample_count"),
         ("mean_sin_theta", "se_mean_sin_theta", "sample_count")),
        ("absolute_time",),
    )


def merge_mass(paths: list[Path]) -> list[dict]:
    return merge_rows(
        [read_csv(path) for path in paths], ("series", "checkpoint"),
        (("mean_total_action", "se_total_action", "sample_count"),
         ("mean_I_left", "se_I_left", "sample_count"),
         ("mean_I_right", "se_I_right", "sample_count")),
        ("absolute_time",),
    )


def merge_differences(paths: list[Path]) -> list[dict]:
    reps = [read_csv(path) for path in paths]
    indexed = [{row["j"]: row for row in rows} for rows in reps]
    result = []
    for j in sorted(indexed[0], key=int):
        rows = [table[j] for table in indexed]
        counts = [int(row["sample_count"]) for row in rows]
        q3, q3se, count = combine(
            [f(row["mean_I_q3"]) for row in rows],
            [0.0 for _ in rows], counts,
        )
        q4, q4se, _ = combine(
            [f(row["mean_I_q4"]) for row in rows],
            [0.0 for _ in rows], counts,
        )
        delta, delta_se, _ = combine(
            [f(row["mean_delta_I_q4_minus_q3"]) for row in rows],
            [f(row["se_delta_I"]) for row in rows], counts,
        )
        denominator = 0.5 * (abs(q3) + abs(q4))
        result.append({
            "j": int(j), "mean_I_q3": q3, "mean_I_q4": q4,
            "mean_delta_I_q4_minus_q3": delta, "se_delta_I": delta_se,
            "relative_delta_I": delta / denominator if denominator else math.nan,
            "abs_relative_delta_I": abs(delta) / denominator if denominator else math.nan,
            "z_delta_I": delta / delta_se if delta_se > 0 else math.nan,
            "sample_count": count,
        })
    return result


def fit_power(ns: np.ndarray, values: np.ndarray, ses: np.ndarray) -> dict:
    x = np.log(ns)
    y = np.log(values)
    design = np.column_stack((np.ones_like(x), x))
    beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    s2 = float(residual @ residual)  # one residual degree of freedom
    covariance = s2 * np.linalg.inv(design.T @ design)
    slope_se = math.sqrt(max(0.0, float(covariance[1, 1])))
    return {
        "amplitude": math.exp(float(beta[0])), "exponent": float(beta[1]),
        "exponent_se_ols": slope_se,
        "ci_low_df1": float(beta[1]) - T975_DF1 * slope_se,
        "ci_high_df1": float(beta[1]) + T975_DF1 * slope_se,
        "log_sse": s2, "residual_dof": 1,
        "relative_measurement_se_max": float(np.max(ses / values)),
    }


def plot_profiles(figures: Path, logical: dict, temp: str, rescaled: bool) -> None:
    colors = {25: "#0072B2", 50: "#E69F00", 100: "#009E73"}
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for n in NS:
        rows = logical[(temp, n)]["profile"]
        x = np.asarray([(int(row["j"]) - 1) / (n - 1) for row in rows])
        mean = np.asarray([f(row["mean_I"]) for row in rows])
        se = np.asarray([f(row["se_mean_I"]) for row in rows])
        factor = math.sqrt(n) if rescaled else 1.0
        ax.plot(x, factor * mean, color=colors[n], lw=1.7, label=fr"$n={n}$")
        ax.fill_between(x, factor * (mean - Z95 * se), factor * (mean + Z95 * se),
                        color=colors[n], alpha=0.16, linewidth=0)
    ax.set_xlabel(r"normalised site $x=(j-1)/(n-1)$")
    ax.set_ylabel(r"$\sqrt{n}\,\langle I_j\rangle$" if rescaled else r"$\langle I_j\rangle$")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    stem = f"open_action_profiles_{'rescaled_' if rescaled else ''}{temp}"
    fig.savefig(figures / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figures / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_quarters(figures: Path, logical: dict, temp: str) -> None:
    colors = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 10.2), sharex=True)
    for ax, n in zip(axes, NS):
        rows = logical[(temp, n)]["quarters"]
        for quarter in range(1, 5):
            selected = [row for row in rows if int(row["quarter"]) == quarter]
            x = [(int(row["j"]) - 1) / (n - 1) for row in selected]
            y = [f(row["mean_I"]) for row in selected]
            se = [f(row["se_mean_I"]) for row in selected]
            ax.plot(x, y, color=colors[quarter - 1], lw=1.4, label=f"Q{quarter}")
            ax.fill_between(x, np.asarray(y) - Z95 * np.asarray(se),
                            np.asarray(y) + Z95 * np.asarray(se),
                            color=colors[quarter - 1], alpha=0.10, linewidth=0)
        ax.set_ylabel(fr"$n={n}$")
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, ncol=4)
    axes[-1].set_xlabel(r"normalised site $x$")
    fig.supylabel(r"quarter-averaged $\langle I_j\rangle$")
    fig.tight_layout()
    fig.savefig(figures / f"open_quarter_profiles_{temp}.pdf", bbox_inches="tight")
    fig.savefig(figures / f"open_quarter_profiles_{temp}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_mass(figures: Path, logical: dict, temp: str) -> None:
    colors = {25: "#0072B2", 50: "#E69F00", 100: "#009E73"}
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3))
    for n in NS:
        rows = [row for row in logical[(temp, n)]["mass"]
                if row["series"] in {"burnin_snapshot", "measurement_snapshot"}]
        rows.sort(key=lambda row: f(row["absolute_time"]))
        time = np.asarray([f(row["absolute_time"]) for row in rows])
        for ax, field, sefield in (
            (axes[0], "mean_total_action", "se_total_action"),
            (axes[1], "mean_I_right", "se_I_right"),
        ):
            y = np.asarray([f(row[field]) for row in rows])
            se = np.asarray([f(row[sefield]) for row in rows])
            ax.errorbar(time, y, yerr=Z95 * se, color=colors[n], marker="o",
                        markersize=3, capsize=2, lw=1.2, label=fr"$n={n}$")
    axes[0].set_ylabel(r"instantaneous ensemble mean $M$")
    axes[1].set_ylabel(r"instantaneous free-end $\langle I_n\rangle$")
    for ax in axes:
        ax.set_xlabel("simulation time")
        ax.grid(alpha=0.22)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures / f"open_stationarity_{temp}.pdf", bbox_inches="tight")
    fig.savefig(figures / f"open_stationarity_{temp}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_dir.resolve()
    raw = root / "raw"
    curated = root / "curated_results"
    analysis = root / "analysis"
    figures = root / "figures"
    for directory in (curated, analysis, figures):
        directory.mkdir(parents=True, exist_ok=True)

    source = root / "src/NLS_open_chain_profiles.cpp"
    binary = root / "bin/NLS_open_chain_profiles"
    if sha256(source) != EXPECTED_SOURCE_SHA or sha256(binary) != EXPECTED_BINARY_SHA:
        raise RuntimeError("source or binary hash differs from frozen identity")

    matrix = list(csv.DictReader((root / "RUN_MATRIX.csv").open()))
    groups: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for row in matrix:
        groups[(row["temperature_label"], int(row["n"]))].append(row)
    logical = {}
    run_rows = []
    input_hash_rows = []
    for key, runs in sorted(groups.items()):
        temp, n = key
        runs.sort(key=lambda row: int(row["replicate"]))
        prefixes = [raw / f"OPEN_{temp}_n{n}_rep{row['replicate']}" for row in runs]
        for prefix, row in zip(prefixes, runs):
            required = (
                "_profile.csv", "_quarter_profiles.csv", "_quarter_differences.csv",
                "_burnin_checkpoints.csv", "_checkpoints.csv", "_state_checkpoints.csv",
                "_mass_checkpoints.csv", "_trajectory_diagnostics.csv", "_summary.csv",
            )
            for suffix in required:
                path = Path(str(prefix) + suffix)
                if not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError(f"missing output {path}")
                input_hash_rows.append({"path": str(path), "sha256": sha256(path)})
            summary = read_csv(Path(str(prefix) + "_summary.csv"))[0]
            if int(summary["valid_trajectories"]) <= 0:
                raise RuntimeError(f"no valid trajectories: {prefix}")
            run_rows.append(summary)

        profile = merge_profile([Path(str(prefix) + "_profile.csv") for prefix in prefixes])
        quarters = merge_quarters([Path(str(prefix) + "_quarter_profiles.csv") for prefix in prefixes])
        differences = merge_differences([Path(str(prefix) + "_quarter_differences.csv") for prefix in prefixes])
        state = merge_state([Path(str(prefix) + "_state_checkpoints.csv") for prefix in prefixes])
        mass = merge_mass([Path(str(prefix) + "_mass_checkpoints.csv") for prefix in prefixes])
        logical[key] = {"profile": profile, "quarters": quarters,
                        "differences": differences, "state": state, "mass": mass,
                        "runs": runs}
        stem = curated / f"OPEN_{temp}_n{n}"
        write_csv(Path(str(stem) + "_profile.csv"), profile)
        write_csv(Path(str(stem) + "_quarter_profiles.csv"), quarters)
        write_csv(Path(str(stem) + "_quarter_differences.csv"), differences)
        write_csv(Path(str(stem) + "_state_checkpoints.csv"), state)
        write_csv(Path(str(stem) + "_mass_checkpoints.csv"), mass)

    stationarity_rows = []
    penetration_rows = []
    scaling_points = []
    curves = []
    spread_rows = []
    for temp in TEMPS:
        for n in NS:
            item = logical[(temp, n)]
            profile = item["profile"]
            differences = item["differences"]
            mass = item["mass"]
            qmass = {int(row["checkpoint"]): row for row in mass
                     if row["series"] == "measurement_quarter"}
            delta = next(row for row in mass
                         if row["series"] == "measurement_quarter_difference_q4_minus_q3")
            mass_den = 0.5 * (abs(f(qmass[3]["mean_total_action"]))
                              + abs(f(qmass[4]["mean_total_action"])))
            right_den = 0.5 * (abs(f(qmass[3]["mean_I_right"]))
                               + abs(f(qmass[4]["mean_I_right"])))
            mass_rel = abs(f(delta["mean_total_action"])) / mass_den
            right_rel = abs(f(delta["mean_I_right"])) / right_den
            mass_ci_zero = abs(f(delta["mean_total_action"])) <= Z95 * f(delta["se_total_action"])
            right_ci_zero = abs(f(delta["mean_I_right"])) <= Z95 * f(delta["se_I_right"])
            interior = [row for row in differences
                        if INTERIOR_LO <= (int(row["j"]) - 1) / (n - 1) <= INTERIOR_HI]
            rels = np.asarray([f(row["abs_relative_delta_I"]) for row in interior])
            median_rel, max_rel = float(np.median(rels)), float(np.max(rels))
            profile_pass = median_rel <= 0.05 and max_rel <= 0.10
            stationary = mass_ci_zero and mass_rel <= 0.05 and right_ci_zero \
                and right_rel <= 0.05 and profile_pass
            stationarity_rows.append({
                "temperature": temp, "n": n,
                "delta_mass_q4_minus_q3": f(delta["mean_total_action"]),
                "se_delta_mass": f(delta["se_total_action"]),
                "z_delta_mass": f(delta["mean_total_action"]) / f(delta["se_total_action"]),
                "mass_ci_contains_zero": int(mass_ci_zero),
                "mass_relative_change": mass_rel,
                "delta_free_end_q4_minus_q3": f(delta["mean_I_right"]),
                "se_delta_free_end": f(delta["se_I_right"]),
                "z_delta_free_end": f(delta["mean_I_right"]) / f(delta["se_I_right"]),
                "free_end_ci_contains_zero": int(right_ci_zero),
                "free_end_relative_change": right_rel,
                "free_end_significantly_growing": int(
                    f(delta["mean_I_right"]) - Z95 * f(delta["se_I_right"]) > 0),
                "interior_median_site_relative_change": median_rel,
                "interior_max_site_relative_change": max_rel,
                "profile_change_gate_pass": int(profile_pass),
                "stationary_gate_pass": int(stationary),
            })

            x = np.asarray([(int(row["j"]) - 1) / (n - 1) for row in profile])
            mean = np.asarray([f(row["mean_I"]) for row in profile])
            se = np.asarray([f(row["se_mean_I"]) for row in profile])
            mid = (n + 1) // 2
            scaling_points.extend((
                {"temperature": temp, "metric": "mid_chain", "n": n,
                 "j": mid, "mean_I": mean[mid - 1], "se_mean_I": se[mid - 1],
                 "stationary": int(stationary)},
                {"temperature": temp, "metric": "free_end", "n": n,
                 "j": n, "mean_I": mean[-1], "se_mean_I": se[-1],
                 "stationary": int(stationary)},
            ))

            threshold = PENETRATION_FRACTION * mean[0]
            candidates = np.flatnonzero(mean <= threshold)
            if candidates.size:
                idx = int(candidates[0])
                if idx == 0:
                    crossing_j = 1.0
                else:
                    y0, y1 = mean[idx - 1], mean[idx]
                    fraction = (threshold - y0) / (y1 - y0) if y1 != y0 else 1.0
                    crossing_j = idx + fraction + 1.0
                first_site = idx + 1
                reached = 1
                crossing_fraction = crossing_j / n
            else:
                first_site = ""
                crossing_j = crossing_fraction = math.nan
                reached = 0
            penetration_rows.append({
                "temperature": temp, "n": n, "fraction_of_left": PENETRATION_FRACTION,
                "left_action": mean[0], "threshold": threshold,
                "threshold_reached": reached, "first_site_at_or_below": first_site,
                "interpolated_site": crossing_j,
                "interpolated_fraction_of_n": crossing_fraction,
                "stationary": int(stationary),
            })

        grid = np.linspace(0.0, 1.0, 101)
        raw_matrix, scaled_matrix = [], []
        for n in NS:
            profile = logical[(temp, n)]["profile"]
            x = np.asarray([(int(row["j"]) - 1) / (n - 1) for row in profile])
            mean = np.asarray([f(row["mean_I"]) for row in profile])
            se = np.asarray([f(row["se_mean_I"]) for row in profile])
            interp = np.interp(grid, x, mean)
            interp_se = np.interp(grid, x, se)
            raw_matrix.append(interp)
            scaled_matrix.append(math.sqrt(n) * interp)
            for gx, value, error in zip(grid, interp, interp_se):
                curves.append({"temperature": temp, "n": n, "x": gx,
                               "mean_I_interpolated": value, "se_I_interpolated": error,
                               "sqrt_n_mean_I": math.sqrt(n) * value,
                               "se_sqrt_n_mean_I": math.sqrt(n) * error})
        for scale, matrix in (("raw", np.asarray(raw_matrix)),
                              ("sqrt_n", np.asarray(scaled_matrix))):
            spread = (np.max(matrix, axis=0) - np.min(matrix, axis=0)) / np.mean(matrix, axis=0)
            for region, mask in (("full", np.ones(grid.size, dtype=bool)),
                                 ("interior_0p10_0p90", (grid >= .10) & (grid <= .90))):
                spread_rows.append({"temperature": temp, "scale": scale, "region": region,
                                    "grid_points": int(np.count_nonzero(mask)),
                                    "median_relative_spread": float(np.median(spread[mask])),
                                    "max_relative_spread": float(np.max(spread[mask]))})

        plot_profiles(figures, logical, temp, False)
        plot_profiles(figures, logical, temp, True)
        plot_quarters(figures, logical, temp)
        plot_mass(figures, logical, temp)

    scaling_rows = []
    for temp in TEMPS:
        for metric in ("mid_chain", "free_end"):
            points = [row for row in scaling_points
                      if row["temperature"] == temp and row["metric"] == metric]
            points.sort(key=lambda row: row["n"])
            fit = fit_power(
                np.asarray([row["n"] for row in points], dtype=float),
                np.asarray([row["mean_I"] for row in points], dtype=float),
                np.asarray([row["se_mean_I"] for row in points], dtype=float),
            )
            fit.update({"temperature": temp, "metric": metric,
                        "all_stationary": int(all(row["stationary"] for row in points)),
                        "n_values": "25;50;100"})
            scaling_rows.append(fit)

    write_csv(analysis / "run_summary.csv", run_rows)
    write_csv(analysis / "input_hashes.csv", input_hash_rows)
    write_csv(analysis / "stationarity_gates.csv", stationarity_rows)
    write_csv(analysis / "scaling_points.csv", scaling_points)
    write_csv(analysis / "scaling_fits.csv", scaling_rows)
    write_csv(analysis / "penetration_depth.csv", penetration_rows)
    write_csv(analysis / "rescaled_curves.csv", curves)
    write_csv(analysis / "collapse_spread.csv", spread_rows)
    all_stationary = all(row["stationary_gate_pass"] for row in stationarity_rows)
    verdict = {
        "all_six_logical_conditions_stationary": bool(all_stationary),
        "stationary_conditions": [f"{r['temperature']}_n{r['n']}" for r in stationarity_rows
                                  if r["stationary_gate_pass"]],
        "nonstationary_conditions": [f"{r['temperature']}_n{r['n']}" for r in stationarity_rows
                                     if not r["stationary_gate_pass"]],
        "profile_claim_boundary": (
            "stationary_profile_scaling_allowed" if all_stationary
            else "scaling_and_collapse_descriptive_for_nonstationary_conditions"
        ),
        "bond_sine_quantitative_use": False,
        "penetration_fraction": PENETRATION_FRACTION,
    }
    (analysis / "VERDICT.json").write_text(json.dumps(verdict, indent=2) + "\n")
    manifest = {
        "analysis_sha256": sha256(Path(__file__)),
        "source_sha256": sha256(source), "binary_sha256": sha256(binary),
        "protocol_sha256": sha256(root / "PROTOCOL.md"),
        "run_matrix_sha256": sha256(root / "RUN_MATRIX.csv"),
        "grid": "0:0.01:1", "interior": [INTERIOR_LO, INTERIOR_HI],
        "stationarity_relative_thresholds": {"mass": .05, "free_end": .05,
                                                "profile_median": .05, "profile_max": .10},
        "inputs": input_hash_rows,
    }
    (analysis / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
