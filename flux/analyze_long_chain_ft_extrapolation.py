#!/usr/bin/env python3
"""Finite-time and finite-population audit for controlled long-chain SCGFs.

The primary estimator is the intercept of

    psi_{t,N}(k) = psi_{infinity,N}(k) + a_N(k) / t

over a fixed pre-coalescence window.  Independent cloning runs are fitted
separately; uncertainty is computed across run-level intercepts.  The script
does not impose Gallavotti--Cohen symmetry in any fit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment_root", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--window-start", type=float, default=20.0)
    parser.add_argument("--window-stop", type=float, default=40.0)
    parser.add_argument("--window-step", type=float, default=5.0)
    return parser.parse_args()


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}")
    return rows[0]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty file {path}")
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mean_se_ci(values: list[float]) -> tuple[float, float, float, float]:
    array = np.asarray(values, dtype=float)
    mean = float(array.mean())
    if array.size < 2:
        return mean, math.nan, math.nan, math.nan
    se = float(array.std(ddof=1) / math.sqrt(array.size))
    critical = float(student_t.ppf(0.975, array.size - 1))
    return mean, se, mean - critical * se, mean + critical * se


def welch_interval(
    mean_a: float, se_a: float, n_a: int,
    mean_b: float, se_b: float, n_b: int,
) -> tuple[float, float, float, float]:
    delta = mean_a - mean_b
    va = se_a * se_a
    vb = se_b * se_b
    se = math.sqrt(va + vb)
    if not math.isfinite(se) or se == 0.0 or n_a < 2 or n_b < 2:
        return delta, se, math.nan, math.nan
    denominator = va * va / (n_a - 1) + vb * vb / (n_b - 1)
    degrees = (va + vb) ** 2 / denominator
    critical = float(student_t.ppf(0.975, degrees))
    return delta, se, delta - critical * se, delta + critical * se


def convergence_gate(delta: float, se_a: float, se_b: float) -> bool:
    combined = math.hypot(se_a, se_b)
    return (
        math.isfinite(delta) and math.isfinite(combined)
        and abs(delta) <= 0.01 and abs(delta) <= 2.0 * combined
    )


def fit_line(times: np.ndarray, values: np.ndarray) -> tuple[float, float, float, float]:
    x = 1.0 / times
    slope, intercept = np.polyfit(x, values, 1)
    fitted = slope * x + intercept
    residual = values - fitted
    rss = float(np.dot(residual, residual))
    centered = values - values.mean()
    tss = float(np.dot(centered, centered))
    r2 = 1.0 if tss <= 1.0e-30 else 1.0 - rss / tss
    rmse = math.sqrt(rss / values.size)
    return float(intercept), float(slope), r2, rmse


def main() -> None:
    args = parse_args()
    if not (
        args.window_start > 0.0
        and args.window_stop > args.window_start
        and args.window_step > 0.0
    ):
        raise SystemExit("invalid extrapolation window")
    count = int(round((args.window_stop - args.window_start) / args.window_step))
    times = args.window_start + args.window_step * np.arange(count + 1)
    if not math.isclose(float(times[-1]), args.window_stop):
        raise SystemExit("window must be divisible by window step")

    runs: list[dict[str, object]] = []
    curves: dict[tuple[object, ...], list[np.ndarray]] = defaultdict(list)
    for summary_path in sorted(args.experiment_root.rglob("*_summary.csv")):
        summary = read_one(summary_path)
        if not summary.get("mode", "").startswith("controlled_exact"):
            continue
        if float(summary["observation_time"]) + 1.0e-12 < args.window_stop:
            continue
        series_path = summary_path.with_name(
            summary_path.name.replace("_summary.csv", "_timeseries.csv")
        )
        series = read_rows(series_path)
        by_time = {round(float(row["time"]), 12): row for row in series}
        selected = [by_time.get(round(float(time), 12)) for time in times]
        if any(row is None for row in selected):
            raise ValueError(f"missing extrapolation time in {series_path}")
        selected_rows = [row for row in selected if row is not None]
        values = np.asarray(
            [
                float(row["cumulative_log_normalizer"]) / float(row["time"])
                for row in selected_rows
            ]
        )
        intercept, slope, r2, rmse = fit_line(times, values)
        prefix = [row for row in series if float(row["time"]) <= args.window_stop]
        final = selected_rows[-1]
        population = int(summary["clone_count"])
        root_weight_ess = float(
            final.get("root_weight_ess", final["root_count_ess"])
        )
        support = (
            int(summary["midpoint_failures"]) == 0
            and min(float(row["weight_ess"]) for row in prefix) >= 0.1 * population
            and int(final["unique_roots"]) >= 8
            and root_weight_ess >= 4.0
        )
        key = (
            summary["mode"], int(summary["n"]), population,
            float(summary["selection_time"]), float(summary["dt"]),
            float(summary["k"]), float(summary["gauge_shift"]),
            float(summary["control_scale"]),
            float(summary.get("resample_threshold", 1.0)),
        )
        curves[key].append(values)
        runs.append(
            {
                "path": str(summary_path),
                "mode": key[0], "n": key[1], "clone_count": key[2],
                "selection_time": key[3], "dt": key[4], "k": key[5],
                "gauge_shift": key[6], "control_scale": key[7],
                "resample_threshold": key[8], "seed": int(summary["seed"]),
                "window_start": args.window_start,
                "window_stop": args.window_stop,
                "window_points": len(times),
                "extrapolated_scgf": intercept,
                "inverse_time_slope": slope,
                "run_fit_r2": r2,
                "run_fit_rmse": rmse,
                "minimum_weight_ess": min(
                    float(row["weight_ess"]) for row in prefix
                ),
                "final_unique_roots": int(final["unique_roots"]),
                "final_root_weight_ess": root_weight_ess,
                "support_gate": int(support),
            }
        )

    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in runs:
        key = (
            row["mode"], row["n"], row["clone_count"],
            row["selection_time"], row["dt"], row["k"],
            row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        )
        grouped[key].append(row)

    aggregate: list[dict[str, object]] = []
    for key, rows in sorted(grouped.items()):
        mean, se, low, high = mean_se_ci(
            [float(row["extrapolated_scgf"]) for row in rows]
        )
        mean_curve = np.mean(np.stack(curves[key]), axis=0)
        _, curve_slope, curve_r2, curve_rmse = fit_line(times, mean_curve)
        curve_range = float(np.ptp(mean_curve))
        fit_shape = curve_r2 >= 0.5 or curve_range <= 0.01
        support = len(rows) >= 4 and all(int(row["support_gate"]) for row in rows)
        aggregate.append(
            {
                "mode": key[0], "n": key[1], "clone_count": key[2],
                "selection_time": key[3], "dt": key[4], "k": key[5],
                "gauge_shift": key[6], "control_scale": key[7],
                "resample_threshold": key[8], "independent_runs": len(rows),
                "window_start": args.window_start,
                "window_stop": args.window_stop,
                "mean_extrapolated_scgf": mean,
                "run_intercept_se": se,
                "intercept_ci_low": low,
                "intercept_ci_high": high,
                "mean_curve_inverse_time_slope": curve_slope,
                "mean_curve_fit_r2": curve_r2,
                "mean_curve_fit_rmse": curve_rmse,
                "mean_curve_range": curve_range,
                "fit_shape_gate": int(fit_shape),
                "support_gate": int(support and fit_shape),
            }
        )

    lookup = {
        (
            row["mode"], row["n"], row["clone_count"],
            row["selection_time"], row["dt"], round(float(row["k"]), 12),
            row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        ): row
        for row in aggregate
    }
    pairs: list[dict[str, object]] = []
    for row in aggregate:
        k = float(row["k"])
        if k >= 0.5:
            continue
        partner = lookup.get(
            (
                row["mode"], row["n"], row["clone_count"],
                row["selection_time"], row["dt"], round(1.0 - k, 12),
                row["gauge_shift"], row["control_scale"],
                row["resample_threshold"],
            )
        )
        if partner is None:
            continue
        delta, se, low, high = welch_interval(
            float(row["mean_extrapolated_scgf"]),
            float(row["run_intercept_se"]), int(row["independent_runs"]),
            float(partner["mean_extrapolated_scgf"]),
            float(partner["run_intercept_se"]),
            int(partner["independent_runs"]),
        )
        support = bool(int(row["support_gate"]) and int(partner["support_gate"]))
        gate = support and abs(delta) <= 0.01 and low <= 0.0 <= high
        pairs.append(
            {
                "mode": row["mode"], "n": row["n"],
                "clone_count": row["clone_count"],
                "selection_time": row["selection_time"], "dt": row["dt"],
                "k": k, "one_minus_k": 1.0 - k,
                "gauge_shift": row["gauge_shift"],
                "control_scale": row["control_scale"],
                "resample_threshold": row["resample_threshold"],
                "extrapolated_residual": delta,
                "residual_se": se,
                "residual_ci_low": low,
                "residual_ci_high": high,
                "support_gate": int(support),
                "extrapolated_gc_gate": int(gate),
            }
        )

    population_members: list[dict[str, object]] = []
    pop_groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in aggregate:
        key = (
            row["mode"], row["n"], row["selection_time"], row["dt"],
            row["k"], row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        )
        pop_groups[key].append(row)
    for key, rows in sorted(pop_groups.items()):
        ordered = sorted(rows, key=lambda item: int(item["clone_count"]))
        for low, high in zip(ordered, ordered[1:]):
            if int(high["clone_count"]) != 2 * int(low["clone_count"]):
                continue
            delta = (
                float(high["mean_extrapolated_scgf"])
                - float(low["mean_extrapolated_scgf"])
            )
            support = bool(int(low["support_gate"]) and int(high["support_gate"]))
            gate = support and convergence_gate(
                delta, float(low["run_intercept_se"]),
                float(high["run_intercept_se"]),
            )
            population_members.append(
                {
                    "mode": key[0], "n": key[1],
                    "selection_time": key[2], "dt": key[3], "k": key[4],
                    "gauge_shift": key[5], "control_scale": key[6],
                    "resample_threshold": key[7],
                    "lower_clone_count": low["clone_count"],
                    "upper_clone_count": high["clone_count"],
                    "intercept_delta_upper_minus_lower": delta,
                    "combined_se": math.hypot(
                        float(low["run_intercept_se"]),
                        float(high["run_intercept_se"]),
                    ),
                    "support_gate": int(support),
                    "population_gate": int(gate),
                }
            )

    pair_populations: list[dict[str, object]] = []
    pair_groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pairs:
        key = (
            row["mode"], row["n"], row["selection_time"], row["dt"],
            row["k"], row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        )
        pair_groups[key].append(row)
    for key, rows in sorted(pair_groups.items()):
        ordered = sorted(rows, key=lambda item: int(item["clone_count"]))
        for low, high in zip(ordered, ordered[1:]):
            if int(high["clone_count"]) != 2 * int(low["clone_count"]):
                continue
            delta = (
                float(high["extrapolated_residual"])
                - float(low["extrapolated_residual"])
            )
            support = bool(int(low["support_gate"]) and int(high["support_gate"]))
            gate = support and convergence_gate(
                delta, float(low["residual_se"]), float(high["residual_se"])
            )
            pair_populations.append(
                {
                    "mode": key[0], "n": key[1],
                    "selection_time": key[2], "dt": key[3], "k": key[4],
                    "gauge_shift": key[5], "control_scale": key[6],
                    "resample_threshold": key[7],
                    "lower_clone_count": low["clone_count"],
                    "upper_clone_count": high["clone_count"],
                    "residual_delta_upper_minus_lower": delta,
                    "combined_se": math.hypot(
                        float(low["residual_se"]), float(high["residual_se"])
                    ),
                    "support_gate": int(support),
                    "pair_population_gate": int(gate),
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "run_extrapolations.csv", runs)
    write_csv(args.output_dir / "aggregate_extrapolations.csv", aggregate)
    write_csv(args.output_dir / "paired_extrapolations.csv", pairs)
    write_csv(args.output_dir / "population_members.csv", population_members)
    write_csv(args.output_dir / "population_pairs.csv", pair_populations)
    audit = {
        "window": [float(time) for time in times],
        "runs": len(runs),
        "groups": len(aggregate),
        "pairs": len(pairs),
        "passing_pairs": sum(int(row["extrapolated_gc_gate"]) for row in pairs),
        "population_member_checks": len(population_members),
        "passing_population_member_checks": sum(
            int(row["population_gate"]) for row in population_members
        ),
        "population_pair_checks": len(pair_populations),
        "passing_population_pair_checks": sum(
            int(row["pair_population_gate"]) for row in pair_populations
        ),
    }
    (args.output_dir / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
