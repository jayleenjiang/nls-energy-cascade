#!/usr/bin/env python3
"""Aggregate controlled-cloning runs for long-chain GC acceptance gates."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as student_t


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "experiment_roots", nargs="+", type=Path,
        help="one or more roots containing controlled-cloning summaries",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--horizons", nargs="+", type=float, default=[20, 40, 80])
    return parser.parse_args()


def read_single(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected one summary row in {path}")
    return rows[0]


def read_timeseries(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty timeseries {path}")
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def exact_row(rows: list[dict[str, str]], time: float) -> dict[str, str]:
    matches = [row for row in rows if math.isclose(float(row["time"]), time)]
    if len(matches) != 1:
        raise ValueError(f"expected one timeseries row at t={time}, got {len(matches)}")
    return matches[0]


def mean_se_ci(values: list[float]) -> tuple[float, float, float, float]:
    array = np.asarray(values, dtype=float)
    mean = float(array.mean())
    if array.size < 2:
        return mean, math.nan, math.nan, math.nan
    se = float(array.std(ddof=1) / math.sqrt(array.size))
    critical = float(student_t.ppf(0.975, array.size - 1))
    return mean, se, mean - critical * se, mean + critical * se


def welch_interval(
    mean_a: float,
    se_a: float,
    n_a: int,
    mean_b: float,
    se_b: float,
    n_b: int,
) -> tuple[float, float, float, float]:
    residual = mean_a - mean_b
    variance_a = se_a * se_a
    variance_b = se_b * se_b
    se = math.sqrt(variance_a + variance_b)
    if not math.isfinite(se) or se == 0.0 or n_a < 2 or n_b < 2:
        return residual, se, math.nan, math.nan
    denominator = (
        variance_a * variance_a / (n_a - 1)
        + variance_b * variance_b / (n_b - 1)
    )
    degrees = (variance_a + variance_b) ** 2 / denominator
    critical = float(student_t.ppf(0.975, degrees))
    return residual, se, residual - critical * se, residual + critical * se


def convergence_gate(delta: float, se_a: float, se_b: float) -> bool:
    """Apply the predeclared absolute and two-combined-SE convergence limits."""
    combined = math.hypot(se_a, se_b)
    return (
        math.isfinite(delta)
        and math.isfinite(combined)
        and abs(delta) <= 0.01
        and abs(delta) <= 2.0 * combined
    )


def main() -> None:
    args = parse_args()
    run_rows: list[dict[str, object]] = []
    summary_paths = sorted(
        {
            path.resolve()
            for root in args.experiment_roots
            for path in root.rglob("*_summary.csv")
        }
    )
    for summary_path in summary_paths:
        summary = read_single(summary_path)
        if "clone_count" not in summary or not summary.get("mode", "").startswith(
            "controlled_exact"
        ):
            continue
        timeseries_path = summary_path.with_name(
            summary_path.name.replace("_summary.csv", "_timeseries.csv")
        )
        timeseries = read_timeseries(timeseries_path)
        observation = float(summary["observation_time"])
        for horizon in args.horizons:
            if horizon > observation + 1.0e-12:
                continue
            current = exact_row(timeseries, horizon)
            half = exact_row(timeseries, horizon / 2.0)
            cumulative = float(current["cumulative_log_normalizer"])
            cumulative_half = float(half["cumulative_log_normalizer"])
            prefix = [row for row in timeseries if float(row["time"]) <= horizon]
            root_weight = float(
                current.get("root_weight_ess", current["root_count_ess"])
            )
            run_rows.append(
                {
                    "path": str(summary_path),
                    "mode": summary["mode"],
                    "n": int(summary["n"]),
                    "clone_count": int(summary["clone_count"]),
                    "horizon": horizon,
                    "selection_time": float(summary["selection_time"]),
                    "dt": float(summary["dt"]),
                    "k": float(summary["k"]),
                    "gauge_shift": float(summary["gauge_shift"]),
                    "control_scale": float(summary["control_scale"]),
                    "resample_threshold": float(
                        summary.get("resample_threshold", 1.0)
                    ),
                    "seed": int(summary["seed"]),
                    "scgf": cumulative / horizon,
                    "late_half_scgf":
                        (cumulative - cumulative_half) / (horizon / 2.0),
                    "minimum_weight_ess_to_horizon": min(
                        float(row["weight_ess"]) for row in prefix
                    ),
                    "final_unique_roots": int(current["unique_roots"]),
                    "final_root_weight_ess": root_weight,
                    "midpoint_failures": int(summary["midpoint_failures"]),
                }
            )
    if not run_rows:
        raise SystemExit("no controlled-cloning runs found")

    groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in run_rows:
        key = (
            row["mode"], row["n"], row["clone_count"], row["horizon"],
            row["selection_time"], row["dt"], row["k"],
            row["gauge_shift"], row["control_scale"], row["resample_threshold"],
        )
        groups[key].append(row)

    aggregate: list[dict[str, object]] = []
    for key, rows in sorted(groups.items()):
        (
            mode, n, clone_count, horizon, selection, dt, k,
            gauge, control, resample,
        ) = key
        mean, se, ci_low, ci_high = mean_se_ci(
            [float(row["scgf"]) for row in rows]
        )
        late_mean, late_se, late_low, late_high = mean_se_ci(
            [float(row["late_half_scgf"]) for row in rows]
        )
        numerical = all(int(row["midpoint_failures"]) == 0 for row in rows)
        weight_support = all(
            float(row["minimum_weight_ess_to_horizon"]) >= 0.1 * int(clone_count)
            for row in rows
        )
        calibration_support = horizon > 20.0 or all(
            int(row["final_unique_roots"]) >= 32
            and float(row["final_root_weight_ess"]) >= 16.0
            for row in rows
        )
        aggregate.append(
            {
                "mode": mode,
                "n": n,
                "clone_count": clone_count,
                "horizon": horizon,
                "selection_time": selection,
                "dt": dt,
                "k": k,
                "gauge_shift": gauge,
                "control_scale": control,
                "resample_threshold": resample,
                "independent_runs": len(rows),
                "mean_scgf": mean,
                "run_scgf_se": se,
                "scgf_ci_low": ci_low,
                "scgf_ci_high": ci_high,
                "mean_late_half_scgf": late_mean,
                "late_half_scgf_se": late_se,
                "late_half_ci_low": late_low,
                "late_half_ci_high": late_high,
                "minimum_weight_ess": min(
                    float(row["minimum_weight_ess_to_horizon"]) for row in rows
                ),
                "minimum_final_unique_roots": min(
                    int(row["final_unique_roots"]) for row in rows
                ),
                "minimum_final_root_weight_ess": min(
                    float(row["final_root_weight_ess"]) for row in rows
                ),
                "support_gate": int(
                    len(rows) >= 4 and numerical and weight_support
                    and calibration_support
                ),
            }
        )

    lookup = {
        (
            row["mode"], row["n"], row["clone_count"], row["horizon"],
            row["selection_time"], row["dt"], round(float(row["k"]), 12),
            row["gauge_shift"], row["control_scale"], row["resample_threshold"],
        ): row
        for row in aggregate
    }
    pairs: list[dict[str, object]] = []
    for row in aggregate:
        k = float(row["k"])
        if k >= 0.5:
            continue
        partner_key = (
            row["mode"], row["n"], row["clone_count"], row["horizon"],
            row["selection_time"], row["dt"], round(1.0 - k, 12),
            row["gauge_shift"], row["control_scale"], row["resample_threshold"],
        )
        partner = lookup.get(partner_key)
        if partner is None:
            continue
        residual, residual_se, ci_low, ci_high = welch_interval(
            float(row["mean_scgf"]), float(row["run_scgf_se"]),
            int(row["independent_runs"]), float(partner["mean_scgf"]),
            float(partner["run_scgf_se"]), int(partner["independent_runs"]),
        )
        late_residual, late_se, late_low, late_high = welch_interval(
            float(row["mean_late_half_scgf"]),
            float(row["late_half_scgf_se"]), int(row["independent_runs"]),
            float(partner["mean_late_half_scgf"]),
            float(partner["late_half_scgf_se"]),
            int(partner["independent_runs"]),
        )
        support = bool(int(row["support_gate"]) and int(partner["support_gate"]))
        gate = (
            support and abs(residual) <= 0.01 and ci_low <= 0.0 <= ci_high
            and abs(late_residual) <= 0.01
            and late_low <= 0.0 <= late_high
        )
        pairs.append(
            {
                "mode": row["mode"],
                "n": row["n"],
                "clone_count": row["clone_count"],
                "horizon": row["horizon"],
                "selection_time": row["selection_time"],
                "dt": row["dt"],
                "k": k,
                "one_minus_k": 1.0 - k,
                "gauge_shift": row["gauge_shift"],
                "control_scale": row["control_scale"],
                "resample_threshold": row["resample_threshold"],
                "residual": residual,
                "residual_se": residual_se,
                "residual_ci_low": ci_low,
                "residual_ci_high": ci_high,
                "late_half_residual": late_residual,
                "late_half_residual_se": late_se,
                "late_half_ci_low": late_low,
                "late_half_ci_high": late_high,
                "paired_support_gate": int(support),
                "paired_gc_gate": int(gate),
            }
        )

    population_rows: list[dict[str, object]] = []
    population_groups: dict[tuple[object, ...], list[dict[str, object]]] = (
        defaultdict(list)
    )
    for row in aggregate:
        key = (
            row["mode"], row["n"], row["horizon"], row["selection_time"],
            row["dt"], row["k"], row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        )
        population_groups[key].append(row)
    for key, rows in sorted(population_groups.items()):
        ordered = sorted(rows, key=lambda item: int(item["clone_count"]))
        for low, high in zip(ordered, ordered[1:]):
            low_n = int(low["clone_count"])
            high_n = int(high["clone_count"])
            if high_n != 2 * low_n:
                continue
            delta = float(high["mean_scgf"]) - float(low["mean_scgf"])
            late_delta = (
                float(high["mean_late_half_scgf"])
                - float(low["mean_late_half_scgf"])
            )
            support = bool(int(low["support_gate"]) and int(high["support_gate"]))
            gate = (
                support
                and convergence_gate(
                    delta, float(low["run_scgf_se"]),
                    float(high["run_scgf_se"]),
                )
                and convergence_gate(
                    late_delta, float(low["late_half_scgf_se"]),
                    float(high["late_half_scgf_se"]),
                )
            )
            population_rows.append(
                {
                    "mode": key[0],
                    "n": key[1],
                    "horizon": key[2],
                    "selection_time": key[3],
                    "dt": key[4],
                    "k": key[5],
                    "gauge_shift": key[6],
                    "control_scale": key[7],
                    "resample_threshold": key[8],
                    "lower_clone_count": low_n,
                    "upper_clone_count": high_n,
                    "scgf_delta_upper_minus_lower": delta,
                    "combined_scgf_se": math.hypot(
                        float(low["run_scgf_se"]),
                        float(high["run_scgf_se"]),
                    ),
                    "late_half_delta_upper_minus_lower": late_delta,
                    "combined_late_half_se": math.hypot(
                        float(low["late_half_scgf_se"]),
                        float(high["late_half_scgf_se"]),
                    ),
                    "both_support_gate": int(support),
                    "population_gate": int(gate),
                }
            )

    pair_population_rows: list[dict[str, object]] = []
    pair_population_groups: dict[tuple[object, ...], list[dict[str, object]]] = (
        defaultdict(list)
    )
    for row in pairs:
        key = (
            row["mode"], row["n"], row["horizon"], row["selection_time"],
            row["dt"], row["k"], row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        )
        pair_population_groups[key].append(row)
    for key, rows in sorted(pair_population_groups.items()):
        ordered = sorted(rows, key=lambda item: int(item["clone_count"]))
        for low, high in zip(ordered, ordered[1:]):
            low_n = int(low["clone_count"])
            high_n = int(high["clone_count"])
            if high_n != 2 * low_n:
                continue
            delta = float(high["residual"]) - float(low["residual"])
            late_delta = (
                float(high["late_half_residual"])
                - float(low["late_half_residual"])
            )
            support = bool(
                int(low["paired_support_gate"])
                and int(high["paired_support_gate"])
            )
            gate = (
                support
                and convergence_gate(
                    delta, float(low["residual_se"]),
                    float(high["residual_se"]),
                )
                and convergence_gate(
                    late_delta, float(low["late_half_residual_se"]),
                    float(high["late_half_residual_se"]),
                )
            )
            pair_population_rows.append(
                {
                    "mode": key[0],
                    "n": key[1],
                    "horizon": key[2],
                    "selection_time": key[3],
                    "dt": key[4],
                    "k": key[5],
                    "gauge_shift": key[6],
                    "control_scale": key[7],
                    "resample_threshold": key[8],
                    "lower_clone_count": low_n,
                    "upper_clone_count": high_n,
                    "residual_delta_upper_minus_lower": delta,
                    "combined_residual_se": math.hypot(
                        float(low["residual_se"]),
                        float(high["residual_se"]),
                    ),
                    "late_half_delta_upper_minus_lower": late_delta,
                    "combined_late_half_se": math.hypot(
                        float(low["late_half_residual_se"]),
                        float(high["late_half_residual_se"]),
                    ),
                    "both_support_gate": int(support),
                    "pair_population_gate": int(gate),
                }
            )

    time_rows: list[dict[str, object]] = []
    time_groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pairs:
        key = (
            row["mode"], row["n"], row["clone_count"],
            row["selection_time"], row["dt"], row["k"],
            row["gauge_shift"], row["control_scale"],
            row["resample_threshold"],
        )
        time_groups[key].append(row)
    for key, rows in sorted(time_groups.items()):
        ordered = sorted(rows, key=lambda item: float(item["horizon"]))
        if len(ordered) < 3:
            continue
        penultimate, final = ordered[-2:]
        delta = float(final["residual"]) - float(penultimate["residual"])
        late_delta = (
            float(final["late_half_residual"])
            - float(penultimate["late_half_residual"])
        )
        support = bool(
            int(penultimate["paired_support_gate"])
            and int(final["paired_support_gate"])
        )
        gate = (
            support
            and convergence_gate(
                delta, float(penultimate["residual_se"]),
                float(final["residual_se"]),
            )
            and convergence_gate(
                late_delta, float(penultimate["late_half_residual_se"]),
                float(final["late_half_residual_se"]),
            )
        )
        time_rows.append(
            {
                "mode": key[0],
                "n": key[1],
                "clone_count": key[2],
                "selection_time": key[3],
                "dt": key[4],
                "k": key[5],
                "gauge_shift": key[6],
                "control_scale": key[7],
                "resample_threshold": key[8],
                "horizons_available": len(ordered),
                "penultimate_horizon": penultimate["horizon"],
                "final_horizon": final["horizon"],
                "residual_delta_final_minus_penultimate": delta,
                "combined_residual_se": math.hypot(
                    float(penultimate["residual_se"]),
                    float(final["residual_se"]),
                ),
                "late_half_delta_final_minus_penultimate": late_delta,
                "combined_late_half_se": math.hypot(
                    float(penultimate["late_half_residual_se"]),
                    float(final["late_half_residual_se"]),
                ),
                "both_support_gate": int(support),
                "time_convergence_gate": int(gate),
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "controlled_runs.csv", run_rows)
    write_csv(args.output_dir / "controlled_aggregate.csv", aggregate)
    if pairs:
        write_csv(args.output_dir / "paired_gc_residuals.csv", pairs)
    if population_rows:
        write_csv(args.output_dir / "population_convergence_members.csv", population_rows)
    if pair_population_rows:
        write_csv(
            args.output_dir / "population_convergence_pairs.csv",
            pair_population_rows,
        )
    if time_rows:
        write_csv(args.output_dir / "time_convergence_pairs.csv", time_rows)

    for n in sorted({int(row["n"]) for row in pairs}):
        selected = [row for row in pairs if int(row["n"]) == n]
        figure, axis = plt.subplots(figsize=(6.8, 4.5))
        for key in sorted({(int(row["clone_count"]), float(row["k"])) for row in selected}):
            population, k = key
            rows = sorted(
                [row for row in selected if int(row["clone_count"]) == population
                 and math.isclose(float(row["k"]), k)],
                key=lambda row: float(row["horizon"]),
            )
            axis.errorbar(
                [float(row["horizon"]) for row in rows],
                [float(row["residual"]) for row in rows],
                yerr=[float(row["residual_se"]) for row in rows],
                marker="o", capsize=3, label=rf"$N_c={population},k={k:g}$",
            )
        axis.axhline(0.0, color="black", ls="--", lw=1)
        axis.set_xlabel("observation horizon")
        axis.set_ylabel(r"$\psi_t(k)-\psi_t(1-k)$")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
        figure.tight_layout()
        figure.savefig(args.output_dir / f"paired_gc_convergence_n{n}.png", dpi=220)
        plt.close(figure)

    audit = {
        "runs": len(run_rows),
        "groups": len(aggregate),
        "pairs": len(pairs),
        "supported_pairs": sum(int(row["paired_support_gate"]) for row in pairs),
        "passing_pairs": sum(int(row["paired_gc_gate"]) for row in pairs),
        "population_member_comparisons": len(population_rows),
        "passing_population_member_comparisons": sum(
            int(row["population_gate"]) for row in population_rows
        ),
        "population_pair_comparisons": len(pair_population_rows),
        "passing_population_pair_comparisons": sum(
            int(row["pair_population_gate"]) for row in pair_population_rows
        ),
        "time_convergence_comparisons": len(time_rows),
        "passing_time_convergence_comparisons": sum(
            int(row["time_convergence_gate"]) for row in time_rows
        ),
    }
    (args.output_dir / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
