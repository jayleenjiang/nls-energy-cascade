#!/usr/bin/env python3
"""Aggregate independent entropy-cloning runs and compare with direct SCGF."""

from __future__ import annotations

import argparse
import csv
import glob
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_one(path: Path) -> dict[str, str] | None:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "clone_count" not in reader.fieldnames:
            return None
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"expected one summary row in {path}")
    row = rows[0]
    row["_path"] = str(path)
    row.setdefault("mode", "naive")
    return row


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


def stable_log_mean(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    return maximum + math.log(float(np.exp(values - maximum).mean()))


def load_direct(path: Path | None) -> dict[tuple[int, float, float], dict[str, str]]:
    if path is None:
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        (int(row["n"]), float(row["tau"]), round(float(row["k"]), 12)): row
        for row in rows
    }


def late_slope(summary_path: Path) -> float:
    timeseries_path = Path(str(summary_path).replace("_summary.csv", "_timeseries.csv"))
    with timeseries_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    time = np.asarray([float(row["time"]) for row in rows])
    log_z = np.asarray([float(row["cumulative_log_normalizer"]) for row in rows])
    keep = time >= 0.5 * time[-1]
    if np.count_nonzero(keep) < 3:
        return math.nan
    return float(np.polyfit(time[keep], log_z[keep], 1)[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patterns", nargs="+", help="glob patterns for summary CSVs")
    parser.add_argument("--direct-scgf", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    paths: list[Path] = []
    for pattern in args.patterns:
        paths.extend(Path(value) for value in glob.glob(pattern, recursive=True))
    paths = sorted(set(paths))
    if not paths:
        raise SystemExit("no cloning summaries matched")

    direct = load_direct(args.direct_scgf)
    raw_rows: list[dict[str, object]] = []
    groups: defaultdict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    skipped: list[Path] = []
    for path in paths:
        row = read_one(path)
        if row is None:
            skipped.append(path)
            continue
        parsed: dict[str, object] = {
            "path": str(path),
            "mode": row["mode"],
            "n": int(row["n"]),
            "clone_count": int(row["clone_count"]),
            "burnin": float(row["burnin"]),
            "observation_time": float(row["observation_time"]),
            "selection_time": float(row["selection_time"]),
            "dt": float(row["dt"]),
            "k": float(row["k"]),
            "seed": int(row["seed"]),
            "scgf": float(row["scgf"]),
            "late_half_logz_slope": late_slope(path),
            "mean_population_entropy_rate": float(row["mean_population_entropy_rate"]),
            "mean_population_action_current": float(row["mean_population_action_current"]),
            "minimum_weight_ess": float(row["minimum_weight_ess"]),
            "minimum_root_count_ess": float(row["minimum_root_count_ess"]),
            "minimum_unique_roots": int(row["minimum_unique_roots"]),
            "midpoint_failures": int(row["midpoint_failures"]),
        }
        raw_rows.append(parsed)
        key = (
            parsed["mode"], parsed["n"], parsed["clone_count"],
            parsed["observation_time"], parsed["selection_time"],
            parsed["dt"], parsed["k"],
        )
        groups[key].append(parsed)

    aggregate: list[dict[str, object]] = []
    for key, rows in sorted(groups.items()):
        mode, n, clone_count, observation, selection, dt, k = key
        values = np.asarray([float(row["scgf"]) for row in rows])
        late = np.asarray([float(row["late_half_logz_slope"]) for row in rows])
        log_normalizers = float(observation) * values
        direct_row = direct.get((int(n), float(observation), round(float(k), 12)))
        mean = float(values.mean())
        se = (
            float(values.std(ddof=1) / math.sqrt(values.size))
            if values.size > 1 else math.nan
        )
        finite_late = late[np.isfinite(late)]
        pooled = stable_log_mean(log_normalizers) / float(observation)
        direct_psi = float(direct_row["psi"]) if direct_row else math.nan
        aggregate.append(
            {
                "mode": mode,
                "n": n,
                "clone_count": clone_count,
                "observation_time": observation,
                "selection_time": selection,
                "dt": dt,
                "k": k,
                "independent_runs": values.size,
                "mean_run_scgf": mean,
                "run_scgf_se": se,
                "pooled_normalizer_scgf": pooled,
                "mean_late_half_logz_slope": (
                    float(finite_late.mean()) if finite_late.size else math.nan
                ),
                "late_half_slope_se": (
                    float(finite_late.std(ddof=1) /
                          math.sqrt(finite_late.size))
                    if finite_late.size > 1 else math.nan
                ),
                "direct_scgf": direct_psi,
                "direct_ci_low": (
                    float(direct_row["psi_ci_low"]) if direct_row else math.nan
                ),
                "direct_ci_high": (
                    float(direct_row["psi_ci_high"]) if direct_row else math.nan
                ),
                "direct_reliability_gate": (
                    int(direct_row["direct_reliability_gate"])
                    if direct_row else 0
                ),
                "mean_minus_direct": mean - direct_psi if direct_row else math.nan,
                "pooled_minus_direct": (
                    pooled - direct_psi if direct_row else math.nan
                ),
                "minimum_weight_ess_across_runs": min(
                    float(row["minimum_weight_ess"]) for row in rows
                ),
                "minimum_root_ess_across_runs": min(
                    float(row["minimum_root_count_ess"]) for row in rows
                ),
                "minimum_unique_roots_across_runs": min(
                    int(row["minimum_unique_roots"]) for row in rows
                ),
                "total_midpoint_failures": sum(
                    int(row["midpoint_failures"]) for row in rows
                ),
            }
        )

    output = args.output_dir
    write_csv(output / "cloning_runs.csv", raw_rows)
    write_csv(output / "cloning_aggregate.csv", aggregate)

    comparable = [
        row for row in aggregate
        if int(row["n"]) == 10
        and float(row["observation_time"]) == 20.0
        and int(row["direct_reliability_gate"]) == 1
    ]
    if comparable:
        figure, axis = plt.subplots(figsize=(7.2, 4.8))
        direct_points: dict[float, float] = {}
        for row in comparable:
            k = float(row["k"])
            direct_points[k] = float(row["direct_scgf"])
            marker = "o" if row["mode"] == "guided" else "s"
            label = f"{row['mode']}, N={row['clone_count']}"
            axis.errorbar(
                [k], [float(row["mean_run_scgf"])],
                yerr=[float(row["run_scgf_se"])], fmt=marker,
                capsize=3, label=label,
            )
        k_sorted = sorted(direct_points)
        axis.plot(k_sorted, [direct_points[k] for k in k_sorted],
                  "k--", lw=1.5, label="direct million-sample")
        handles, labels = axis.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        axis.legend(unique.values(), unique.keys(), fontsize=8)
        axis.set_xlabel(r"tilt $k$")
        axis.set_ylabel(r"finite-time $\psi_{20}(k)$")
        axis.grid(alpha=0.25)
        figure.tight_layout()
        figure.savefig(output / "direct_cloning_crosscheck.png", dpi=220)
        figure.savefig(output / "direct_cloning_crosscheck.pdf")
        plt.close(figure)

    population = sorted(
        (
            row for row in aggregate
            if row["mode"] == "guided" and int(row["n"]) == 10
            and float(row["observation_time"]) == 20.0
            and float(row["selection_time"]) == 0.5
            and float(row["dt"]) == 0.0005
            and abs(float(row["k"]) - 0.3) < 1.0e-12
            and int(row["independent_runs"]) >= 4
        ),
        key=lambda row: int(row["clone_count"]),
    )
    if population:
        figure, axis = plt.subplots(figsize=(6.6, 4.5))
        clone_counts = [int(row["clone_count"]) for row in population]
        means = [float(row["mean_run_scgf"]) for row in population]
        errors = [float(row["run_scgf_se"]) for row in population]
        axis.errorbar(clone_counts, means, yerr=errors, fmt="o-", capsize=4,
                      label="guided cloning (mean $\\pm$ SE)")
        direct_value = float(population[-1]["direct_scgf"])
        direct_low = float(population[-1]["direct_ci_low"])
        direct_high = float(population[-1]["direct_ci_high"])
        axis.axhline(direct_value, color="black", ls="--",
                     label="direct million-sample estimate")
        axis.axhspan(direct_low, direct_high, color="black", alpha=0.12)
        axis.set_xscale("log", base=2)
        axis.set_xticks(clone_counts, labels=[str(value) for value in clone_counts])
        axis.set_xlabel(r"clone population $N_c$")
        axis.set_ylabel(r"finite-time $\psi_{20}(0.3)$")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
        figure.tight_layout()
        figure.savefig(output / "population_convergence_k0p3.png", dpi=220)
        figure.savefig(output / "population_convergence_k0p3.pdf")
        plt.close(figure)

    print(f"Aggregated {len(raw_rows)} runs into {len(aggregate)} groups")
    if skipped:
        print(f"Skipped {len(skipped)} non-cloning summary file(s)")


if __name__ == "__main__":
    main()
