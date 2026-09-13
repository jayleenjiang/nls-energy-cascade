#!/usr/bin/env python3
"""Prospective population-bias analysis for controlled NLS cloning."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
MATRIX = ROOT / "RUN_MATRIX.csv"
ANALYSIS = ROOT / "analysis"
PRIMARY_POPULATIONS = (512, 1024, 2048, 4096)
BOOTSTRAPS = 20_000


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def single_csv_row(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def close(a: float, b: float, tolerance: float = 5.0e-13) -> bool:
    return abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))


def collect_runs(require_primary: bool) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    seen_seeds: set[int] = set()
    for row in read_csv(MATRIX):
        prefix = ROOT / row["prefix"]
        summary_path = Path(f"{prefix}_summary.csv")
        timeseries_path = Path(f"{prefix}_timeseries.csv")
        command_path = Path(f"{prefix}.command.txt")
        exists = summary_path.exists() or timeseries_path.exists()
        complete = summary_path.stat().st_size > 0 if summary_path.exists() else False
        complete = complete and timeseries_path.exists() and timeseries_path.stat().st_size > 0
        if not complete:
            if exists:
                raise RuntimeError(f"incomplete run artifact at {prefix}")
            if require_primary and row["stage"] == "primary":
                raise RuntimeError(f"missing primary run at {prefix}")
            continue

        summary = single_csv_row(summary_path)
        expected = {
            "mode": "controlled_exact",
            "n": int(row["n"]),
            "clone_count": int(row["clone_count"]),
            "seed": int(row["seed"]),
        }
        for key, value in expected.items():
            actual: object = summary[key]
            if isinstance(value, int):
                actual = int(actual)
            if actual != value:
                raise RuntimeError(f"{prefix}: {key}={actual}, expected {value}")
        for key in ("burnin", "observation_time", "selection_time", "dt", "k",
                    "gauge_shift", "control_scale"):
            matrix_key = "horizon" if key == "observation_time" else key
            if not close(float(summary[key]), float(row[matrix_key])):
                raise RuntimeError(f"{prefix}: parameter mismatch for {key}")
        if not close(float(summary["observable_left_heat_coefficient"]), 0.0):
            raise RuntimeError(f"{prefix}: left heat coefficient is not zero")
        if not close(float(summary["observable_right_heat_coefficient"]), -0.4):
            raise RuntimeError(f"{prefix}: right heat coefficient is not -0.4")

        final_ts: dict[str, str] | None = None
        checkpoint: dict[str, str] | None = None
        with timeseries_path.open(newline="") as handle:
            for value in csv.DictReader(handle):
                final_ts = value
                if close(float(value["time"]), 20.0):
                    checkpoint = value
        if final_ts is None:
            raise RuntimeError(f"empty timeseries {timeseries_path}")
        if not close(float(final_ts["time"]), float(row["horizon"])):
            raise RuntimeError(f"wrong final time in {timeseries_path}")
        if not close(float(final_ts["scgf"]), float(summary["scgf"]), 2.0e-12):
            raise RuntimeError(f"summary/timeseries SCGF mismatch at {prefix}")

        population = int(row["clone_count"])
        midpoint_failures = int(summary["midpoint_failures"])
        minimum_weight_ess = float(summary["minimum_weight_ess"])
        finite = all(math.isfinite(float(summary[key])) for key in (
            "scgf", "minimum_weight_ess", "mean_midpoint_iterations",
            "mean_hamiltonian_energy_error_rate"))
        support_pass = (
            finite and midpoint_failures == 0
            and minimum_weight_ess >= 0.1 * population
        )
        seed = int(row["seed"])
        if seed in seen_seeds:
            raise RuntimeError(f"duplicate seed {seed}")
        seen_seeds.add(seed)
        records.append({
            "study": row["study"],
            "stage": row["stage"],
            "n": int(row["n"]),
            "clone_count": population,
            "k": float(row["k"]),
            "run": int(row["run"]),
            "seed": seed,
            "horizon": float(row["horizon"]),
            "selection_time": float(row["selection_time"]),
            "dt": float(row["dt"]),
            "psi_final": float(summary["scgf"]),
            "psi_t20": float(checkpoint["scgf"]) if checkpoint else math.nan,
            "minimum_weight_ess": minimum_weight_ess,
            "minimum_weight_ess_fraction": minimum_weight_ess / population,
            "minimum_unique_roots": int(summary["minimum_unique_roots"]),
            "minimum_root_weight_ess": float(summary["minimum_root_weight_ess"]),
            "midpoint_failures": midpoint_failures,
            "support_pass": int(support_pass),
            "summary_path": str(summary_path.relative_to(REPO)),
            "timeseries_path": str(timeseries_path.relative_to(REPO)),
            "command_path": str(command_path.relative_to(REPO)),
        })
    return records


def mean_se(values: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(values))
    se = float(np.std(values, ddof=1) / math.sqrt(values.size)) if values.size > 1 else math.nan
    return mean, se


def fit_hc3(populations: np.ndarray, values: np.ndarray) -> dict[str, float]:
    if values.size < 8 or np.unique(populations).size < 4:
        raise RuntimeError("intercept fit needs at least four populations and eight runs")
    x = 1.0 / populations.astype(float)
    design = np.column_stack((np.ones(values.size), x))
    inv = np.linalg.inv(design.T @ design)
    beta = inv @ design.T @ values
    residual = values - design @ beta
    leverage = np.einsum("ij,jk,ik->i", design, inv, design)
    scaled = residual / (1.0 - leverage)
    meat = design.T @ (design * np.square(scaled)[:, None])
    covariance = inv @ meat @ inv
    se = np.sqrt(np.diag(covariance))
    df = values.size - 2
    critical = float(student_t.ppf(0.975, df))
    fitted = design @ beta
    total = float(np.sum(np.square(values - np.mean(values))))
    r2 = 1.0 - float(np.sum(np.square(residual))) / total if total > 0 else math.nan
    return {
        "intercept": float(beta[0]),
        "slope": float(beta[1]),
        "intercept_se": float(se[0]),
        "slope_se": float(se[1]),
        "intercept_ci_low": float(beta[0] - critical * se[0]),
        "intercept_ci_high": float(beta[0] + critical * se[0]),
        "df": float(df),
        "r2": r2,
        "hc3_var_intercept": float(covariance[0, 0]),
        "max_abs_residual": float(np.max(np.abs(values - fitted))),
    }


def analysis_series(records: list[dict[str, object]]) -> list[tuple[str, list[dict[str, object]], str]]:
    series: list[tuple[str, list[dict[str, object]], str]] = []
    n2 = [r for r in records if r["study"] == "n2_known_answer"]
    series.append(("n2_t0p1", n2, "psi_final"))
    for n in (10, 20, 30, 40):
        subset = [r for r in records if r["study"] == "long_chain" and r["n"] == n]
        series.append((f"n{n}_final", subset, "psi_final"))
    n10_checkpoint = [
        r for r in records
        if r["study"] == "long_chain" and r["n"] == 10 and close(float(r["k"]), 0.3)
    ]
    series.append(("n10_t20_k0p3", n10_checkpoint, "psi_t20"))
    return series


def build_population_tables(records: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[tuple[str, float], dict[str, float]]]:
    levels: list[dict[str, object]] = []
    fits: list[dict[str, object]] = []
    fit_map: dict[tuple[str, float], dict[str, float]] = {}
    for context, context_records, value_key in analysis_series(records):
        by_k: dict[float, list[dict[str, object]]] = defaultdict(list)
        for record in context_records:
            by_k[float(record["k"])].append(record)
        for k, members in sorted(by_k.items()):
            populations = np.asarray([int(r["clone_count"]) for r in members])
            values = np.asarray([float(r[value_key]) for r in members])
            for population in sorted(np.unique(populations)):
                selected = values[populations == population]
                mean, se = mean_se(selected)
                support = all(int(r["support_pass"]) for r in members if int(r["clone_count"]) == population)
                levels.append({
                    "context": context, "k": k, "clone_count": int(population),
                    "independent_runs": int(selected.size), "mean_psi": mean,
                    "run_se": se, "support_pass": int(support),
                })
            fit = fit_hc3(populations, values)
            fit_map[(context, k)] = fit
            fits.append({"context": context, "k": k, "runs": values.size,
                         "population_levels": np.unique(populations).size, **fit})
    return levels, fits, fit_map


def plateau_rows(levels: list[dict[str, object]], fit_map: dict[tuple[str, float], dict[str, float]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float], list[dict[str, object]]] = defaultdict(list)
    for row in levels:
        grouped[(str(row["context"]), float(row["k"]))].append(row)
    output: list[dict[str, object]] = []
    for key, rows in sorted(grouped.items()):
        context, k = key
        rows = sorted(rows, key=lambda value: int(value["clone_count"]))
        if len(rows) < 4:
            raise RuntimeError(f"too few population levels for {key}")
        low, high = rows[-2], rows[-1]
        combined_se = math.sqrt(float(low["run_se"]) ** 2 + float(high["run_se"]) ** 2)
        adjacent_delta = float(high["mean_psi"]) - float(low["mean_psi"])
        fit = fit_map[key]
        high_population = int(high["clone_count"])
        correction = abs(fit["slope"]) / high_population
        support_pass = all(int(row["support_pass"]) and int(row["independent_runs"]) == 4 for row in rows)
        adjacent_pass = abs(adjacent_delta) <= max(0.002, 2.0 * combined_se)
        correction_pass = correction <= 0.002
        precision_pass = fit["intercept_se"] <= 0.0025
        passed = support_pass and adjacent_pass and correction_pass and precision_pass
        output.append({
            "context": context, "k": k,
            "lower_population": int(low["clone_count"]),
            "upper_population": high_population,
            "adjacent_delta": adjacent_delta,
            "combined_se": combined_se,
            "adjacent_tolerance": max(0.002, 2.0 * combined_se),
            "estimated_correction_at_upper": correction,
            "intercept_se": fit["intercept_se"],
            "support_pass": int(support_pass),
            "adjacent_pass": int(adjacent_pass),
            "correction_pass": int(correction_pass),
            "precision_pass": int(precision_pass),
            "plateau_pass": int(passed),
        })
    return output


def optional_requests(plateaus: list[dict[str, object]]) -> list[dict[str, object]]:
    requests: list[dict[str, object]] = []
    groups = [
        ("n2_known_answer", 2, [row for row in plateaus if row["context"] == "n2_t0p1"]),
        ("long_chain", 10, [row for row in plateaus if row["context"] in ("n10_final", "n10_t20_k0p3")]),
        ("long_chain", 20, [row for row in plateaus if row["context"] == "n20_final"]),
        ("long_chain", 30, [row for row in plateaus if row["context"] == "n30_final"]),
        ("long_chain", 40, [row for row in plateaus if row["context"] == "n40_final"]),
    ]
    for study, n, rows in groups:
        failed = [f"{row['context']}:k={row['k']}" for row in rows if not int(row["plateau_pass"])]
        requests.append({
            "study": study, "n": n, "run_8192": int(bool(failed)),
            "failed_member_only_gates": ";".join(failed),
        })
    return requests


def direct_reference_from_samples(stream: np.ndarray, sigma: np.ndarray,
                                  tau: float, tilts: tuple[float, ...],
                                  label: str, seed: int) -> list[dict[str, object]]:
    unique, inverse = np.unique(stream.astype(int), return_inverse=True)
    counts = np.bincount(inverse)
    output: list[dict[str, object]] = []
    for index, k in enumerate(tilts):
        log_weight = -k * sigma
        maximum = float(np.max(log_weight))
        sums = np.bincount(inverse, weights=np.exp(log_weight - maximum))
        psi = float((maximum + math.log(float(np.sum(sums)) / sigma.size)) / tau)
        rng = np.random.default_rng(seed + index * 1009)
        boot = np.empty(BOOTSTRAPS)
        for draw in range(BOOTSTRAPS):
            selected = rng.integers(0, unique.size, unique.size)
            boot[draw] = (
                maximum + math.log(float(np.sum(sums[selected])) / float(np.sum(counts[selected])))
            ) / tau
        ess = float(math.exp(2.0 * logsumexp(log_weight) - logsumexp(2.0 * log_weight)))
        output.append({
            "reference": label, "tau": tau, "k": k, "psi": psi,
            "se": float(np.std(boot, ddof=1)),
            "ci_low": float(np.quantile(boot, 0.025)),
            "ci_high": float(np.quantile(boot, 0.975)),
            "samples": sigma.size, "independent_streams": unique.size,
            "bootstrap_draws": BOOTSTRAPS, "weight_ess": ess,
        })
    return output


def n2_gauge_reference() -> tuple[list[dict[str, object]], Path]:
    path = ROOT / "raw/n2_direct_reference/n2_t0p1_blocks.csv"
    if not path.exists():
        raise RuntimeError(f"missing new n=2 direct reference {path}")
    data = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 11))
    stream = data[:, 0].astype(int)
    sigma_r = -0.4 * data[:, 1]
    return direct_reference_from_samples(
        stream, sigma_r, 0.1, (0.3, 0.5, 0.7),
        "new_n2_direct_sigma_R", 20260831001), path


def accepted_n2_total_entropy_reference() -> tuple[list[dict[str, object]], Path]:
    path = REPO / "experiments/entropy_ft_scgf_2026-08-27/total_entropy_n2_short/production/driven_blocks.csv"
    sys.path.insert(0, str(REPO / "flux"))
    from analyze_total_entropy_n2 import load_blocks  # type: ignore
    from analyze_total_entropy_parametric_n2 import crossfit_driven  # type: ignore
    blocks = load_blocks(path)
    ratio, valid, _ = crossfit_driven(blocks)
    sigma = blocks.entropy_medium[valid] + ratio[valid]
    return direct_reference_from_samples(
        blocks.stream[valid], sigma, 0.1, (0.3, 0.5, 0.7),
        "accepted_n2_parametric_total_entropy", 20260831021), path


def n10_reference() -> dict[str, object]:
    path = REPO / "experiments/long_chain_ft_2026-08-28/direct_gauge_analysis/scgf_points.csv"
    rows = [row for row in read_csv(path)
            if int(row["n"]) == 10 and close(float(row["tau"]), 20.0)
            and close(float(row["k"]), 0.3)
            and row["observable"] == "right_bath_entropy_current"]
    if len(rows) != 1:
        raise RuntimeError("cannot identify frozen n=10 direct reference")
    row = rows[0]
    low, high = float(row["psi_ci_low"]), float(row["psi_ci_high"])
    return {
        "reference": "accepted_n10_direct_sigma_R", "tau": 20.0, "k": 0.3,
        "psi": float(row["psi"]), "se": (high - low) / (2.0 * 1.96),
        "ci_low": low, "ci_high": high, "samples": int(row["sample_count"]),
        "independent_streams": int(row["n_streams"]),
        "bootstrap_draws": 300, "weight_ess": float(row["sample_weight_ess"]),
        "source_path": str(path.relative_to(REPO)), "source_sha256": sha256(path),
    }


def welch_interval(delta: float, variance_a: float, df_a: float,
                   variance_b: float, df_b: float) -> tuple[float, float, float, float]:
    variance = variance_a + variance_b
    se = math.sqrt(variance)
    denominator = variance_a * variance_a / df_a + variance_b * variance_b / df_b
    df = variance * variance / denominator if denominator > 0 else min(df_a, df_b)
    critical = float(student_t.ppf(0.975, df))
    return se, df, delta - critical * se, delta + critical * se


def extrapolated_gc_rows(fit_map: dict[tuple[str, float], dict[str, float]]) -> list[dict[str, object]]:
    pairs = (("n10_final", 0.3, 0.7), ("n20_final", 0.4, 0.6),
             ("n30_final", 0.4, 0.6), ("n40_final", 0.4, 0.6))
    rows: list[dict[str, object]] = []
    for context, low_k, high_k in pairs:
        low, high = fit_map[(context, low_k)], fit_map[(context, high_k)]
        delta = low["intercept"] - high["intercept"]
        se, df, ci_low, ci_high = welch_interval(
            delta, low["hc3_var_intercept"], low["df"],
            high["hc3_var_intercept"], high["df"])
        half = 0.5 * (ci_high - ci_low)
        passed = ci_low <= 0.0 <= ci_high and abs(delta) <= 0.005 and half <= 0.005
        rows.append({
            "context": context, "k": low_k, "one_minus_k": high_k,
            "psi_inf_k": low["intercept"], "psi_inf_one_minus_k": high["intercept"],
            "extrapolated_residual": delta, "residual_se": se, "df": df,
            "ci_low": ci_low, "ci_high": ci_high, "ci_half_width": half,
            "contains_zero": int(ci_low <= 0.0 <= ci_high),
            "absolute_gate": int(abs(delta) <= 0.005),
            "precision_gate": int(half <= 0.005), "gc_gate": int(passed),
        })
    return rows


def reference_gap_rows(levels: list[dict[str, object]],
                       fit_map: dict[tuple[str, float], dict[str, float]],
                       references: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    reference_map = {(str(row["reference"]), float(row["k"])): row for row in references}
    specifications = [
        ("n2_t0p1", "new_n2_direct_sigma_R", k, 0.005) for k in (0.3, 0.5, 0.7)
    ] + [("n10_t20_k0p3", "accepted_n10_direct_sigma_R", 0.3, 0.003)]
    by_population: list[dict[str, object]] = []
    extrapolated: list[dict[str, object]] = []
    for context, reference_name, k, margin in specifications:
        reference = reference_map[(reference_name, k)]
        matching = [row for row in levels if row["context"] == context and close(float(row["k"]), k)]
        for row in sorted(matching, key=lambda value: int(value["clone_count"])):
            gap = float(row["mean_psi"]) - float(reference["psi"])
            se = math.sqrt(float(row["run_se"]) ** 2 + float(reference["se"]) ** 2)
            by_population.append({
                "context": context, "k": k, "clone_count": row["clone_count"],
                "cloning_mean": row["mean_psi"], "cloning_run_se": row["run_se"],
                "direct_psi": reference["psi"], "direct_se": reference["se"],
                "gap": gap, "gap_se": se, "ci_low": gap - 1.96 * se,
                "ci_high": gap + 1.96 * se,
            })
        fit = fit_map[(context, k)]
        gap = fit["intercept"] - float(reference["psi"])
        variance = fit["hc3_var_intercept"] + float(reference["se"]) ** 2
        se = math.sqrt(variance)
        ci_low, ci_high = gap - 1.96 * se, gap + 1.96 * se
        passed = ci_low <= 0.0 <= ci_high and abs(gap) <= margin
        extrapolated.append({
            "context": context, "k": k, "cloning_intercept": fit["intercept"],
            "direct_psi": reference["psi"], "extrapolated_gap": gap,
            "gap_se": se, "ci_low": ci_low, "ci_high": ci_high,
            "contains_zero": int(ci_low <= 0.0 <= ci_high),
            "absolute_margin": margin, "absolute_gate": int(abs(gap) <= margin),
            "reference_gate": int(passed),
        })
    return by_population, extrapolated


def generate_command_audit(records: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = []
    for record in records:
        path = REPO / str(record["command_path"])
        if not path.exists():
            raise RuntimeError(f"missing command record {path}")
        rows.append({
            "study": record["study"], "stage": record["stage"], "n": record["n"],
            "clone_count": record["clone_count"], "k": record["k"],
            "run": record["run"], "seed": record["seed"],
            "command": path.read_text().strip(),
        })
    direct = ROOT / "raw/n2_direct_reference/n2_t0p1.command.txt"
    if direct.exists():
        rows.insert(0, {"study": "n2_direct_reference", "stage": "primary", "n": 2,
                        "clone_count": "", "k": "", "run": 1, "seed": 210000001,
                        "command": direct.read_text().strip()})
    write_csv(ANALYSIS / "executed_commands.csv", rows,
              ["study", "stage", "n", "clone_count", "k", "run", "seed", "command"])


def raw_manifest(records: list[dict[str, object]], extra: list[Path]) -> None:
    paths = []
    for record in records:
        paths.extend((REPO / str(record["summary_path"]), REPO / str(record["timeseries_path"])))
    paths.extend(extra)
    with (ANALYSIS / "RAW_MANIFEST.sha256").open("w") as handle:
        for path in sorted(set(paths)):
            handle.write(f"{sha256(path)}  {path.relative_to(REPO)}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("plateau", "final"), required=True)
    args = parser.parse_args()
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    records = collect_runs(require_primary=True)
    raw_fields = [
        "study", "stage", "n", "clone_count", "k", "run", "seed", "horizon",
        "selection_time", "dt", "psi_final", "psi_t20", "minimum_weight_ess",
        "minimum_weight_ess_fraction", "minimum_unique_roots",
        "minimum_root_weight_ess", "midpoint_failures", "support_pass",
        "summary_path", "timeseries_path", "command_path",
    ]
    write_csv(ANALYSIS / "raw_psi_per_seed.csv", records, raw_fields)
    levels, fits, fit_map = build_population_tables(records)
    level_fields = ["context", "k", "clone_count", "independent_runs", "mean_psi", "run_se", "support_pass"]
    fit_fields = ["context", "k", "runs", "population_levels", "intercept", "slope",
                  "intercept_se", "slope_se", "intercept_ci_low", "intercept_ci_high",
                  "df", "r2", "hc3_var_intercept", "max_abs_residual"]
    write_csv(ANALYSIS / "population_levels.csv", levels, level_fields)
    write_csv(ANALYSIS / "population_member_fits.csv", fits, fit_fields)
    plateaus = plateau_rows(levels, fit_map)
    plateau_fields = [
        "context", "k", "lower_population", "upper_population", "adjacent_delta",
        "combined_se", "adjacent_tolerance", "estimated_correction_at_upper",
        "intercept_se", "support_pass", "adjacent_pass", "correction_pass",
        "precision_pass", "plateau_pass",
    ]
    write_csv(ANALYSIS / "population_plateau_gates.csv", plateaus, plateau_fields)
    requests = optional_requests(plateaus)
    request_name = (
        "optional_population_requests_primary.csv"
        if args.stage == "plateau" else "optional_population_requests_final.csv"
    )
    write_csv(ANALYSIS / request_name, requests,
              ["study", "n", "run_8192", "failed_member_only_gates"])

    if args.stage == "plateau":
        print(json.dumps({"stage": "plateau", "runs": len(records), "requests": requests}, indent=2))
        return

    for request in requests:
        if int(request["run_8192"]):
            expected = [r for r in records if int(r["n"]) == int(request["n"])
                        and r["study"] == request["study"] and int(r["clone_count"]) == 8192]
            expected_count = 12 if request["study"] == "n2_known_answer" else 8
            if len(expected) != expected_count:
                raise RuntimeError(f"required N=8192 group incomplete for {request}")

    n2_refs, n2_path = n2_gauge_reference()
    for row in n2_refs:
        row["source_path"] = str(n2_path.relative_to(REPO))
        row["source_sha256"] = sha256(n2_path)
    total_refs, total_path = accepted_n2_total_entropy_reference()
    for row in total_refs:
        row["source_path"] = str(total_path.relative_to(REPO))
        row["source_sha256"] = sha256(total_path)
    n10_ref = n10_reference()
    references = [*n2_refs, n10_ref]
    reference_fields = ["reference", "tau", "k", "psi", "se", "ci_low", "ci_high",
                        "samples", "independent_streams", "bootstrap_draws", "weight_ess",
                        "source_path", "source_sha256"]
    write_csv(ANALYSIS / "same_observable_direct_references.csv", references, reference_fields)
    write_csv(ANALYSIS / "accepted_total_entropy_reference_separate.csv", total_refs, reference_fields)

    gc_rows = extrapolated_gc_rows(fit_map)
    gc_fields = [
        "context", "k", "one_minus_k", "psi_inf_k", "psi_inf_one_minus_k",
        "extrapolated_residual", "residual_se", "df", "ci_low", "ci_high",
        "ci_half_width", "contains_zero", "absolute_gate", "precision_gate", "gc_gate",
    ]
    write_csv(ANALYSIS / "extrapolated_gc_residuals.csv", gc_rows, gc_fields)
    gap_levels, gap_intercepts = reference_gap_rows(levels, fit_map, references)
    write_csv(ANALYSIS / "reference_gaps_by_population.csv", gap_levels,
              ["context", "k", "clone_count", "cloning_mean", "cloning_run_se",
               "direct_psi", "direct_se", "gap", "gap_se", "ci_low", "ci_high"])
    write_csv(ANALYSIS / "extrapolated_reference_gaps.csv", gap_intercepts,
              ["context", "k", "cloning_intercept", "direct_psi", "extrapolated_gap",
               "gap_se", "ci_low", "ci_high", "contains_zero", "absolute_margin",
               "absolute_gate", "reference_gate"])
    generate_command_audit(records)
    raw_manifest(records, [n2_path, total_path])

    final_plateau_pass = all(int(row["plateau_pass"]) for row in plateaus)
    reference_pass = all(int(row["reference_gate"]) for row in gap_intercepts)
    symmetry_pass = all(int(row["gc_gate"]) for row in gc_rows)
    support_pass = all(int(record["support_pass"]) for record in records)
    overall = final_plateau_pass and reference_pass and symmetry_pass and support_pass
    verdict = "ESTABLISHED_WITHIN_FROZEN_NUMERICAL_SCOPE" if overall else "NOT_NUMERICALLY_ESTABLISHED"
    commit = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
    metadata = {
        "verdict": verdict,
        "git_commit_at_analysis": commit,
        "source_sha256": sha256(ROOT / "source_archive/NLS_entropy_cloning.cpp"),
        "binary_sha256": sha256(ROOT / "bin/entropy_cloning_v2"),
        "support_pass": support_pass,
        "population_plateau_pass": final_plateau_pass,
        "physical_known_answer_and_driven_crosscheck_pass": reference_pass,
        "extrapolated_gc_pass": symmetry_pass,
        "raw_runs": len(records),
        "claim_boundary": "finite n, frozen tilts/horizons only; not a theorem or n->infinity result",
        "observable_note": "cloning gaps use Sigma_R references; accepted finite-time total entropy is separate",
    }
    (ANALYSIS / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    lines = [
        "# Final cloning-bias validation verdict", "", f"**{verdict}**", "",
        f"- support/integrity: {'PASS' if support_pass else 'FAIL'}",
        f"- population plateau: {'PASS' if final_plateau_pass else 'FAIL'}",
        f"- physical known-answer plus n=10 direct gap: {'PASS' if reference_pass else 'FAIL'}",
        f"- extrapolated GC residuals: {'PASS' if symmetry_pass else 'FAIL'}", "",
        "The finite-time accepted n=2 total-entropy result is reported separately and is not",
        "subtracted from the Sigma_R cloning estimator.  A failed item leaves the long-chain",
        "GC claim numerically unresolved under the frozen protocol.", "",
    ]
    (ANALYSIS / "FINAL_VERDICT.md").write_text("\n".join(lines))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
