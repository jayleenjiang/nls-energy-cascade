#!/usr/bin/env python3
"""Frozen multi-affinity timestep and equilibrium estimator-floor analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TIMES = (5, 10, 20, 25, 40, 80, 160, 320, 640)
SELECTIONS = ("all_admissible", "tau_ge_10")
MIN_NEGATIVE = 500
MIN_SYMMETRIC_PAIRS = 10
MIN_BOOTSTRAP = 800
BOOTSTRAP = 1000
STREAMS = 128
BLOCKS_PER_STREAM = 31252
REFERENCE_DBETA = 0.05714285714285716
SHORT_ANALYSIS_SHA = "100bc49fc4849d3e5f248600a7bb4000af6aadc9161415e1ce2b301de6f9863c"
POOLED_ANALYSIS_SHA = "dc9ed26c662a327459486f245ecf57e245cf4eb407601027de5684afc68d9730"
MID_REFERENCE_SHA = "91f4413cb2275fcda19e20cadeaa18e88ed113f5e922bdafdf3552db2078d159"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(name: str, path: Path, expected_sha: str):
    actual = sha256(path)
    if actual != expected_sha:
        raise RuntimeError(f"{path}: expected SHA {expected_sha}, found {actual}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and fields is None:
        raise ValueError(f"empty rows require fields: {path}")
    if fields is None:
        fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: np.ndarray) -> tuple[int, float, float, float]:
    valid = np.asarray(values, dtype=np.float64)
    valid = valid[np.isfinite(valid)]
    if valid.size < MIN_BOOTSTRAP:
        return int(valid.size), math.nan, math.nan, math.nan
    low, high = np.percentile(valid, [2.5, 97.5])
    return int(valid.size), float(np.std(valid, ddof=1)), float(low), float(high)


def load_cases(path: Path, root: Path, has_affinity: bool) -> list[dict]:
    result = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            item = dict(row)
            if not has_affinity:
                item["affinity"] = "equilibrium"
            for key in ("T_left", "T_right", "delta_beta", "dt"):
                item[key] = float(item[key])
            for key in ("seed", "bootstrap_seed"):
                item[key] = int(item[key])
            item["path"] = (root / item["blocks_file"]).resolve()
            result.append(item)
    return result


def weighted_operator(x: np.ndarray, se: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack((np.ones_like(x), x))
    weight = 1.0 / np.square(se)
    operator = np.linalg.solve(
        design.T @ (weight[:, None] * design), design.T * weight
    )
    return design, operator


def fit_dt_dependence(
    affinity: str,
    selection: str,
    cases: list[dict],
    products: dict[str, dict],
) -> tuple[dict, list[dict]]:
    ordered = sorted(cases, key=lambda item: item["dt"], reverse=True)
    dts = np.asarray([case["dt"] for case in ordered], dtype=np.float64)
    summaries = [products[case["case"]]["summaries"][selection] for case in ordered]
    values = np.asarray([row["a_inf_over_delta_beta"] for row in summaries], dtype=float)
    times_used = [row["times_used"] for row in summaries]
    all_finite = bool(np.all(np.isfinite(values)))
    all_adequate = bool(all(int(row["linear_adequate"]) == 1 for row in summaries))
    if not all_finite:
        row = {
            "affinity": affinity,
            "selection": selection,
            "timesteps": ";".join(f"{value:.8g}" for value in dts),
            "ratios": ";".join(str(value) for value in values),
            "times_used_by_dt": "|".join(
                f"{dt:.8g}:{times}" for dt, times in zip(dts, times_used)
            ),
            "same_times_used": int(len(set(times_used)) == 1),
            "all_finite_time_fits_adequate": int(all_adequate),
            "bootstrap_accepted": 0,
            "bootstrap_total": BOOTSTRAP,
            "dt_zero_ratio": math.nan,
            "dt_zero_ci_low": math.nan,
            "dt_zero_ci_high": math.nan,
            "dt_slope": math.nan,
            "dt_slope_ci_low": math.nan,
            "dt_slope_ci_high": math.nan,
            "linear_reduced_chi2": math.nan,
            "max_abs_standardized_residual": math.nan,
            "linear_Odt_compatible": 0,
            "dt_zero_ci_contains_one": "",
            "interpretable": 0,
        }
        return row, []

    boot_matrix = np.vstack([
        np.asarray([
            row["a_inf_over_delta_beta"]
            for row in products[case["case"]]["extrapolation_bootstrap"][selection]
        ], dtype=float)
        for case in ordered
    ])
    ses = np.asarray([
        np.std(row[np.isfinite(row)], ddof=1) for row in boot_matrix
    ])
    design, operator = weighted_operator(dts, ses)
    beta = operator @ values
    fitted = design @ beta
    standardized = (values - fitted) / ses
    chi2 = float(np.sum(np.square(standardized)))
    reduced = chi2 / (len(dts) - 2)
    max_abs = float(np.max(np.abs(standardized)))
    valid = np.all(np.isfinite(boot_matrix), axis=0)
    beta_boot = np.full((2, BOOTSTRAP), np.nan)
    beta_boot[:, valid] = operator @ boot_matrix[:, valid]
    accepted, _, intercept_low, intercept_high = percentile(beta_boot[0])
    _, _, slope_low, slope_high = percentile(beta_boot[1])
    compatible = reduced <= 2.0 and max_abs <= 3.0
    same_times = len(set(times_used)) == 1
    row = {
        "affinity": affinity,
        "selection": selection,
        "timesteps": ";".join(f"{value:.8g}" for value in dts),
        "ratios": ";".join(f"{value:.17g}" for value in values),
        "times_used_by_dt": "|".join(
            f"{dt:.8g}:{times}" for dt, times in zip(dts, times_used)
        ),
        "same_times_used": int(same_times),
        "all_finite_time_fits_adequate": int(all_adequate),
        "bootstrap_accepted": accepted,
        "bootstrap_total": BOOTSTRAP,
        "dt_zero_ratio": float(beta[0]),
        "dt_zero_ci_low": intercept_low,
        "dt_zero_ci_high": intercept_high,
        "dt_slope": float(beta[1]),
        "dt_slope_ci_low": slope_low,
        "dt_slope_ci_high": slope_high,
        "linear_reduced_chi2": reduced,
        "max_abs_standardized_residual": max_abs,
        "linear_Odt_compatible": int(compatible),
        "dt_zero_ci_contains_one": int(intercept_low <= 1.0 <= intercept_high),
        "interpretable": int(all_adequate and same_times and accepted >= MIN_BOOTSTRAP),
    }
    boot_rows = [{
        "affinity": affinity,
        "selection": selection,
        "replicate": replicate,
        "valid": int(valid[replicate]),
        "dt_zero_ratio": beta_boot[0, replicate],
        "dt_slope": beta_boot[1, replicate],
    } for replicate in range(BOOTSTRAP)]
    return row, boot_rows


def analyze_case(case: dict, short, pooled) -> tuple[dict, list[dict], list[dict]]:
    matrix = short.load_heat_matrix(case["path"], BLOCKS_PER_STREAM)
    rng = np.random.default_rng(case["bootstrap_seed"])
    multiplicities = rng.multinomial(
        STREAMS, np.full(STREAMS, 1.0 / STREAMS), size=BOOTSTRAP
    )
    per_time: dict[int, dict] = {}
    per_rows: list[dict] = []
    symmetric_rows: list[dict] = []
    for tau in TIMES:
        by_stream = short.aggregate(matrix, tau)
        values = by_stream.ravel()
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        standard = (values - mean) / sd
        fit = short.matched_fit(values)
        slope_boot, intercept_boot = short.bootstrap_fit(
            by_stream, fit, multiplicities
        )
        accepted, boot_se, ci_low, ci_high = percentile(slope_boot)
        n_negative = int(np.count_nonzero(values < 0.0))
        pairs = int(fit.positive_k.size)
        admissible = bool(
            fit.resolved and n_negative >= MIN_NEGATIVE
            and pairs >= MIN_SYMMETRIC_PAIRS
        )
        reference = case["delta_beta"]
        row = {
            "affinity": case["affinity"],
            "case": case["case"],
            "kind": case["kind"],
            "T_left": case["T_left"],
            "T_right": case["T_right"],
            "delta_beta": reference,
            "dt": case["dt"],
            "seed": case["seed"],
            "bootstrap_seed": case["bootstrap_seed"],
            "tau": tau,
            "inverse_tau": 1.0 / tau,
            "N_windows": int(values.size),
            "mean_Q": mean,
            "std_Q": sd,
            "skew_Q": float(np.mean(standard ** 3)),
            "excess_kurtosis_Q": float(np.mean(standard ** 4) - 3.0),
            "n_negative": n_negative,
            "p_negative": n_negative / values.size,
            "dQ": float(fit.dx),
            "n_symmetric_pairs": pairs,
            "full_sample_resolved": int(fit.resolved),
            "admissible": int(admissible),
            "a_fit": float(fit.slope),
            "a_bootstrap_se": boot_se,
            "a_ci_low": ci_low,
            "a_ci_high": ci_high,
            "ci_contains_reference": (
                int(ci_low <= reference <= ci_high)
                if math.isfinite(ci_low) else ""
            ),
            "a_over_delta_beta": (
                float(fit.slope / reference)
                if reference and fit.resolved else math.nan
            ),
            "intercept": float(fit.intercept),
            "R2": float(fit.r2),
            "quadratic_linear": float(fit.quadratic_linear),
            "quadratic_coefficient": float(fit.quadratic),
            "quadratic_intercept": float(fit.quadratic_intercept),
            "quadratic_R2": float(fit.quadratic_r2),
            "bootstrap_accepted": accepted,
            "bootstrap_total": BOOTSTRAP,
        }
        per_rows.append(row)
        per_time[tau] = {
            "row": row,
            "slope_boot": slope_boot,
            "intercept_boot": intercept_boot,
        }
        if fit.resolved:
            for k, pos, neg in zip(
                fit.positive_k, fit.positive_indices, fit.negative_indices
            ):
                plus = int(fit.counts[pos])
                minus = int(fit.counts[neg])
                symmetric_rows.append({
                    "affinity": case["affinity"],
                    "case": case["case"],
                    "dt": case["dt"],
                    "tau": tau,
                    "k": int(k),
                    "Q_center": float(k * fit.dx),
                    "count_plus": plus,
                    "count_minus": minus,
                    "log_ratio": float(math.log(plus / minus)),
                    "log_ratio_se": float(math.sqrt(1.0 / plus + 1.0 / minus)),
                })
    del matrix

    summaries = {}
    summary_bootstrap = {}
    residual_rows = []
    extrapolation_rows = []
    extrapolation_bootstrap_rows = []
    admissible = [tau for tau in TIMES if per_time[tau]["row"]["admissible"]]
    selections = {
        "all_admissible": admissible,
        "tau_ge_10": [tau for tau in admissible if tau >= 10],
    }
    for selection, selected in selections.items():
        summary, residual, bootstrap = pooled.fit_selection(
            case["case"], case["delta_beta"], selection, selected, per_time
        )
        summary.update({
            "affinity": case["affinity"],
            "kind": case["kind"],
            "dt": case["dt"],
            "T_left": case["T_left"],
            "T_right": case["T_right"],
            "seed": case["seed"],
        })
        for item in residual:
            item.update({"affinity": case["affinity"], "dt": case["dt"]})
        for item in bootstrap:
            item.update({"affinity": case["affinity"], "dt": case["dt"]})
        summaries[selection] = summary
        summary_bootstrap[selection] = bootstrap
        extrapolation_rows.append(summary)
        residual_rows.extend(residual)
        extrapolation_bootstrap_rows.extend(bootstrap)
    product = {
        "case": case,
        "per_time": per_time,
        "summaries": summaries,
        "extrapolation_bootstrap": summary_bootstrap,
        "extrapolation_rows": extrapolation_rows,
        "residual_rows": residual_rows,
        "extrapolation_bootstrap_rows": extrapolation_bootstrap_rows,
    }
    return product, per_rows, symmetric_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_dir.resolve()
    analysis_dir = root / "analysis"
    figures_dir = root / "figures"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    baseline_dir = root.parent / "entropy_ft_window_length_2026-09-10"
    pooled_dir = root.parent / "entropy_ft_pooled_extrapolation_2026-09-10"
    timestep_dir = root.parent / "entropy_ft_timestep_bias_2026-09-10"
    short_path = baseline_dir / "analyze_window_dependence.py"
    pooled_path = pooled_dir / "analyze_pooled_extrapolation.py"
    mid_reference_path = timestep_dir / "analysis/extrapolation_summary.csv"
    if sha256(mid_reference_path) != MID_REFERENCE_SHA:
        raise RuntimeError("completed (7,5) reference summary hash changed")
    short = load_module("frozen_short_multiaffinity", short_path, SHORT_ANALYSIS_SHA)
    pooled = load_module("frozen_pooled_multiaffinity", pooled_path, POOLED_ANALYSIS_SHA)

    driven_cases = load_cases(root / "CASES.tsv", root, True)
    equilibrium_cases = load_cases(root / "EQUILIBRIUM_CASES.tsv", root, False)
    all_cases = driven_cases + equilibrium_cases
    products: dict[str, dict] = {}
    per_tau_rows: list[dict] = []
    symmetric_rows: list[dict] = []
    slope_bootstrap_rows: list[dict] = []
    extrapolation_rows: list[dict] = []
    extrapolation_bootstrap_rows: list[dict] = []
    residual_rows: list[dict] = []
    input_rows: list[dict] = []

    for case in all_cases:
        print(f"loading {case['case']}: {case['path']}", flush=True)
        product, case_rows, case_symmetric = analyze_case(case, short, pooled)
        products[case["case"]] = product
        per_tau_rows.extend(case_rows)
        symmetric_rows.extend(case_symmetric)
        extrapolation_rows.extend(product["extrapolation_rows"])
        extrapolation_bootstrap_rows.extend(product["extrapolation_bootstrap_rows"])
        residual_rows.extend(product["residual_rows"])
        for replicate in range(BOOTSTRAP):
            row = {
                "affinity": case["affinity"],
                "case": case["case"],
                "dt": case["dt"],
                "replicate": replicate,
            }
            for tau in TIMES:
                row[f"a_tau_{tau}"] = product["per_time"][tau]["slope_boot"][replicate]
            slope_bootstrap_rows.append(row)
        input_rows.append({
            "affinity": case["affinity"],
            "case": case["case"],
            "kind": case["kind"],
            "path": str(case["path"]),
            "rows": STREAMS * BLOCKS_PER_STREAM,
            "sha256": sha256(case["path"]),
        })

    write_csv(analysis_dir / "per_tau_results.csv", per_tau_rows)
    write_csv(analysis_dir / "symmetric_bin_raw_counts.csv", symmetric_rows)
    write_csv(analysis_dir / "slope_bootstrap.csv", slope_bootstrap_rows)
    write_csv(analysis_dir / "extrapolation_summary.csv", extrapolation_rows)
    write_csv(analysis_dir / "extrapolation_bootstrap.csv", extrapolation_bootstrap_rows)
    write_csv(analysis_dir / "fit_residuals.csv", residual_rows)
    write_csv(analysis_dir / "input_hashes.csv", input_rows)

    contrast_rows = []
    dt_rows = []
    dt_bootstrap_rows = []
    by_affinity = {
        affinity: [case for case in driven_cases if case["affinity"] == affinity]
        for affinity in ("weak", "moderate")
    }
    for affinity, cases in by_affinity.items():
        ordered = sorted(cases, key=lambda item: item["dt"], reverse=True)
        base = ordered[0]
        for fine in ordered[1:]:
            for tau in TIMES:
                base_item = products[base["case"]]["per_time"][tau]
                fine_item = products[fine["case"]]["per_time"][tau]
                diff_boot = fine_item["slope_boot"] - base_item["slope_boot"]
                accepted, diff_se, low, high = percentile(diff_boot)
                base_mean = base_item["row"]["mean_Q"]
                fine_mean = fine_item["row"]["mean_Q"]
                contrast_rows.append({
                    "affinity": affinity,
                    "quantity": "per_tau_slope",
                    "selection": "",
                    "tau": tau,
                    "baseline_case": base["case"],
                    "finer_case": fine["case"],
                    "baseline_dt": base["dt"],
                    "finer_dt": fine["dt"],
                    "baseline_value": base_item["row"]["a_fit"],
                    "finer_value": fine_item["row"]["a_fit"],
                    "finer_minus_baseline": fine_item["row"]["a_fit"] - base_item["row"]["a_fit"],
                    "difference_bootstrap_se": diff_se,
                    "difference_ci_low": low,
                    "difference_ci_high": high,
                    "difference_ci_contains_zero": int(low <= 0 <= high),
                    "bootstrap_accepted": accepted,
                    "mean_Q_baseline": base_mean,
                    "mean_Q_finer": fine_mean,
                    "mean_Q_relative_change": (fine_mean - base_mean) / base_mean,
                })
            for selection in SELECTIONS:
                base_summary = products[base["case"]]["summaries"][selection]
                fine_summary = products[fine["case"]]["summaries"][selection]
                base_boot = np.asarray([
                    row["a_inf_over_delta_beta"]
                    for row in products[base["case"]]["extrapolation_bootstrap"][selection]
                ])
                fine_boot = np.asarray([
                    row["a_inf_over_delta_beta"]
                    for row in products[fine["case"]]["extrapolation_bootstrap"][selection]
                ])
                accepted, diff_se, low, high = percentile(fine_boot - base_boot)
                contrast_rows.append({
                    "affinity": affinity,
                    "quantity": "intercept_ratio",
                    "selection": selection,
                    "tau": "",
                    "baseline_case": base["case"],
                    "finer_case": fine["case"],
                    "baseline_dt": base["dt"],
                    "finer_dt": fine["dt"],
                    "baseline_value": base_summary["a_inf_over_delta_beta"],
                    "finer_value": fine_summary["a_inf_over_delta_beta"],
                    "finer_minus_baseline": float(fine_summary["a_inf_over_delta_beta"]) - float(base_summary["a_inf_over_delta_beta"]),
                    "difference_bootstrap_se": diff_se,
                    "difference_ci_low": low,
                    "difference_ci_high": high,
                    "difference_ci_contains_zero": int(low <= 0 <= high) if math.isfinite(low) else "",
                    "bootstrap_accepted": accepted,
                    "mean_Q_baseline": "",
                    "mean_Q_finer": "",
                    "mean_Q_relative_change": "",
                })
        for selection in SELECTIONS:
            row, boot = fit_dt_dependence(affinity, selection, cases, products)
            dt_rows.append(row)
            dt_bootstrap_rows.extend(boot)
    write_csv(analysis_dir / "timestep_contrasts.csv", contrast_rows)
    write_csv(analysis_dir / "dt_zero_extrapolation.csv", dt_rows)
    write_csv(
        analysis_dir / "dt_zero_bootstrap.csv",
        dt_bootstrap_rows,
        fields=["affinity", "selection", "replicate", "valid", "dt_zero_ratio", "dt_slope"],
    )

    eq_rows = []
    eq_boot_values: dict[str, np.ndarray] = {}
    for case in equilibrium_cases:
        product = products[case["case"]]
        slopes = np.asarray([product["per_time"][tau]["row"]["a_fit"] for tau in TIMES])
        boot = np.vstack([product["per_time"][tau]["slope_boot"] for tau in TIMES])
        ses = np.asarray([
            np.std(row[np.isfinite(row)], ddof=1) for row in boot
        ])
        weights = 1.0 / np.square(ses)
        combined = float(np.sum(weights * slopes) / np.sum(weights))
        combined_boot = np.sum(weights[:, None] * boot, axis=0) / np.sum(weights)
        accepted, combined_se, low, high = percentile(combined_boot)
        eq_boot_values[case["case"]] = combined_boot
        n_positive = int(np.count_nonzero(slopes > 0))
        n_negative = int(np.count_nonzero(slopes < 0))
        eq_rows.append({
            "case": case["case"],
            "dt": case["dt"],
            "weighted_mean_slope": combined,
            "bootstrap_se": combined_se,
            "ci_low": low,
            "ci_high": high,
            "ci_contains_zero": int(low <= 0 <= high),
            "z_from_zero": combined / combined_se,
            "positive_windows": n_positive,
            "negative_windows": n_negative,
            "zero_windows": len(TIMES) - n_positive - n_negative,
            "sign_pattern_tau_order": ";".join("+" if value > 0 else "-" if value < 0 else "0" for value in slopes),
            "fraction_of_dbeta_0p057143": combined / REFERENCE_DBETA,
            "bootstrap_accepted": accepted,
            "bootstrap_total": BOOTSTRAP,
            "weight_rule": "fixed_inverse_bootstrap_variance",
        })
    eq_rows.sort(key=lambda row: float(row["dt"]), reverse=True)
    coarse_case = next(case for case in equilibrium_cases if case["dt"] == 0.0005)
    fine_case = next(case for case in equilibrium_cases if case["dt"] == 0.00025)
    eq_diff_boot = eq_boot_values[fine_case["case"]] - eq_boot_values[coarse_case["case"]]
    accepted, diff_se, diff_low, diff_high = percentile(eq_diff_boot)
    coarse_row = next(row for row in eq_rows if row["case"] == coarse_case["case"])
    fine_row = next(row for row in eq_rows if row["case"] == fine_case["case"])
    eq_difference = [{
        "coarse_case": coarse_case["case"],
        "fine_case": fine_case["case"],
        "fine_minus_coarse": float(fine_row["weighted_mean_slope"]) - float(coarse_row["weighted_mean_slope"]),
        "bootstrap_se": diff_se,
        "ci_low": diff_low,
        "ci_high": diff_high,
        "ci_contains_zero": int(diff_low <= 0 <= diff_high),
        "bootstrap_accepted": accepted,
        "bootstrap_total": BOOTSTRAP,
    }]
    write_csv(analysis_dir / "equilibrium_weighted_baseline.csv", eq_rows)
    write_csv(analysis_dir / "equilibrium_baseline_difference.csv", eq_difference)

    mid_reference_rows = read_csv(mid_reference_path)
    mid_finest = {
        row["selection"]: row
        for row in mid_reference_rows
        if row["case"] == "driven_dt1p25e4"
    }
    affinity_rows = []
    for selection in SELECTIONS:
        for affinity in ("weak", "moderate"):
            finest_case = min(by_affinity[affinity], key=lambda item: item["dt"])
            summary = products[finest_case["case"]]["summaries"][selection]
            adequate = int(summary["linear_adequate"]) == 1
            contains = int(summary["reference_in_ci"]) == 1
            affinity_rows.append({
                "affinity": affinity,
                "selection": selection,
                "source": finest_case["case"],
                "dt": finest_case["dt"],
                "a_inf_over_delta_beta": summary["a_inf_over_delta_beta"],
                "ratio_ci_low": summary["ratio_ci_low"],
                "ratio_ci_high": summary["ratio_ci_high"],
                "linear_reduced_chi2": summary["linear_reduced_chi2"],
                "max_abs_standardized_residual": summary["max_abs_standardized_residual"],
                "adequate": int(adequate),
                "contains_one": int(contains),
                "finest_step_pass": int(adequate and contains),
            })
        mid = mid_finest[selection]
        adequate = int(mid["linear_adequate"]) == 1
        contains = int(mid["reference_in_ci"]) == 1
        affinity_rows.append({
            "affinity": "mid",
            "selection": selection,
            "source": "entropy_ft_timestep_bias_2026-09-10/driven_dt1p25e4",
            "dt": mid["dt"],
            "a_inf_over_delta_beta": mid["a_inf_over_delta_beta"],
            "ratio_ci_low": mid["ratio_ci_low"],
            "ratio_ci_high": mid["ratio_ci_high"],
            "linear_reduced_chi2": mid["linear_reduced_chi2"],
            "max_abs_standardized_residual": mid["max_abs_standardized_residual"],
            "adequate": int(adequate),
            "contains_one": int(contains),
            "finest_step_pass": int(adequate and contains),
        })
    write_csv(analysis_dir / "three_affinity_finest_summary.csv", affinity_rows)

    nominal_floor = max(abs(float(row["weighted_mean_slope"])) for row in eq_rows)
    mid_all = mid_finest["all_admissible"]
    mid_residual = abs((1.0 - float(mid_all["a_inf_over_delta_beta"])) * REFERENCE_DBETA)
    floor_ratio = nominal_floor / mid_residual
    primary_affinity = [row for row in affinity_rows if row["selection"] == "all_admissible"]
    all_three = all(int(row["finest_step_pass"]) == 1 for row in primary_affinity)
    moderate = next(row for row in primary_affinity if row["affinity"] == "moderate")
    if all_three:
        verdict = "THREE_AFFINITIES_NUMERICALLY_CONSISTENT_WITH_FR"
    elif int(moderate["adequate"]) == 0:
        verdict = "MODERATE_AFFINITY_ASYMPTOTE_UNRESOLVED"
    else:
        verdict = "AT_LEAST_ONE_FINEST_AFFINITY_EXCLUDES_FR"
    verdict_payload = {
        "verdict": verdict,
        "all_three_primary_finest_pass": all_three,
        "nominal_estimator_floor_absolute": nominal_floor,
        "nominal_estimator_floor_fraction_of_dbeta_0p057143": nominal_floor / REFERENCE_DBETA,
        "mid_affinity_finest_residual_absolute": mid_residual,
        "floor_to_mid_residual_ratio": floor_ratio,
        "same_order_by_predeclared_factor_three_rule": bool(1.0 / 3.0 <= floor_ratio <= 3.0),
        "equilibrium_baselines_statistically_consistent": bool(eq_difference[0]["ci_contains_zero"]),
        "primary_rows": primary_affinity,
    }
    (analysis_dir / "VERDICT.json").write_text(json.dumps(verdict_payload, indent=2) + "\n")

    colors = {0.0005: "#0072B2", 0.00025: "#E69F00", 0.000125: "#009E73"}
    for affinity, cases in by_affinity.items():
        fig, ax = plt.subplots(figsize=(7.3, 4.7))
        for case in sorted(cases, key=lambda item: item["dt"], reverse=True):
            rows = [products[case["case"]]["per_time"][tau]["row"] for tau in TIMES]
            admissible = [row for row in rows if row["admissible"]]
            ax.errorbar(
                [1.0 / row["tau"] for row in admissible],
                [row["a_fit"] for row in admissible],
                yerr=[
                    [row["a_fit"] - row["a_ci_low"] for row in admissible],
                    [row["a_ci_high"] - row["a_fit"] for row in admissible],
                ],
                marker="o", capsize=2.5, color=colors[case["dt"]],
                label=fr"$dt={case['dt']:.2g}$",
            )
            summary = products[case["case"]]["summaries"]["all_admissible"]
            if math.isfinite(float(summary["a_inf"])):
                xline = np.linspace(0, max(1.0 / row["tau"] for row in admissible), 150)
                ax.plot(xline, float(summary["a_inf"]) + float(summary["c"]) * xline,
                        color=colors[case["dt"]], linewidth=1)
        ax.axhline(cases[0]["delta_beta"], color="black", linestyle="--", linewidth=1)
        ax.set_xlabel(r"$1/\tau$")
        ax.set_ylabel(r"$a(\tau)$")
        ax.grid(alpha=0.22)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(figures_dir / f"finite_time_{affinity}.pdf", bbox_inches="tight")
        fig.savefig(figures_dir / f"finite_time_{affinity}.png", dpi=240, bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    for affinity, marker in (("weak", "o"), ("moderate", "s")):
        cases = sorted(by_affinity[affinity], key=lambda item: item["dt"], reverse=True)
        rows = [products[case["case"]]["summaries"]["all_admissible"] for case in cases]
        axes[0].errorbar(
            [case["dt"] for case in cases],
            [float(row["a_inf_over_delta_beta"]) for row in rows],
            yerr=[
                [float(row["a_inf_over_delta_beta"]) - float(row["ratio_ci_low"]) for row in rows],
                [float(row["ratio_ci_high"]) - float(row["a_inf_over_delta_beta"]) for row in rows],
            ], marker=marker, capsize=3, label=affinity,
        )
        base_case = max(cases, key=lambda item: item["dt"])
        for case in cases:
            for tau in (5, 25):
                base_mean = products[base_case["case"]]["per_time"][tau]["row"]["mean_Q"]
                mean = products[case["case"]]["per_time"][tau]["row"]["mean_Q"]
                axes[1].plot(case["dt"], (mean - base_mean) / base_mean,
                             marker=marker, color=colors[case["dt"]])
    axes[0].axhline(1, color="black", linestyle="--", linewidth=1)
    axes[0].set_xlabel(r"$dt$")
    axes[0].set_ylabel(r"$a_\infty/\Delta\beta$")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.22)
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel(r"$dt$")
    axes[1].set_ylabel(r"relative shift in $\langle Q\rangle$")
    axes[1].grid(alpha=0.22)
    fig.tight_layout()
    fig.savefig(figures_dir / "timestep_summary.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "timestep_summary.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.3, 4.7))
    for case in equilibrium_cases:
        rows = [products[case["case"]]["per_time"][tau]["row"] for tau in TIMES]
        ax.errorbar(
            TIMES, [row["a_fit"] for row in rows],
            yerr=[
                [row["a_fit"] - row["a_ci_low"] for row in rows],
                [row["a_ci_high"] - row["a_fit"] for row in rows],
            ], marker="o", capsize=2.5, color=colors[case["dt"]],
            label=fr"$dt={case['dt']:.2g}$",
        )
        eq = next(row for row in eq_rows if row["case"] == case["case"])
        ax.axhline(float(eq["weighted_mean_slope"]), color=colors[case["dt"]], linewidth=1)
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xscale("log", base=2)
    ax.set_xlabel(r"window $\tau$")
    ax.set_ylabel(r"equilibrium $a(\tau)$")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures_dir / "equilibrium_estimator_floor.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "equilibrium_estimator_floor.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "analysis_source_sha256": sha256(Path(__file__)),
        "protocol_sha256": sha256(root / "PROTOCOL.md"),
        "cases_sha256": sha256(root / "CASES.tsv"),
        "equilibrium_cases_sha256": sha256(root / "EQUILIBRIUM_CASES.tsv"),
        "short_analysis_sha256": sha256(short_path),
        "pooled_analysis_sha256": sha256(pooled_path),
        "mid_reference_sha256": sha256(mid_reference_path),
        "times": TIMES,
        "min_negative": MIN_NEGATIVE,
        "min_symmetric_pairs": MIN_SYMMETRIC_PAIRS,
        "bootstrap": BOOTSTRAP,
        "streams": STREAMS,
        "blocks_per_stream": BLOCKS_PER_STREAM,
        "inputs": input_rows,
    }
    (analysis_dir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

