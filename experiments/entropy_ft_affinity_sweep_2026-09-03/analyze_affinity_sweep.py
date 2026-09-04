#!/usr/bin/env python3
"""Frozen matched-bin/CLT analysis for the n=10 heat affinity sweep."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_T = 20
TIMES = (20, 40, 80, 160, 320, 640)
BINS_PER_SD = 20
MIN_COUNT = 10
MIN_PAIRS = 3
MIN_BOOTSTRAP_RESOLVED = 800
EXPECTED_STREAMS = 128
EXPECTED_BLOCKS_PER_STREAM = 7813
EXPECTED_ROWS = EXPECTED_STREAMS * EXPECTED_BLOCKS_PER_STREAM


@dataclass(frozen=True)
class Case:
    name: str
    t_left: float
    t_right: float
    path: Path

    @property
    def delta_beta(self) -> float:
        return 1.0 / self.t_right - 1.0 / self.t_left


@dataclass
class FitResult:
    resolved: bool
    dx: float
    kmin: int
    kmax: int
    k_lo: int | None
    k_hi: int | None
    positive_k: np.ndarray
    positive_indices: np.ndarray
    negative_indices: np.ndarray
    counts: np.ndarray
    a_fit: float
    intercept: float
    a_se: float
    r2: float


def parse_case(text: str) -> Case:
    parts = text.split(":", 3)
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            "case must be NAME:T_LEFT:T_RIGHT:PATH"
        )
    return Case(parts[0], float(parts[1]), float(parts[2]), Path(parts[3]))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_heat_matrix(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(path)
    data = np.loadtxt(
        path,
        delimiter=",",
        skiprows=1,
        usecols=(0, 1, 2, 3),
        dtype=np.float64,
    )
    if data.shape != (EXPECTED_ROWS, 4):
        raise RuntimeError(
            f"{path}: expected {(EXPECTED_ROWS, 4)}, got {data.shape}"
        )
    if not np.isfinite(data).all():
        raise RuntimeError(f"{path}: non-finite selected fields")
    order = np.lexsort((data[:, 1], data[:, 0]))
    ordered = data[order]
    expected_stream = np.repeat(
        np.arange(EXPECTED_STREAMS), EXPECTED_BLOCKS_PER_STREAM
    )
    expected_block = np.tile(
        np.arange(EXPECTED_BLOCKS_PER_STREAM), EXPECTED_STREAMS
    )
    if not np.array_equal(ordered[:, 0].astype(np.int64), expected_stream):
        raise RuntimeError(f"{path}: stream IDs/order are incomplete")
    if not np.array_equal(ordered[:, 1].astype(np.int64), expected_block):
        raise RuntimeError(f"{path}: block IDs/order are incomplete")
    q = 0.5 * (ordered[:, 2] - ordered[:, 3])
    return q.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM)


def aggregate_windows(q_by_stream: np.ndarray, t: int) -> np.ndarray:
    multiple = t // BASE_T
    if t % BASE_T != 0:
        raise ValueError(f"window {t} is not a multiple of {BASE_T}")
    n_window = q_by_stream.shape[1] // multiple
    trimmed = q_by_stream[:, : n_window * multiple]
    return trimmed.reshape(EXPECTED_STREAMS, n_window, multiple).sum(axis=2)


def contiguous_true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    padded = np.concatenate(([False], mask, [False])).astype(np.int8)
    delta = np.diff(padded)
    starts = np.flatnonzero(delta == 1)
    ends = np.flatnonzero(delta == -1) - 1
    return list(zip(starts.tolist(), ends.tolist()))


def fit_weighted_ratio(
    counts: np.ndarray,
    positive_indices: np.ndarray,
    negative_indices: np.ndarray,
    x_positive: np.ndarray,
) -> tuple[float, float, float, float] | None:
    c_plus = counts[positive_indices].astype(np.float64)
    c_minus = counts[negative_indices].astype(np.float64)
    if np.any(c_plus < MIN_COUNT) or np.any(c_minus < MIN_COUNT):
        return None
    ratio = np.log(c_plus / c_minus)
    variance = 1.0 / c_plus + 1.0 / c_minus
    weight = 1.0 / variance
    design = np.column_stack((x_positive, np.ones_like(x_positive)))
    normal = design.T @ (weight[:, None] * design)
    rhs = design.T @ (weight * ratio)
    try:
        coefficient = np.linalg.solve(normal, rhs)
        covariance = np.linalg.inv(normal)
    except np.linalg.LinAlgError:
        return None
    fitted = design @ coefficient
    residual = ratio - fitted
    weighted_mean = np.sum(weight * ratio) / np.sum(weight)
    ss_residual = np.sum(weight * residual * residual)
    ss_total = np.sum(weight * (ratio - weighted_mean) ** 2)
    r2 = 1.0 - ss_residual / ss_total if ss_total > 0.0 else math.nan
    return (
        float(coefficient[0]),
        float(coefficient[1]),
        float(math.sqrt(covariance[0, 0])),
        float(r2),
    )


def matched_bin_fit(values: np.ndarray) -> FitResult:
    sample_sd = float(np.std(values, ddof=1))
    dx = sample_sd / BINS_PER_SD
    if not math.isfinite(dx) or dx <= 0.0:
        return FitResult(
            False, dx, 0, 0, None, None,
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            math.nan, math.nan, math.nan, math.nan,
        )
    kmin = int(math.floor(float(np.min(values)) / dx) - 1)
    kmax = int(math.ceil(float(np.max(values)) / dx) + 1)
    integer_centers = np.arange(kmin, kmax + 1, dtype=np.int64)
    edges = np.arange(kmin - 0.5, kmax + 1.5, dtype=np.float64) * dx
    counts, _ = np.histogram(values, bins=edges)
    reliable = counts >= MIN_COUNT
    eligible = []
    for lo, hi in contiguous_true_runs(reliable):
        if integer_centers[lo] < 0 and integer_centers[hi] > 0:
            eligible.append((lo, hi))
    if not eligible:
        return FitResult(
            False, dx, kmin, kmax, None, None,
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64), counts,
            math.nan, math.nan, math.nan, math.nan,
        )
    lo, hi = sorted(eligible, key=lambda pair: (-(pair[1] - pair[0] + 1), pair[0]))[0]
    k_lo = int(integer_centers[lo])
    k_hi = int(integer_centers[hi])
    k_symmetric = min(-k_lo, k_hi)
    positive_k = np.arange(1, k_symmetric + 1, dtype=np.int64)
    positive_indices = positive_k - kmin
    negative_indices = -positive_k - kmin
    valid = (
        (positive_indices >= 0)
        & (positive_indices < counts.size)
        & (negative_indices >= 0)
        & (negative_indices < counts.size)
    )
    positive_k = positive_k[valid]
    positive_indices = positive_indices[valid]
    negative_indices = negative_indices[valid]
    keep = (
        (counts[positive_indices] >= MIN_COUNT)
        & (counts[negative_indices] >= MIN_COUNT)
    )
    positive_k = positive_k[keep]
    positive_indices = positive_indices[keep]
    negative_indices = negative_indices[keep]
    if positive_k.size < MIN_PAIRS:
        return FitResult(
            False, dx, kmin, kmax, k_lo, k_hi,
            positive_k, positive_indices, negative_indices, counts,
            math.nan, math.nan, math.nan, math.nan,
        )
    fitted = fit_weighted_ratio(
        counts, positive_indices, negative_indices, positive_k * dx
    )
    if fitted is None:
        return FitResult(
            False, dx, kmin, kmax, k_lo, k_hi,
            positive_k, positive_indices, negative_indices, counts,
            math.nan, math.nan, math.nan, math.nan,
        )
    a_fit, intercept, a_se, r2 = fitted
    return FitResult(
        True, dx, kmin, kmax, k_lo, k_hi,
        positive_k, positive_indices, negative_indices, counts,
        a_fit, intercept, a_se, r2,
    )


def bootstrap_slopes(
    values_by_stream: np.ndarray,
    fit: FitResult,
    multiplicities: np.ndarray,
) -> np.ndarray:
    result = np.full(multiplicities.shape[0], np.nan, dtype=np.float64)
    if not fit.resolved:
        return result
    edges = np.arange(
        fit.kmin - 0.5, fit.kmax + 1.5, dtype=np.float64
    ) * fit.dx
    counts_by_stream = np.vstack(
        [np.histogram(row, bins=edges)[0] for row in values_by_stream]
    )
    bootstrap_counts = multiplicities @ counts_by_stream
    x_positive = fit.positive_k * fit.dx
    for index, counts in enumerate(bootstrap_counts):
        fitted = fit_weighted_ratio(
            counts, fit.positive_indices, fit.negative_indices, x_positive
        )
        if fitted is not None:
            result[index] = fitted[0]
    return result


def finite_or_nan(value: float) -> float:
    return float(value) if math.isfinite(float(value)) else math.nan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026090499)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    symmetric_rows: list[dict] = []
    extrapolation_rows: list[dict] = []
    crossover_rows: list[dict] = []
    hash_rows: list[dict] = []
    case_products: dict[str, dict] = {}

    for case_index, case in enumerate(args.case):
        print(f"loading {case.name}: {case.path}", flush=True)
        q_matrix = load_heat_matrix(case.path)
        hash_rows.append(
            {
                "case": case.name,
                "path": str(case.path.resolve()),
                "rows": EXPECTED_ROWS,
                "sha256": sha256(case.path),
            }
        )
        rng = np.random.default_rng(args.seed + 1009 * case_index)
        multiplicities = rng.multinomial(
            EXPECTED_STREAMS,
            np.full(EXPECTED_STREAMS, 1.0 / EXPECTED_STREAMS),
            size=args.bootstrap,
        )
        per_time: dict[int, dict] = {}
        for t in TIMES:
            values_by_stream = aggregate_windows(q_matrix, t)
            values = values_by_stream.ravel()
            mean = float(np.mean(values))
            sample_std = float(np.std(values, ddof=1))
            standardized = (values - mean) / sample_std
            skewness = float(np.mean(standardized**3))
            excess_kurtosis = float(np.mean(standardized**4) - 3.0)
            n_negative = int(np.count_nonzero(values < 0.0))
            a_gauss = float(2.0 * mean / (sample_std * sample_std))
            if case.delta_beta != 0.0 and mean != 0.0:
                gauss_ft = float(
                    sample_std * sample_std
                    / (2.0 * mean / case.delta_beta)
                )
            else:
                gauss_ft = math.nan
            fit = matched_bin_fit(values)
            boot = bootstrap_slopes(values_by_stream, fit, multiplicities)
            valid_boot = boot[np.isfinite(boot)]
            if valid_boot.size >= MIN_BOOTSTRAP_RESOLVED:
                boot_mean = float(np.mean(valid_boot))
                boot_se = float(np.std(valid_boot, ddof=1))
                boot_lo, boot_hi = np.percentile(valid_boot, [2.5, 97.5])
                boot_lo = float(boot_lo)
                boot_hi = float(boot_hi)
            else:
                boot_mean = boot_se = boot_lo = boot_hi = math.nan

            row = {
                "case": case.name,
                "T_left": case.t_left,
                "T_right": case.t_right,
                "delta_beta": case.delta_beta,
                "t": t,
                "N_windows": values.size,
                "mean_Q": mean,
                "std_Q": sample_std,
                "skew_Q": skewness,
                "excess_kurtosis_Q": excess_kurtosis,
                "n_negative": n_negative,
                "dx": fit.dx,
                "reliable_k_lo": fit.k_lo if fit.k_lo is not None else "",
                "reliable_k_hi": fit.k_hi if fit.k_hi is not None else "",
                "n_symmetric_pairs": fit.positive_k.size,
                "resolved": int(fit.resolved),
                "a_fit": finite_or_nan(fit.a_fit),
                "a_fit_wls_se": finite_or_nan(fit.a_se),
                "a_fit_intercept": finite_or_nan(fit.intercept),
                "a_fit_R2": finite_or_nan(fit.r2),
                "bootstrap_resolved": valid_boot.size,
                "bootstrap_total": args.bootstrap,
                "a_bootstrap_mean": finite_or_nan(boot_mean),
                "a_bootstrap_se": finite_or_nan(boot_se),
                "a_bootstrap_ci_low": finite_or_nan(boot_lo),
                "a_bootstrap_ci_high": finite_or_nan(boot_hi),
                "a_Gauss": a_gauss,
                "gaussFT": finite_or_nan(gauss_ft),
            }
            summary_rows.append(row)
            per_time[t] = {"row": row, "fit": fit, "boot": boot}

            if fit.positive_k.size:
                for k, pos_index, neg_index in zip(
                    fit.positive_k,
                    fit.positive_indices,
                    fit.negative_indices,
                ):
                    c_plus = int(fit.counts[pos_index])
                    c_minus = int(fit.counts[neg_index])
                    symmetric_rows.append(
                        {
                            "case": case.name,
                            "delta_beta": case.delta_beta,
                            "t": t,
                            "k": int(k),
                            "Q_center": float(k * fit.dx),
                            "count_plus": c_plus,
                            "count_minus": c_minus,
                            "log_ratio": float(math.log(c_plus / c_minus)),
                            "log_ratio_se": float(
                                math.sqrt(1.0 / c_plus + 1.0 / c_minus)
                            ),
                        }
                    )

        resolved_times = [
            t for t in TIMES if per_time[t]["row"]["resolved"] == 1
        ]
        a_inf = a_inf_lo = a_inf_hi = ratio = ratio_lo = ratio_hi = math.nan
        extrapolation_r2 = math.nan
        valid_intercepts = np.array([], dtype=np.float64)
        if len(resolved_times) >= 3:
            x = np.array([1.0 / t for t in resolved_times], dtype=np.float64)
            y = np.array(
                [per_time[t]["row"]["a_fit"] for t in resolved_times],
                dtype=np.float64,
            )
            design = np.column_stack((x, np.ones_like(x)))
            coefficient = np.linalg.lstsq(design, y, rcond=None)[0]
            fitted = design @ coefficient
            a_inf = float(coefficient[1])
            ss_residual = float(np.sum((y - fitted) ** 2))
            ss_total = float(np.sum((y - np.mean(y)) ** 2))
            extrapolation_r2 = (
                1.0 - ss_residual / ss_total if ss_total > 0.0 else math.nan
            )
            boot_matrix = np.vstack([per_time[t]["boot"] for t in resolved_times])
            valid_columns = np.all(np.isfinite(boot_matrix), axis=0)
            if int(np.count_nonzero(valid_columns)) >= MIN_BOOTSTRAP_RESOLVED:
                pseudoinverse = np.linalg.pinv(design)
                boot_coefficients = pseudoinverse @ boot_matrix[:, valid_columns]
                valid_intercepts = boot_coefficients[1]
                a_inf_lo, a_inf_hi = np.percentile(valid_intercepts, [2.5, 97.5])
                a_inf_lo = float(a_inf_lo)
                a_inf_hi = float(a_inf_hi)
                if case.delta_beta != 0.0:
                    ratio = a_inf / case.delta_beta
                    ratio_lo = a_inf_lo / case.delta_beta
                    ratio_hi = a_inf_hi / case.delta_beta

        if case.delta_beta == 0.0:
            plateau_change = math.nan
            plateau = "NOT_APPLICABLE"
        elif len(resolved_times) >= 2:
            previous_t, last_t = resolved_times[-2], resolved_times[-1]
            if last_t == 2 * previous_t:
                previous_a = per_time[previous_t]["row"]["a_fit"]
                last_a = per_time[last_t]["row"]["a_fit"]
                plateau_change = abs(last_a - previous_a) / abs(previous_a)
                plateau = "YES" if plateau_change < 0.05 else "NO"
            else:
                plateau_change = math.nan
                plateau = "UNAVAILABLE"
        else:
            plateau_change = math.nan
            plateau = "UNAVAILABLE"

        if len(resolved_times) < 3 or valid_intercepts.size < MIN_BOOTSTRAP_RESOLVED:
            ft_status = "UNRESOLVED"
        elif case.delta_beta == 0.0:
            if a_inf_lo <= 0.0 <= a_inf_hi:
                ft_status = "CONSISTENT_WITH_EQUILIBRIUM"
            else:
                ft_status = "FAIL"
        elif ratio_lo <= 1.0 <= ratio_hi:
            ft_status = "CONSISTENT_WITH_FT"
        else:
            ft_status = "FAIL"

        extrapolation_rows.append(
            {
                "case": case.name,
                "delta_beta": case.delta_beta,
                "resolved_times": ";".join(map(str, resolved_times)),
                "n_resolved_times": len(resolved_times),
                "a_inf": finite_or_nan(a_inf),
                "a_inf_ci_low": finite_or_nan(a_inf_lo),
                "a_inf_ci_high": finite_or_nan(a_inf_hi),
                "a_inf_over_delta_beta": finite_or_nan(ratio),
                "ratio_ci_low": finite_or_nan(ratio_lo),
                "ratio_ci_high": finite_or_nan(ratio_hi),
                "joint_bootstrap_resolved": int(valid_intercepts.size),
                "joint_bootstrap_total": args.bootstrap,
                "extrapolation_R2": finite_or_nan(extrapolation_r2),
                "plateau": plateau,
                "last_doubling_relative_change": finite_or_nan(plateau_change),
                "FT_status": ft_status,
            }
        )

        if resolved_times:
            largest_resolved_t = resolved_times[-1]
            selected_a_gauss = per_time[largest_resolved_t]["row"]["a_Gauss"]
            if case.delta_beta != 0.0:
                a_gauss_ratio = selected_a_gauss / case.delta_beta
            else:
                a_gauss_ratio = math.nan
        else:
            largest_resolved_t = ""
            selected_a_gauss = a_gauss_ratio = math.nan
        crossover_rows.append(
            {
                "case": case.name,
                "delta_beta": case.delta_beta,
                "availability": "AVAILABLE",
                "FT_status": ft_status,
                "a_inf": finite_or_nan(a_inf),
                "a_inf_ci_low": finite_or_nan(a_inf_lo),
                "a_inf_ci_high": finite_or_nan(a_inf_hi),
                "a_inf_over_delta_beta": finite_or_nan(ratio),
                "ratio_ci_low": finite_or_nan(ratio_lo),
                "ratio_ci_high": finite_or_nan(ratio_hi),
                "largest_resolved_t": largest_resolved_t,
                "a_Gauss_at_largest_resolved_t": finite_or_nan(selected_a_gauss),
                "a_Gauss_over_delta_beta": finite_or_nan(a_gauss_ratio),
                "gaussFT_t640": per_time[640]["row"]["gaussFT"],
                "skew_t160": per_time[160]["row"]["skew_Q"],
                "excess_kurtosis_t160": per_time[160]["row"]["excess_kurtosis_Q"],
                "n_negative_t160": per_time[160]["row"]["n_negative"],
            }
        )
        case_products[case.name] = {
            "case": case,
            "per_time": per_time,
            "extrapolation": extrapolation_rows[-1],
            "crossover": crossover_rows[-1],
        }

    if not any(case.delta_beta == 0.0 for case in args.case):
        crossover_rows.append(
            {
                "case": "dbeta_0p000000",
                "delta_beta": 0.0,
                "availability": "UNAVAILABLE_N10_PRODUCTION_NOT_FOUND",
                "FT_status": "UNAVAILABLE",
                "a_inf": math.nan,
                "a_inf_ci_low": math.nan,
                "a_inf_ci_high": math.nan,
                "a_inf_over_delta_beta": math.nan,
                "ratio_ci_low": math.nan,
                "ratio_ci_high": math.nan,
                "largest_resolved_t": "",
                "a_Gauss_at_largest_resolved_t": math.nan,
                "a_Gauss_over_delta_beta": math.nan,
                "gaussFT_t640": math.nan,
                "skew_t160": math.nan,
                "excess_kurtosis_t160": math.nan,
                "n_negative_t160": "",
            }
        )
    crossover_rows.sort(key=lambda row: float(row["delta_beta"]))

    summary_fields = [
        "case", "T_left", "T_right", "delta_beta", "t", "N_windows",
        "mean_Q", "std_Q", "skew_Q", "excess_kurtosis_Q", "n_negative",
        "dx", "reliable_k_lo", "reliable_k_hi", "n_symmetric_pairs",
        "resolved", "a_fit", "a_fit_wls_se", "a_fit_intercept", "a_fit_R2",
        "bootstrap_resolved", "bootstrap_total", "a_bootstrap_mean",
        "a_bootstrap_se", "a_bootstrap_ci_low", "a_bootstrap_ci_high",
        "a_Gauss", "gaussFT",
    ]
    symmetric_fields = [
        "case", "delta_beta", "t", "k", "Q_center", "count_plus",
        "count_minus", "log_ratio", "log_ratio_se",
    ]
    extrapolation_fields = [
        "case", "delta_beta", "resolved_times", "n_resolved_times", "a_inf",
        "a_inf_ci_low", "a_inf_ci_high", "a_inf_over_delta_beta",
        "ratio_ci_low", "ratio_ci_high", "joint_bootstrap_resolved",
        "joint_bootstrap_total", "extrapolation_R2", "plateau",
        "last_doubling_relative_change", "FT_status",
    ]
    crossover_fields = [
        "case", "delta_beta", "availability", "FT_status", "a_inf",
        "a_inf_ci_low", "a_inf_ci_high", "a_inf_over_delta_beta",
        "ratio_ci_low", "ratio_ci_high", "largest_resolved_t",
        "a_Gauss_at_largest_resolved_t", "a_Gauss_over_delta_beta",
        "gaussFT_t640", "skew_t160", "excess_kurtosis_t160",
        "n_negative_t160",
    ]
    write_csv(args.output_dir / "window_summary.csv", summary_rows, summary_fields)
    write_csv(
        args.output_dir / "symmetric_bin_raw_counts.csv",
        symmetric_rows,
        symmetric_fields,
    )
    write_csv(
        args.output_dir / "infinite_time_extrapolation.csv",
        extrapolation_rows,
        extrapolation_fields,
    )
    write_csv(
        args.output_dir / "crossover_summary.csv",
        crossover_rows,
        crossover_fields,
    )
    write_csv(
        args.output_dir / "input_hashes.csv",
        hash_rows,
        ["case", "path", "rows", "sha256"],
    )

    available_cross = [row for row in crossover_rows if row["availability"] == "AVAILABLE"]
    beta = np.array([float(row["delta_beta"]) for row in available_cross])
    order = np.argsort(beta)
    beta = beta[order]
    ordered_rows = [available_cross[i] for i in order]

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    direct_label_used = False
    unresolved_label_used = False
    gaussian_label_used = False
    for row in ordered_rows:
        x_value = float(row["delta_beta"])
        if math.isfinite(float(row["a_inf"])):
            y_value = float(row["a_inf"])
            lo = float(row["a_inf_ci_low"])
            hi = float(row["a_inf_ci_high"])
            if math.isfinite(lo) and math.isfinite(hi):
                ax.errorbar(
                    x_value, y_value,
                    yerr=[[y_value - lo], [hi - y_value]],
                    fmt="o", color="#1f77b4", capsize=3,
                    label="direct-tail $a_\\infty$ with accepted CI"
                    if not direct_label_used else None,
                )
                direct_label_used = True
            else:
                ax.plot(
                    x_value, y_value, "o", markerfacecolor="none",
                    markeredgecolor="#1f77b4",
                    label="full-sample $a_\\infty$; CI gate unresolved"
                    if not unresolved_label_used else None,
                )
                unresolved_label_used = True
        a_gauss = float(row["a_Gauss_at_largest_resolved_t"])
        if math.isfinite(a_gauss):
            ax.plot(
                x_value, a_gauss, "s", color="#d95f02",
                label="$a_{\\mathrm{Gauss}}$ at largest resolved $t$"
                if not gaussian_label_used else None,
            )
            gaussian_label_used = True
    upper = max(0.42, float(np.max(beta)) * 1.05)
    ax.plot([0.0, upper], [0.0, upper], "k--", linewidth=1.2, label="$y=\\Delta\\beta$")
    ax.set_xlim(0.0, upper)
    ax.set_xlabel(r"affinity $\Delta\beta$")
    ax.set_ylabel("slope")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.figures_dir / "slope_crossover_vs_affinity.png", dpi=240)
    fig.savefig(args.figures_dir / "slope_crossover_vs_affinity.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    gauss_values = np.array([float(row["gaussFT_t640"]) for row in ordered_rows])
    ax.plot(beta, gauss_values, "o-", color="#2ca02c")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.2)
    ax.set_xlabel(r"affinity $\Delta\beta$")
    ax.set_ylabel(r"Gaussian FT ratio at $t=640$")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.figures_dir / "gaussFT_vs_affinity.png", dpi=240)
    fig.savefig(args.figures_dir / "gaussFT_vs_affinity.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for row in extrapolation_rows:
        case = case_products[row["case"]]["case"]
        per_time = case_products[row["case"]]["per_time"]
        times = [t for t in TIMES if per_time[t]["row"]["resolved"] == 1]
        if not times:
            continue
        x = np.array([1.0 / t for t in times])
        y = np.array([per_time[t]["row"]["a_fit"] for t in times])
        ax.plot(x, y, "o-", label=rf"$\Delta\beta={case.delta_beta:.3f}$")
    ax.set_xlabel(r"$1/t$")
    ax.set_ylabel(r"direct-tail slope $a_{\mathrm{fit}}$")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.figures_dir / "a_fit_vs_inverse_time.png", dpi=240)
    fig.savefig(args.figures_dir / "a_fit_vs_inverse_time.pdf")
    plt.close(fig)

    print(f"analysis complete: {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
