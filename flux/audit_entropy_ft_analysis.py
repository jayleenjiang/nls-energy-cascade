#!/usr/bin/env python3
"""Independently recompute and audit entropy/action FT analysis products.

This checker reads the raw block CSV files rather than importing aggregation
or fitting routines from ``analyze_entropy_ft.py``.  It is intended as a
second implementation of the numerical bookkeeping, not merely a schema
check of files written by the primary analyzer.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import norm


EXPECTED_HEADER = [
    "stream_id",
    "block_id",
    "q_left",
    "q_right",
    "delta_energy",
    "entropy_medium",
    "entropy_rate",
    "action_current",
    "energy_balance_error",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--supplement-dir", required=True, type=Path)
    parser.add_argument("--adaptive-dir", type=Path)
    parser.add_argument("--time-scaling-dir", type=Path)
    parser.add_argument("--expected-n", default="10,20,30,40")
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def one_row(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise ValueError(f"{path}: expected one row, found {len(rows)}")
    return rows[0]


def grouped(rows: list[dict[str, str]]) -> dict[tuple[int, float], list[dict[str, str]]]:
    result: dict[tuple[int, float], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        result[(int(row["n"]), float(row["tau"]))].append(row)
    return result


def indexed(rows: list[dict[str, str]]) -> dict[tuple[int, float], dict[str, str]]:
    result: dict[tuple[int, float], dict[str, str]] = {}
    for row in rows:
        key = (int(row["n"]), float(row["tau"]))
        if key in result:
            raise ValueError(f"duplicate row for n={key[0]}, tau={key[1]}")
        result[key] = row
    return result


def indexed_observable(
    rows: list[dict[str, str]],
) -> dict[tuple[int, float, str], dict[str, str]]:
    result: dict[tuple[int, float, str], dict[str, str]] = {}
    for row in rows:
        key = (int(row["n"]), float(row["tau"]), row["observable"])
        if key in result:
            raise ValueError(
                f"duplicate row for n={key[0]}, tau={key[1]}, observable={key[2]}"
            )
        result[key] = row
    return result


def grouped_observable(
    rows: list[dict[str, str]],
) -> dict[tuple[int, float, str], list[dict[str, str]]]:
    result: dict[tuple[int, float, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        result[(int(row["n"]), float(row["tau"]), row["observable"])].append(row)
    return result


def close(a: float, b: float, atol: float = 5.0e-11, rtol: float = 5.0e-9) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= atol + rtol * max(abs(a), abs(b))


class Audit:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.checks = 0

    def require(self, condition: bool, message: str) -> None:
        self.checks += 1
        if not condition:
            self.errors.append(message)

    def number(self, actual: float, expected: float, message: str) -> None:
        self.require(close(actual, expected), f"{message}: got {actual:.17g}, expected {expected:.17g}")


def weighted_line_fit(
    x: np.ndarray, y: np.ndarray, variance: np.ndarray
) -> tuple[float, float, float]:
    weights = 1.0 / np.sqrt(np.maximum(variance, np.finfo(float).tiny))
    design = np.column_stack([x, np.ones_like(x)])
    coefficients, *_ = np.linalg.lstsq(
        design * weights[:, None], y * weights, rcond=None
    )
    slope, intercept = coefficients
    fitted = slope * x + intercept
    weighted_mean = np.average(y, weights=weights * weights)
    residual = np.sum((weights * (y - fitted)) ** 2)
    total = np.sum((weights * (y - weighted_mean)) ** 2)
    r_squared = float(1.0 - residual / total) if total > 0.0 else float("nan")
    return float(slope), float(intercept), r_squared


def audit_symmetry(
    audit: Audit,
    label: str,
    values: np.ndarray,
    tau: float,
    rows: list[dict[str, str]],
    summary: dict[str, str],
) -> None:
    if not rows:
        audit.require(int(summary["ft_bins_used"]) == 0, f"{label}: missing bins but fit reports use")
        return
    rows = sorted(rows, key=lambda row: float(row["a_low"]))
    edges = np.asarray([float(rows[0]["a_low"])] + [float(row["a_high"]) for row in rows])
    plus, _ = np.histogram(values.ravel(), bins=edges)
    minus, _ = np.histogram(-values.ravel(), bins=edges)
    sample_count = values.size
    raw = np.full(len(rows), np.nan)
    variance = np.full(len(rows), np.nan)
    used = np.zeros(len(rows), dtype=bool)
    iat = float(summary["autocorrelation_time"])
    for index, row in enumerate(rows):
        audit.require(int(row["plus_count"]) == int(plus[index]), f"{label}: plus bin {index} count")
        audit.require(int(row["minus_count"]) == int(minus[index]), f"{label}: minus bin {index} count")
        plus_four = (plus[index] + 2.0) / (sample_count + 4.0)
        minus_four = (minus[index] + 2.0) / (sample_count + 4.0)
        expected_plot = math.log(plus_four / minus_four) / tau
        audit.number(float(row["symmetry_plus_four"]), expected_plot, f"{label}: plus-four bin {index}")
        if plus[index] > 0 and minus[index] > 0:
            raw[index] = math.log(plus[index] / minus[index]) / tau
            variance[index] = iat * (1.0 / plus[index] + 1.0 / minus[index]) / (tau * tau)
            audit.number(float(row["symmetry_raw"]), raw[index], f"{label}: raw bin {index}")
        else:
            audit.require(math.isnan(float(row["symmetry_raw"])), f"{label}: zero-count raw bin {index} must be NaN")
        used[index] = int(row["fit_used"]) == 1
        if used[index]:
            audit.require(plus[index] > 0 and minus[index] > 0, f"{label}: fitted zero-count bin {index}")
    audit.require(int(np.count_nonzero(used)) == int(summary["ft_bins_used"]), f"{label}: fit-used count")
    if np.count_nonzero(used) >= 3:
        centers = 0.5 * (edges[:-1] + edges[1:])
        slope, intercept, r_squared = weighted_line_fit(
            centers[used], raw[used], variance[used]
        )
        audit.number(float(summary["ft_slope"]), slope, f"{label}: slope")
        audit.number(float(summary["ft_intercept"]), intercept, f"{label}: intercept")
        audit.number(float(summary["ft_r_squared"]), r_squared, f"{label}: R2")


def audit_tail_rows(
    audit: Audit,
    label: str,
    values: np.ndarray,
    rows: list[dict[str, str]],
    columns: dict[str, str],
) -> None:
    flattened = values.ravel()
    ordered = np.sort(flattened)
    count = flattened.size
    mu = float(np.mean(flattened))
    sigma = float(np.std(flattened, ddof=0))
    previous_plus = count + 1
    previous_minus = count + 1
    for index, row in enumerate(sorted(rows, key=lambda item: float(item["A"]))):
        threshold = float(row["A"])
        plus = int(count - np.searchsorted(ordered, threshold, side="left"))
        minus = int(np.searchsorted(ordered, -threshold, side="right"))
        audit.require(int(row[columns["plus_count"]]) == plus, f"{label}: plus tail count at A={threshold}")
        audit.require(int(row[columns["minus_count"]]) == minus, f"{label}: minus tail count at A={threshold}")
        audit.require(plus <= previous_plus and minus <= previous_minus, f"{label}: nonmonotone tail counts at row {index}")
        previous_plus, previous_minus = plus, minus
        plus_raw = plus / count
        minus_raw = minus / count
        audit.number(float(row[columns["plus_raw"]]), plus_raw, f"{label}: plus raw at A={threshold}")
        audit.number(float(row[columns["minus_raw"]]), minus_raw, f"{label}: minus raw at A={threshold}")
        audit.number(float(row[columns["plus_four"]]), (plus + 2.0) / (count + 4.0), f"{label}: plus-four at A={threshold}")
        audit.number(float(row[columns["minus_four"]]), (minus + 2.0) / (count + 4.0), f"{label}: minus-four at A={threshold}")
        if "plus_normal" in columns:
            audit.number(float(row[columns["plus_normal"]]), float(norm.sf(threshold, loc=mu, scale=sigma)), f"{label}: plus normal at A={threshold}")
            audit.number(float(row[columns["minus_normal"]]), float(norm.cdf(-threshold, loc=mu, scale=sigma)), f"{label}: minus normal at A={threshold}")


def audit_adaptive_symmetry(
    audit: Audit,
    label: str,
    values: np.ndarray,
    tau: float,
    rows: list[dict[str, str]],
    summary: dict[str, str],
) -> None:
    flattened = values.ravel()
    positive = flattened[flattened > 0.0]
    negative = -flattened[flattened < 0.0]
    used_reported = int(summary["adaptive_bins_used"])
    if positive.size == 0 or negative.size == 0:
        audit.require(not rows and used_reported == 0, f"{label}: one-sided sample")
        return

    quantile = float(summary["range_quantile"])
    audit.require(0.0 < quantile < 1.0, f"{label}: range quantile")
    overlap = min(float(np.max(positive)), float(np.max(negative)))
    a_max = min(
        float(np.quantile(positive, quantile)),
        float(np.quantile(negative, quantile)),
        overlap,
    )
    iat = float(summary["autocorrelation_time"])
    minimum = float(summary["min_effective_count"])
    raw_target = max(1, int(math.ceil(minimum * iat)))
    minority = min(
        int(np.count_nonzero(positive <= a_max)),
        int(np.count_nonzero(negative <= a_max)),
    )
    bin_count = min(int(summary["max_bins"]), minority // raw_target)
    if bin_count < 3:
        audit.require(
            not rows
            and used_reported == 0
            and math.isnan(float(summary["adaptive_a_max"])),
            f"{label}: insufficient support",
        )
        return

    audit.number(float(summary["adaptive_a_max"]), a_max, f"{label}: a_max")
    audit.require(len(rows) == bin_count, f"{label}: adaptive bin count")
    rows = sorted(rows, key=lambda row: float(row["a_low"]))
    expected_edges = np.linspace(0.0, a_max, bin_count + 1)
    reported_edges = np.asarray(
        [float(rows[0]["a_low"])] + [float(row["a_high"]) for row in rows]
    )
    audit.require(
        np.allclose(reported_edges, expected_edges, atol=5.0e-11, rtol=5.0e-9),
        f"{label}: equal-width range edges",
    )
    plus, _ = np.histogram(flattened, bins=expected_edges)
    minus, _ = np.histogram(-flattened, bins=expected_edges)
    centers = 0.5 * (expected_edges[:-1] + expected_edges[1:])
    raw = np.full(bin_count, np.nan)
    variance = np.full(bin_count, np.nan)
    nonzero = (plus > 0) & (minus > 0)
    raw[nonzero] = np.log(plus[nonzero] / minus[nonzero]) / tau
    variance[nonzero] = iat * (
        1.0 / plus[nonzero] + 1.0 / minus[nonzero]
    ) / (tau * tau)
    expected_used = nonzero & (plus / iat >= minimum) & (minus / iat >= minimum)
    for index, row in enumerate(rows):
        audit.require(int(row["plus_count"]) == int(plus[index]), f"{label}: plus bin {index}")
        audit.require(int(row["minus_count"]) == int(minus[index]), f"{label}: minus bin {index}")
        audit.number(float(row["a_center"]), float(centers[index]), f"{label}: center {index}")
        if nonzero[index]:
            audit.number(float(row["symmetry_raw"]), float(raw[index]), f"{label}: raw bin {index}")
        else:
            audit.require(math.isnan(float(row["symmetry_raw"])), f"{label}: zero-count raw bin {index}")
        audit.require(
            (int(row["fit_used"]) == 1) == bool(expected_used[index]),
            f"{label}: fit mask {index}",
        )
    audit.require(
        int(np.count_nonzero(expected_used)) == used_reported,
        f"{label}: fit-used count",
    )
    if np.count_nonzero(expected_used) >= 3:
        slope, intercept, r_squared = weighted_line_fit(
            centers[expected_used], raw[expected_used], variance[expected_used]
        )
        audit.number(float(summary["adaptive_slope"]), slope, f"{label}: slope")
        audit.number(float(summary["adaptive_intercept"]), intercept, f"{label}: intercept")
        audit.number(float(summary["adaptive_r_squared"]), r_squared, f"{label}: R2")


def stationarity(values: np.ndarray) -> dict[str, float]:
    quarter = max(1, values.shape[1] // 4)
    first = np.mean(values[:, :quarter], axis=1)
    last = np.mean(values[:, -quarter:], axis=1)
    difference = last - first
    se = float(np.std(difference, ddof=1) / math.sqrt(difference.size))
    return {
        "first_quarter_mean": float(np.mean(first)),
        "last_quarter_mean": float(np.mean(last)),
        "last_minus_first": float(np.mean(difference)),
        "paired_se": se,
        "paired_z": float(np.mean(difference) / se) if se > 0.0 else float("nan"),
    }


def audit_tail_time_scaling(
    audit: Audit,
    survival_rows: list[dict[str, str]],
    fit_rows: list[dict[str, str]],
) -> None:
    survival_by_threshold: dict[tuple[int, float], list[dict[str, str]]] = defaultdict(list)
    for row in survival_rows:
        survival_by_threshold[(int(row["n"]), float(row["A"]))].append(row)

    observed: set[tuple[int, float, str]] = set()
    for row in fit_rows:
        key = (int(row["n"]), float(row["A"]), row["tail"])
        audit.require(key not in observed, f"tail-time duplicate {key}")
        observed.add(key)
        tail = row["tail"]
        audit.require(tail in {"plus", "minus"}, f"tail-time {key}: tail label")
        count_column = "plus_count" if tail == "plus" else "minus_count"
        probability_column = "p_plus_raw" if tail == "plus" else "p_minus_raw"
        minimum_count = int(row["minimum_raw_count"])
        maximum_probability = float(row["maximum_probability"])
        qualified = sorted(
            (
                item
                for item in survival_by_threshold[key[:2]]
                if int(item[count_column]) >= minimum_count
                and 0.0 < float(item[probability_column]) <= maximum_probability
            ),
            key=lambda item: float(item["tau"]),
        )
        audit.require(
            len(qualified) == int(row["time_points"]),
            f"tail-time {key}: qualified time count",
        )
        audit.require(len(qualified) >= 3, f"tail-time {key}: minimum fit support")
        time = np.asarray([float(item["tau"]) for item in qualified])
        probability = np.asarray(
            [float(item[probability_column]) for item in qualified]
        )
        slope, intercept = np.polyfit(time, np.log(probability), 1)
        fitted = slope * time + intercept
        residual = float(np.sum((np.log(probability) - fitted) ** 2))
        total = float(np.sum((np.log(probability) - np.mean(np.log(probability))) ** 2))
        r_squared = 1.0 - residual / total if total > 0.0 else float("nan")
        audit.number(float(row["rate_proxy"]), -float(slope), f"tail-time {key}: rate")
        audit.number(float(row["intercept"]), float(intercept), f"tail-time {key}: intercept")
        audit.number(float(row["r_squared"]), r_squared, f"tail-time {key}: R2")
        audit.number(float(row["t_min"]), float(np.min(time)), f"tail-time {key}: t_min")
        audit.number(float(row["t_max"]), float(np.max(time)), f"tail-time {key}: t_max")
        audit.require(
            int(row["minimum_observed_count"])
            == min(int(item[count_column]) for item in qualified),
            f"tail-time {key}: minimum observed count",
        )
        audit.number(
            float(row["maximum_included_probability"]),
            float(np.max(probability)),
            f"tail-time {key}: maximum included probability",
        )

    if fit_rows:
        minimum_count = int(fit_rows[0]["minimum_raw_count"])
        maximum_probability = float(fit_rows[0]["maximum_probability"])
        minimum_points = min(int(row["time_points"]) for row in fit_rows)
        expected: set[tuple[int, float, str]] = set()
        for key, items in survival_by_threshold.items():
            for tail, count_column, probability_column in [
                ("plus", "plus_count", "p_plus_raw"),
                ("minus", "minus_count", "p_minus_raw"),
            ]:
                count = sum(
                    int(item[count_column]) >= minimum_count
                    and 0.0 < float(item[probability_column]) <= maximum_probability
                    for item in items
                )
                if count >= minimum_points:
                    expected.add((key[0], key[1], tail))
        audit.require(observed == expected, "tail-time fit key set")


def main() -> None:
    args = parse_args()
    expected_n = [int(item) for item in args.expected_n.split(",") if item]
    audit = Audit()

    entropy_summary = indexed(read_rows(args.analysis_dir / "ft_summary.csv"))
    action_summary = indexed(read_rows(args.analysis_dir / "action_symmetry_summary.csv"))
    heat_summary = indexed(read_rows(args.analysis_dir / "heat_symmetry_summary.csv"))
    coupling_summary = indexed(read_rows(args.analysis_dir / "coupling_summary.csv"))
    entropy_bins = grouped(read_rows(args.analysis_dir / "ft_symmetric_bins.csv"))
    action_bins = grouped(read_rows(args.analysis_dir / "action_symmetric_bins.csv"))
    heat_bins = grouped(read_rows(args.analysis_dir / "heat_symmetric_bins.csv"))
    core_tails = grouped(read_rows(args.analysis_dir / "action_tail_normal_fit.csv"))
    stationarity_rows = {
        (int(row["n"]), float(row["tau"]), row["observable"]): row
        for row in read_rows(args.analysis_dir / "stationarity_summary.csv")
    }
    supplement_tail_rows = read_rows(args.supplement_dir / "action_two_tail_survival.csv")
    supplement_tails = grouped(supplement_tail_rows)
    supplement_metrics = indexed(read_rows(args.supplement_dir / "action_normal_tail_fit_metrics.csv"))
    adaptive_summary: dict[tuple[int, float, str], dict[str, str]] = {}
    adaptive_bins: dict[tuple[int, float, str], list[dict[str, str]]] = {}
    if args.adaptive_dir is not None:
        adaptive_summary = indexed_observable(
            read_rows(args.adaptive_dir / "adaptive_symmetry_summary.csv")
        )
        adaptive_bins = grouped_observable(
            read_rows(args.adaptive_dir / "adaptive_symmetry_bins.csv")
        )
    time_scaling_rows: list[dict[str, str]] = []
    if args.time_scaling_dir is not None:
        time_scaling_rows = read_rows(
            args.time_scaling_dir / "action_tail_time_scaling.csv"
        )

    expected_keys = set(entropy_summary)
    audit.require(expected_keys == set(action_summary) == set(heat_summary) == set(coupling_summary), "summary key sets differ")

    chain_reports: list[dict[str, object]] = []
    for n in expected_n:
        summary = one_row(args.run_dir / f"n{n}_summary.csv")
        streams = int(summary["n_streams"])
        blocks_per_stream = int(summary["blocks_per_stream"])
        base_tau = float(summary["block_time"])
        blocks_path = args.run_dir / f"n{n}_blocks.csv"
        with blocks_path.open() as stream:
            header = stream.readline().strip().split(",")
        audit.require(header == EXPECTED_HEADER, f"n={n}: raw header")
        raw = np.loadtxt(blocks_path, delimiter=",", skiprows=1)
        if raw.ndim == 1:
            raw = raw.reshape(1, -1)
        audit.require(raw.shape == (streams * blocks_per_stream, len(EXPECTED_HEADER)), f"n={n}: raw shape")
        raw = raw.reshape(streams, blocks_per_stream, len(EXPECTED_HEADER))

        taus = sorted(tau for chain, tau in expected_keys if chain == n)
        chain_reports.append({"n": n, "raw_blocks": streams * blocks_per_stream, "taus": taus})
        for tau in taus:
            key = (n, tau)
            factor_float = tau / base_tau
            factor = int(round(factor_float))
            audit.require(factor > 0 and math.isclose(factor_float, factor, abs_tol=1.0e-12), f"n={n}, t={tau}: aggregation factor")
            usable = (blocks_per_stream // factor) * factor
            grouped_raw = raw[:, :usable, :].reshape(streams, usable // factor, factor, -1)
            q_left = grouped_raw[:, :, :, 2].sum(axis=2)
            q_right = grouped_raw[:, :, :, 3].sum(axis=2)
            entropy = grouped_raw[:, :, :, 5].sum(axis=2) / tau
            action = grouped_raw[:, :, :, 7].mean(axis=2)
            heat = (q_left - q_right) / (2.0 * tau)
            balance = grouped_raw[:, :, :, 8].sum(axis=2) / tau
            values_by_observable = {
                "entropy_rate": entropy,
                "action_current": action,
                "heat_current": heat,
            }

            entropy_row = entropy_summary[key]
            action_row = action_summary[key]
            heat_row = heat_summary[key]
            sample_count = entropy.size
            audit.require(int(entropy_row["n_samples"]) == sample_count, f"n={n}, t={tau}: entropy N")
            audit.number(float(entropy_row["mean_entropy_rate"]), float(np.mean(entropy)), f"n={n}, t={tau}: entropy mean")
            audit.number(float(entropy_row["std_entropy_rate"]), float(np.std(entropy, ddof=1)), f"n={n}, t={tau}: entropy std")
            audit.require(int(entropy_row["negative_count"]) == int(np.count_nonzero(entropy < 0.0)), f"n={n}, t={tau}: entropy negatives")
            audit.number(float(entropy_row["mean_balance_error_rate"]), float(np.mean(balance)), f"n={n}, t={tau}: balance mean")
            audit.number(float(entropy_row["rms_balance_error_rate"]), float(np.sqrt(np.mean(balance * balance))), f"n={n}, t={tau}: balance RMS")

            for label, values, row in [
                ("action", action, action_row),
                ("heat", heat, heat_row),
            ]:
                audit.require(int(row["n_samples"]) == values.size, f"n={n}, t={tau}: {label} N")
                audit.number(float(row["mean"]), float(np.mean(values)), f"n={n}, t={tau}: {label} mean")
                audit.number(float(row["std"]), float(np.std(values, ddof=1)), f"n={n}, t={tau}: {label} std")
                audit.require(int(row["negative_count"]) == int(np.count_nonzero(values < 0.0)), f"n={n}, t={tau}: {label} negatives")

            affinity = 1.0 / float(summary["Tn"]) - 1.0 / float(summary["T1"])
            audit.number(float(heat_row["target_symmetry_slope"]), affinity, f"n={n}, t={tau}: heat affinity")
            audit_symmetry(audit, f"n={n}, t={tau} entropy", entropy, tau, entropy_bins.get(key, []), entropy_row)
            audit_symmetry(audit, f"n={n}, t={tau} action", action, tau, action_bins.get(key, []), action_row)
            audit_symmetry(audit, f"n={n}, t={tau} heat", heat, tau, heat_bins.get(key, []), heat_row)
            if args.adaptive_dir is not None:
                for observable, values in values_by_observable.items():
                    adaptive_key = (n, tau, observable)
                    audit.require(
                        adaptive_key in adaptive_summary,
                        f"n={n}, t={tau}, {observable}: adaptive summary present",
                    )
                    if adaptive_key in adaptive_summary:
                        audit_adaptive_symmetry(
                            audit,
                            f"n={n}, t={tau} adaptive {observable}",
                            values,
                            tau,
                            adaptive_bins.get(adaptive_key, []),
                            adaptive_summary[adaptive_key],
                        )

            coupling = coupling_summary[key]
            x = action.ravel()
            y = heat.ravel()
            design = np.column_stack([x, np.ones_like(x)])
            slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
            residual = y - slope * x - intercept
            audit.number(float(coupling["mean_action_current"]), float(np.mean(x)), f"n={n}, t={tau}: coupling action mean")
            audit.number(float(coupling["mean_heat_current"]), float(np.mean(y)), f"n={n}, t={tau}: coupling heat mean")
            audit.number(float(coupling["pearson_correlation"]), float(np.corrcoef(x, y)[0, 1]), f"n={n}, t={tau}: correlation")
            audit.number(float(coupling["heat_on_action_slope"]), float(slope), f"n={n}, t={tau}: regression slope")
            audit.number(float(coupling["heat_on_action_intercept"]), float(intercept), f"n={n}, t={tau}: regression intercept")
            expected_fraction = float(np.var(residual, ddof=1) / np.var(y, ddof=1))
            audit.number(float(coupling["residual_variance_fraction"]), expected_fraction, f"n={n}, t={tau}: residual fraction")

            for observable, values in values_by_observable.items():
                row = stationarity_rows[(n, tau, observable)]
                expected = stationarity(values)
                for field, value in expected.items():
                    audit.number(float(row[field]), value, f"n={n}, t={tau}: {observable} {field}")

            audit_tail_rows(
                audit,
                f"n={n}, t={tau} core tails",
                action,
                core_tails[key],
                {
                    "plus_count": "plus_count",
                    "minus_count": "minus_count",
                    "plus_raw": "p_x_ge_A_raw",
                    "minus_raw": "p_x_le_minus_A_raw",
                    "plus_four": "p_x_ge_A_plus_four",
                    "minus_four": "p_x_le_minus_A_plus_four",
                    "plus_normal": "p_x_ge_A_normal",
                    "minus_normal": "p_x_le_minus_A_normal",
                },
            )
            audit_tail_rows(
                audit,
                f"n={n}, t={tau} supplement tails",
                action,
                supplement_tails[key],
                {
                    "plus_count": "plus_count",
                    "minus_count": "minus_count",
                    "plus_raw": "p_plus_raw",
                    "minus_raw": "p_minus_raw",
                    "plus_four": "p_plus_plus_four",
                    "minus_four": "p_minus_plus_four",
                    "plus_normal": "p_plus_normal",
                    "minus_normal": "p_minus_normal",
                },
            )
            metric = supplement_metrics[key]
            audit.require(int(metric["sample_count"]) == action.size, f"n={n}, t={tau}: supplement metric N")
            audit.number(float(metric["sample_mu"]), float(np.mean(action)), f"n={n}, t={tau}: supplement metric mean")
            audit.number(float(metric["sample_sigma"]), float(np.std(action, ddof=0)), f"n={n}, t={tau}: supplement metric std")
            if int(metric["joint_tail_fit_success"]) == 1:
                audit.require(float(metric["joint_tail_fit_sigma"]) > 0.0, f"n={n}, t={tau}: fitted sigma positive")
                audit.require(int(metric["joint_tail_fit_plus_points"]) >= 2 and int(metric["joint_tail_fit_minus_points"]) >= 2, f"n={n}, t={tau}: fitted tail points")

    if args.time_scaling_dir is not None:
        audit_tail_time_scaling(audit, supplement_tail_rows, time_scaling_rows)

    result = {
        "status": "PASS" if not audit.errors else "FAIL",
        "checks": audit.checks,
        "error_count": len(audit.errors),
        "errors": audit.errors,
        "chains": chain_reports,
        "run_dir": str(args.run_dir),
        "analysis_dir": str(args.analysis_dir),
        "supplement_dir": str(args.supplement_dir),
        "adaptive_dir": str(args.adaptive_dir) if args.adaptive_dir else None,
        "time_scaling_dir": str(args.time_scaling_dir)
        if args.time_scaling_dir
        else None,
    }
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    args.output_prefix.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    markdown = [
        "# Entropy/action FT analysis audit",
        "",
        f"**Status:** {result['status']}",
        "",
        f"Checks performed: {audit.checks}",
        f"Errors: {len(audit.errors)}",
        "",
    ]
    if audit.errors:
        markdown.extend(["## Errors", ""] + [f"- {error}" for error in audit.errors])
    else:
        markdown.append(
            "All reported aggregates, tail counts, probabilities, symmetry-bin counts, "
            "weighted fits, adaptive-range robustness fits, action-tail time-scaling "
            "fits, stationarity statistics, "
            "and heat--action coupling metrics "
            "match an independent recomputation from the raw block files."
        )
    args.output_prefix.with_suffix(".md").write_text("\n".join(markdown) + "\n")
    print(json.dumps(result, indent=2))
    if audit.errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
