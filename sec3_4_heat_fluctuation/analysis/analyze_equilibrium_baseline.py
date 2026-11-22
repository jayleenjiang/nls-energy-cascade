#!/usr/bin/env python3
"""Frozen finest-timestep equilibrium estimator-baseline analysis."""

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
MIN_NEGATIVE = 500
MIN_SYMMETRIC_PAIRS = 10
MIN_BOOTSTRAP = 800
BOOTSTRAP = 1000
STREAMS = 128
BLOCKS_PER_STREAM = 31252
REFERENCE_DBETA = 0.05714285714285716
SHORT_ANALYSIS_SHA = "100bc49fc4849d3e5f248600a7bb4000af6aadc9161415e1ce2b301de6f9863c"
POOLED_ANALYSIS_SHA = "dc9ed26c662a327459486f245ecf57e245cf4eb407601027de5684afc68d9730"


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


def read_csv(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"empty rows require explicit fields: {path}")
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


def load_cases(path: Path, root: Path) -> list[dict]:
    result = []
    for row in read_csv(path, delimiter="\t"):
        item = dict(row)
        for key in ("T_left", "T_right", "delta_beta", "dt"):
            item[key] = float(item[key])
        for key in ("seed", "bootstrap_seed"):
            item[key] = int(item[key])
        item["path"] = (root / item["blocks_file"]).resolve()
        result.append(item)
    return result


def analyze_case(case: dict, short, pooled) -> dict:
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
        fit = short.matched_fit(values)
        slope_boot, intercept_boot = short.bootstrap_fit(by_stream, fit, multiplicities)
        accepted, boot_se, ci_low, ci_high = percentile(slope_boot)
        n_negative = int(np.count_nonzero(values < 0.0))
        n_positive = int(np.count_nonzero(values > 0.0))
        pairs = int(fit.positive_k.size)
        admissible = bool(
            fit.resolved and n_negative >= MIN_NEGATIVE
            and pairs >= MIN_SYMMETRIC_PAIRS
        )
        row = {
            "case": case["case"],
            "kind": case["kind"],
            "T_left": case["T_left"],
            "T_right": case["T_right"],
            "delta_beta": case["delta_beta"],
            "dt": case["dt"],
            "seed": case["seed"],
            "bootstrap_seed": case["bootstrap_seed"],
            "tau": tau,
            "inverse_tau": 1.0 / tau,
            "N_windows": int(values.size),
            "mean_Q": float(np.mean(values)),
            "std_Q": float(np.std(values, ddof=1)),
            "n_positive": n_positive,
            "n_negative": n_negative,
            "n_zero": int(values.size - n_positive - n_negative),
            "dQ": float(fit.dx),
            "n_symmetric_pairs": pairs,
            "full_sample_resolved": int(fit.resolved),
            "admissible": int(admissible),
            "a_fit": float(fit.slope),
            "a_bootstrap_se": boot_se,
            "a_ci_low": ci_low,
            "a_ci_high": ci_high,
            "ci_contains_zero": int(ci_low <= 0.0 <= ci_high) if math.isfinite(ci_low) else "",
            "intercept": float(fit.intercept),
            "R2": float(fit.r2),
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

    admissible = [tau for tau in TIMES if per_time[tau]["row"]["admissible"]]
    selections = {
        "all_admissible": admissible,
        "tau_ge_10": [tau for tau in admissible if tau >= 10],
    }
    extrapolation_rows = []
    extrapolation_bootstrap_rows = []
    residual_rows = []
    for selection, selected in selections.items():
        summary, residual, bootstrap = pooled.fit_selection(
            case["case"], 0.0, selection, selected, per_time
        )
        summary.update({
            "kind": case["kind"],
            "dt": case["dt"],
            "T_left": case["T_left"],
            "T_right": case["T_right"],
            "seed": case["seed"],
            "bootstrap_seed": case["bootstrap_seed"],
        })
        extrapolation_rows.append(summary)
        for item in residual:
            item.update({"dt": case["dt"]})
        for item in bootstrap:
            item.update({"dt": case["dt"]})
        residual_rows.extend(residual)
        extrapolation_bootstrap_rows.extend(bootstrap)

    slopes = np.asarray([per_time[tau]["row"]["a_fit"] for tau in TIMES])
    boot = np.vstack([per_time[tau]["slope_boot"] for tau in TIMES])
    if not np.all(np.isfinite(slopes)):
        raise RuntimeError(f"{case['case']}: unresolved full-sample ladder slope")
    valid_per_tau = np.sum(np.isfinite(boot), axis=1)
    if np.any(valid_per_tau < MIN_BOOTSTRAP):
        raise RuntimeError(f"{case['case']}: fewer than {MIN_BOOTSTRAP} valid bootstraps")
    ses = np.asarray([np.std(row[np.isfinite(row)], ddof=1) for row in boot])
    weights = 1.0 / np.square(ses)
    combined = float(np.sum(weights * slopes) / np.sum(weights))
    joint_valid = np.all(np.isfinite(boot), axis=0)
    combined_boot = np.full(BOOTSTRAP, np.nan)
    combined_boot[joint_valid] = (
        np.sum(weights[:, None] * boot[:, joint_valid], axis=0) / np.sum(weights)
    )
    accepted, combined_se, low, high = percentile(combined_boot)
    if accepted < MIN_BOOTSTRAP:
        raise RuntimeError(f"{case['case']}: only {accepted} joint ladder replicates")
    n_plus = int(np.count_nonzero(slopes > 0))
    n_minus = int(np.count_nonzero(slopes < 0))
    weighted_row = {
        "case": case["case"],
        "kind": case["kind"],
        "dt": case["dt"],
        "weighted_mean_slope": combined,
        "bootstrap_se": combined_se,
        "ci_low": low,
        "ci_high": high,
        "ci_contains_zero": int(low <= 0.0 <= high),
        "z_from_zero": combined / combined_se,
        "positive_windows": n_plus,
        "negative_windows": n_minus,
        "zero_windows": len(TIMES) - n_plus - n_minus,
        "sign_pattern_tau_order": ";".join(
            "+" if value > 0 else "-" if value < 0 else "0" for value in slopes
        ),
        "fraction_of_dbeta_0p057143": combined / REFERENCE_DBETA,
        "bootstrap_accepted": accepted,
        "bootstrap_total": BOOTSTRAP,
        "weight_rule": "fixed_inverse_bootstrap_variance_all_nine_tau",
    }
    slope_boot_rows = []
    for replicate in range(BOOTSTRAP):
        row = {
            "case": case["case"],
            "dt": case["dt"],
            "replicate": replicate,
            "weighted_mean_slope": combined_boot[replicate],
        }
        for index, tau in enumerate(TIMES):
            row[f"a_tau_{tau}"] = boot[index, replicate]
        slope_boot_rows.append(row)

    return {
        "case": case,
        "per_time": per_time,
        "per_rows": per_rows,
        "symmetric_rows": symmetric_rows,
        "extrapolation_rows": extrapolation_rows,
        "extrapolation_bootstrap_rows": extrapolation_bootstrap_rows,
        "residual_rows": residual_rows,
        "weighted_row": weighted_row,
        "combined_boot": combined_boot,
        "slope_boot_rows": slope_boot_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_dir.resolve()
    analysis = root / "analysis"
    figures = root / "figures"
    analysis.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    baseline_dir = root.parent / "entropy_ft_window_length_2026-09-10"
    pooled_dir = root.parent / "entropy_ft_pooled_extrapolation_2026-09-10"
    short_path = baseline_dir / "analyze_window_dependence.py"
    pooled_path = pooled_dir / "analyze_pooled_extrapolation.py"
    short = load_module("frozen_short_finest_eq", short_path, SHORT_ANALYSIS_SHA)
    pooled = load_module("frozen_pooled_finest_eq", pooled_path, POOLED_ANALYSIS_SHA)

    cases = load_cases(root / "CASES.tsv", root)
    products = []
    input_rows = []
    for case in cases:
        print(f"loading {case['case']}: {case['path']}", flush=True)
        product = analyze_case(case, short, pooled)
        products.append(product)
        input_rows.append({
            "case": case["case"],
            "kind": case["kind"],
            "path": str(case["path"]),
            "rows": STREAMS * BLOCKS_PER_STREAM,
            "sha256": sha256(case["path"]),
        })

    per_rows = [row for product in products for row in product["per_rows"]]
    symmetric_rows = [row for product in products for row in product["symmetric_rows"]]
    extrapolation_rows = [
        row for product in products for row in product["extrapolation_rows"]
    ]
    extrapolation_bootstrap_rows = [
        row for product in products for row in product["extrapolation_bootstrap_rows"]
    ]
    residual_rows = [row for product in products for row in product["residual_rows"]]
    weighted_rows = [product["weighted_row"] for product in products]
    slope_boot_rows = [row for product in products for row in product["slope_boot_rows"]]

    write_csv(analysis / "per_tau_results.csv", per_rows)
    write_csv(analysis / "symmetric_bin_raw_counts.csv", symmetric_rows)
    write_csv(analysis / "extrapolation_summary.csv", extrapolation_rows)
    write_csv(analysis / "extrapolation_bootstrap.csv", extrapolation_bootstrap_rows)
    write_csv(analysis / "fit_residuals.csv", residual_rows)
    write_csv(analysis / "equilibrium_weighted_baseline.csv", weighted_rows)
    write_csv(analysis / "slope_bootstrap.csv", slope_boot_rows)
    write_csv(analysis / "input_hashes.csv", input_rows)

    difference_rows = []
    for left_index in range(len(products)):
        for right_index in range(left_index + 1, len(products)):
            left = products[left_index]
            right = products[right_index]
            diff = right["combined_boot"] - left["combined_boot"]
            accepted, se, low, high = percentile(diff)
            difference_rows.append({
                "left_case": left["case"]["case"],
                "right_case": right["case"]["case"],
                "left_dt": left["case"]["dt"],
                "right_dt": right["case"]["dt"],
                "right_minus_left": (
                    right["weighted_row"]["weighted_mean_slope"]
                    - left["weighted_row"]["weighted_mean_slope"]
                ),
                "bootstrap_se": se,
                "ci_low": low,
                "ci_high": high,
                "ci_contains_zero": int(low <= 0.0 <= high),
                "bootstrap_accepted": accepted,
                "bootstrap_total": BOOTSTRAP,
            })
    write_csv(analysis / "baseline_pairwise_differences.csv", difference_rows)

    externally_quoted = read_csv(root / "EXTERNALLY_QUOTED_BASELINES.csv")
    comparison_rows = []
    for quoted in externally_quoted:
        dt = float(quoted["dt"])
        local = next(row for row in weighted_rows if math.isclose(float(row["dt"]), dt))
        comparison_rows.append({
            "dt": dt,
            "externally_quoted_slope": quoted["weighted_mean_slope"],
            "externally_quoted_ci_low": quoted["ci_low"],
            "externally_quoted_ci_high": quoted["ci_high"],
            "locally_reprocessed_slope": local["weighted_mean_slope"],
            "locally_reprocessed_ci_low": local["ci_low"],
            "locally_reprocessed_ci_high": local["ci_high"],
            "local_minus_quoted": (
                float(local["weighted_mean_slope"])
                - float(quoted["weighted_mean_slope"])
            ),
            "same_point_estimate_within_1e-12": int(
                abs(float(local["weighted_mean_slope"])
                    - float(quoted["weighted_mean_slope"])) <= 1e-12
            ),
        })
    write_csv(analysis / "quoted_vs_local_coarse_baselines.csv", comparison_rows)

    local_values = np.asarray(
        [float(row["weighted_mean_slope"]) for row in weighted_rows], dtype=float
    )
    magnitudes = np.abs(local_values)
    fine = min(weighted_rows, key=lambda row: float(row["dt"]))
    verdict = {
        "fine_case": fine["case"],
        "fine_dt": float(fine["dt"]),
        "fine_weighted_mean_slope": float(fine["weighted_mean_slope"]),
        "fine_ci": [float(fine["ci_low"]), float(fine["ci_high"])],
        "fine_differs_from_zero": not bool(int(fine["ci_contains_zero"])),
        "fine_fraction_of_dbeta_0p057143": float(fine["fraction_of_dbeta_0p057143"]),
        "magnitude_strictly_shrinks_with_dt": bool(np.all(np.diff(magnitudes) < 0)),
        "signs_in_dt_order": [int(np.sign(value)) for value in local_values],
        "quoted_coarse_values_reproduced": bool(
            all(int(row["same_point_estimate_within_1e-12"]) for row in comparison_rows)
        ),
        "interpretation_rule": (
            "compare driven residuals to fine baseline only after reporting the "
            "quoted-versus-local coarse-baseline discrepancy"
        ),
    }
    (analysis / "VERDICT.json").write_text(json.dumps(verdict, indent=2) + "\n")

    colors = {0.0005: "#0072B2", 0.00025: "#E69F00", 0.000125: "#009E73"}
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for product in products:
        rows = product["per_rows"]
        dt = float(product["case"]["dt"])
        ax.errorbar(
            [int(row["tau"]) for row in rows],
            [float(row["a_fit"]) for row in rows],
            yerr=[
                [float(row["a_fit"]) - float(row["a_ci_low"]) for row in rows],
                [float(row["a_ci_high"]) - float(row["a_fit"]) for row in rows],
            ],
            marker="o", capsize=2.5, color=colors[dt], label=fr"$dt={dt:.3g}$",
        )
        ax.axhline(
            float(product["weighted_row"]["weighted_mean_slope"]),
            color=colors[dt], linewidth=1.0, alpha=0.8,
        )
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_xscale("log", base=2)
    ax.set_xlabel(r"window $\tau$")
    ax.set_ylabel(r"equilibrium asymmetry slope $a(\tau)$")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures / "equilibrium_baseline_ladder.pdf", bbox_inches="tight")
    fig.savefig(figures / "equilibrium_baseline_ladder.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.7, 4.5))
    ordered = sorted(weighted_rows, key=lambda row: float(row["dt"]), reverse=True)
    x = np.arange(len(ordered))
    y = np.asarray([float(row["weighted_mean_slope"]) for row in ordered])
    low = np.asarray([float(row["ci_low"]) for row in ordered])
    high = np.asarray([float(row["ci_high"]) for row in ordered])
    ax.errorbar(x, y, yerr=[y - low, high - y], fmt="o", capsize=4, color="#0072B2")
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_xticks(x, [f"{float(row['dt']):.3g}" for row in ordered])
    ax.set_xlabel(r"timestep $dt$")
    ax.set_ylabel(r"weighted equilibrium baseline $a_{\rm eq}$")
    ax.grid(alpha=0.22)
    fig.tight_layout()
    fig.savefig(figures / "equilibrium_baseline_vs_dt.pdf", bbox_inches="tight")
    fig.savefig(figures / "equilibrium_baseline_vs_dt.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "analysis_source_sha256": sha256(Path(__file__)),
        "protocol_sha256": sha256(root / "PROTOCOL.md"),
        "cases_sha256": sha256(root / "CASES.tsv"),
        "externally_quoted_baselines_sha256": sha256(root / "EXTERNALLY_QUOTED_BASELINES.csv"),
        "short_analysis_path": str(short_path),
        "short_analysis_sha256": sha256(short_path),
        "pooled_analysis_path": str(pooled_path),
        "pooled_analysis_sha256": sha256(pooled_path),
        "times": TIMES,
        "min_negative": MIN_NEGATIVE,
        "min_symmetric_pairs": MIN_SYMMETRIC_PAIRS,
        "bootstrap": BOOTSTRAP,
        "streams": STREAMS,
        "blocks_per_stream": BLOCKS_PER_STREAM,
        "inputs": input_rows,
    }
    (analysis / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

