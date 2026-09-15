#!/usr/bin/env python3
"""Gate the fixed parametric n=2 total-entropy FT analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def row(path: Path) -> dict[str, str]:
    values = rows(path)
    if len(values) != 1:
        raise ValueError(f"expected one row in {path}")
    return values[0]


def numerical_gate(blocks_path: Path) -> tuple[bool, str]:
    summary_path = blocks_path.with_name(
        blocks_path.name.replace("_blocks.csv", "_summary.csv"))
    summary = row(summary_path)
    expected = int(summary["samples"])
    midpoint_failures = int(summary["midpoint_failures"])
    T1 = float(summary["T1"])
    Tn = float(summary["Tn"])
    count = 0
    nonfinite = 0
    formula_max = 0.0
    square_balance = 0.0
    with blocks_path.open(newline="") as stream:
        for value in csv.DictReader(stream):
            numeric = {key: float(value[key]) for key in (
                "energy_start", "energy_end", "q_left", "q_right",
                "delta_energy", "entropy_medium", "energy_balance_error")}
            count += 1
            if not all(math.isfinite(number) for number in numeric.values()):
                nonfinite += 1
                continue
            expected_delta = numeric["energy_end"] - numeric["energy_start"]
            expected_entropy = -numeric["q_left"] / T1 - numeric["q_right"] / Tn
            expected_balance = (
                numeric["q_left"] + numeric["q_right"]
                - numeric["delta_energy"])
            formula_max = max(
                formula_max,
                abs(numeric["delta_energy"] - expected_delta),
                abs(numeric["entropy_medium"] - expected_entropy),
                abs(numeric["energy_balance_error"] - expected_balance),
            )
            square_balance += numeric["energy_balance_error"] ** 2
    rms_balance = math.sqrt(square_balance / count)
    passed = (
        count == expected and nonfinite == 0 and midpoint_failures == 0
        and formula_max <= 5.0e-10 and rms_balance <= 5.0e-5)
    return passed, (
        f"rows={count}/{expected}, nonfinite={nonfinite}, "
        f"midpoint failures={midpoint_failures}, formula max={formula_max:.2e}, "
        f"absolute balance RMS={rms_balance:.2e}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata = json.loads((args.analysis_dir / "analysis_metadata.json").read_text())
    validation = row(args.analysis_dir / "equilibrium_density_validation.csv")
    metrics = {value["dataset"]: value for value in
               rows(args.analysis_dir / "total_entropy_metrics.csv")}
    symmetry = row(args.analysis_dir / "total_entropy_symmetry_summary.csv")
    diagnostics = rows(args.analysis_dir / "density_model_diagnostics.csv")
    stationarity = rows(args.analysis_dir / "endpoint_stationarity.csv")
    checks: list[dict[str, str]] = []

    def add(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name,
                       "status": "PASS" if passed else "BLOCKED",
                       "detail": detail})

    for name, key in (("equilibrium_numerical", "equilibrium_blocks"),
                      ("driven_numerical", "driven_blocks")):
        passed, detail = numerical_gate(Path(metadata[key]))
        add(name, passed, detail)

    support = float(validation["support_fraction"])
    rmse = float(validation["density_log_ratio_rmse"])
    slope = float(validation["density_log_ratio_slope"])
    correlation = float(validation["density_log_ratio_correlation"])
    add(
        "independent_equilibrium_density_validation",
        support >= 0.99 and rmse <= 0.05 and 0.98 <= slope <= 1.02
        and correlation >= 0.995,
        f"support={support:.5f}, RMSE={rmse:.5f}, slope={slope:.5f}, "
        f"correlation={correlation:.5f}",
    )
    optimizer_pass = all(
        int(value["optimizer_success"]) == 1
        and float(value["optimizer_gradient_norm"]) <= 1.0e-5
        for value in diagnostics)
    add(
        "density_optimizer",
        optimizer_pass,
        "; ".join(
            f"{value['dataset']}/{value['fold']}: success="
            f"{value['optimizer_success']}, |grad|="
            f"{float(value['optimizer_gradient_norm']):.2e}"
            for value in diagnostics),
    )

    for dataset, label in (("equilibrium_parametric",
                            "equilibrium_endpoint_error_budget"),
                           ("driven_parametric", "driven_ift")):
        value = metrics[dataset]
        low = float(value["log_ift_stream_bootstrap_low"])
        high = float(value["log_ift_stream_bootstrap_high"])
        ess = float(value["exponential_weight_ess"])
        supported = float(value["support_fraction"])
        estimate = float(value["log_mean_exp_minus_total_entropy"])
        relation_pass = (
            abs(estimate) <= 1.0e-3
            if dataset == "equilibrium_parametric"
            else low <= 0.0 <= high
        )
        add(
            label,
            relation_pass and ess >= 1000 and supported >= 0.99
            and int(value["unique_streams"]) >= 32,
            f"log-IFT={estimate:.6g}, "
            f"CI=[{low:.6g},{high:.6g}], ESS={ess:.1f}, "
            f"support={supported:.5f}",
        )

    exact = metrics["equilibrium_exact_gibbs"]
    exact_value = float(exact["log_mean_exp_minus_total_entropy"])
    add("exact_gibbs_control", abs(exact_value) <= 1.0e-3,
        f"log-IFT={exact_value:.3e}")

    bins = int(symmetry["usable_bins"])
    symmetry_slope = float(symmetry["weighted_slope"])
    symmetry_se = float(symmetry["weighted_slope_se"])
    intercept = float(symmetry["weighted_intercept"])
    add(
        "driven_detailed_ft",
        bins >= 8 and abs(symmetry_slope - 1.0) <= max(0.1, 2.0 * symmetry_se)
        and abs(intercept) <= 0.1,
        f"bins={bins}, slope={symmetry_slope:.6g}+/-{symmetry_se:.3g}, "
        f"intercept={intercept:.6g}",
    )
    sensitivity_path = args.analysis_dir / "model_sensitivity.json"
    if sensitivity_path.exists():
        sensitivity = json.loads(sensitivity_path.read_text())
        sensitivity_pass = len(sensitivity) >= 3 and all(
            float(value["ci_low"]) <= 0.0 <= float(value["ci_high"])
            and abs(float(value["slope"]) - 1.0) <= 0.1
            and int(value["bins"]) >= 8
            for value in sensitivity
        )
        sensitivity_detail = "; ".join(
            f"D{value['D']}K{value['K']}: IFT CI="
            f"[{float(value['ci_low']):.4g},{float(value['ci_high']):.4g}], "
            f"slope={float(value['slope']):.4f}"
            for value in sensitivity
        )
    else:
        sensitivity_pass = False
        sensitivity_detail = "missing model_sensitivity.json"
    add("density_model_sensitivity", sensitivity_pass, sensitivity_detail)
    max_z = max(abs(float(value["z_score"])) for value in stationarity)
    add("endpoint_stationarity", max_z <= 3.0,
        f"maximum |stream-level z|={max_z:.3f}")

    passed = all(value["status"] == "PASS" for value in checks)
    overall = "PASS_SMALL_CHAIN_NESS_TOTAL_ENTROPY_FT" if passed else "BLOCKED"
    payload = {
        "overall": overall,
        "checks": checks,
        "scope": (
            "Finite-time n=2 NESS total-entropy numerical verification; "
            "not a long-chain asymptotic GC proof."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = [
        "# Parametric n=2 total-entropy FT audit", "",
        f"Overall: **{overall}**", "",
        "| Check | Status | Detail |", "|---|---:|---|",
        *(f"| {value['check']} | {value['status']} | {value['detail']} |"
          for value in checks), "",
        "The accepted scope is a finite-time, n=2 NESS total-entropy numerical",
        "verification. It is not a mathematical proof or a long-chain GC result.",
        "",
    ]
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output}: {overall}")


if __name__ == "__main__":
    main()
