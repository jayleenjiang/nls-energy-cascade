#!/usr/bin/env python3
"""Audit the validation-first n=2 total-entropy endpoint pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def one_row(path: Path) -> dict[str, str]:
    values = rows(path)
    if len(values) != 1:
        raise ValueError(f"expected one summary row in {path}")
    return values[0]


def endpoint_numerical_audit(blocks_path: Path) -> tuple[bool, str]:
    summary_path = blocks_path.with_name(
        blocks_path.name.replace("_blocks.csv", "_summary.csv")
    )
    if not blocks_path.exists() or not summary_path.exists():
        return False, f"missing blocks or summary for {blocks_path}"
    summary = one_row(summary_path)
    temperature_left = float(summary["T1"])
    temperature_right = float(summary["Tn"])
    block_time = float(summary["block_time"])
    expected = int(summary["samples"])
    midpoint_failures = int(summary["midpoint_failures"])
    count = 0
    nonfinite = 0
    formula_max = {"delta_energy": 0.0, "entropy": 0.0, "balance": 0.0}
    balance_square_sum = 0.0
    with blocks_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            numeric = {
                key: float(row[key]) for key in (
                    "energy_start", "energy_end", "q_left", "q_right",
                    "delta_energy", "entropy_medium", "energy_balance_error",
                )
            }
            count += 1
            if not all(math.isfinite(value) for value in numeric.values()):
                nonfinite += 1
                continue
            expected_delta = numeric["energy_end"] - numeric["energy_start"]
            expected_entropy = (
                -numeric["q_left"] / temperature_left
                - numeric["q_right"] / temperature_right
            )
            expected_balance = (
                numeric["q_left"] + numeric["q_right"]
                - numeric["delta_energy"]
            )
            formula_max["delta_energy"] = max(
                formula_max["delta_energy"],
                abs(numeric["delta_energy"] - expected_delta),
            )
            formula_max["entropy"] = max(
                formula_max["entropy"],
                abs(numeric["entropy_medium"] - expected_entropy),
            )
            formula_max["balance"] = max(
                formula_max["balance"],
                abs(numeric["energy_balance_error"] - expected_balance),
            )
            balance_square_sum += numeric["energy_balance_error"] ** 2
    rms_balance_rate = (
        math.sqrt(balance_square_sum / count) / block_time if count else math.inf
    )
    reported_rms = float(summary["rms_energy_balance_error_rate"])
    passed = (
        count == expected and nonfinite == 0 and midpoint_failures == 0
        and max(formula_max.values()) <= 5.0e-10
        and abs(rms_balance_rate - reported_rms) <= 5.0e-11
        and rms_balance_rate <= 5.0e-5
    )
    detail = (
        f"rows={count}/{expected}, nonfinite={nonfinite}, "
        f"midpoint failures={midpoint_failures}, formula max="
        f"{max(formula_max.values()):.2e}, balance RMS rate={rms_balance_rate:.2e}"
    )
    return passed, detail


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis_dir", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    validation = rows(args.analysis_dir / "equilibrium_density_validation.csv")
    metrics = {row["dataset"]: row for row in
               rows(args.analysis_dir / "total_entropy_metrics.csv")}
    if not validation or "equilibrium_learned" not in metrics:
        raise SystemExit("incomplete total-entropy analysis")
    heldout = [row for row in validation if row.get("stage") == "heldout_validation"]
    if len(heldout) != 1:
        raise SystemExit("expected one held-out equilibrium validation row")
    selected = heldout[0]
    equilibrium = metrics["equilibrium_learned"]
    exact = metrics["equilibrium_exact_gibbs"]
    driven = metrics["driven_learned"]
    stationarity = rows(args.analysis_dir / "endpoint_stationarity.csv")
    metadata = json.loads(
        (args.analysis_dir / "analysis_metadata.json").read_text()
    )

    checks: list[dict[str, str]] = []

    def add(name: str, passed: bool, detail: str) -> None:
        checks.append({
            "check": name,
            "status": "PASS" if passed else "BLOCKED",
            "detail": detail,
        })

    for name, key in (
        ("equilibrium_numerical_integrity", "equilibrium_blocks"),
        ("driven_numerical_integrity", "driven_blocks"),
    ):
        passed, detail = endpoint_numerical_audit(Path(metadata[key]))
        add(name, passed, detail)

    support = float(selected["support_fraction"])
    slope = float(selected["density_log_ratio_slope"])
    correlation = float(selected["density_log_ratio_correlation"])
    rmse = float(selected["density_log_ratio_rmse"])
    density_pass = (
        support >= 0.95 and 0.9 <= slope <= 1.1
        and correlation >= 0.9 and rmse <= 0.5
    )
    add(
        "equal_temperature_density_ratio",
        density_pass,
        f"support={support:.4f}, slope={slope:.4f}, "
        f"correlation={correlation:.4f}, RMSE={rmse:.4f}",
    )

    ci_low = float(equilibrium["log_ift_stream_bootstrap_low"])
    ci_high = float(equilibrium["log_ift_stream_bootstrap_high"])
    equilibrium_ift_pass = (
        ci_low <= 0.0 <= ci_high
        and float(equilibrium["exponential_weight_ess"]) >= 1000.0
        and int(equilibrium["unique_streams"]) >= 32
    )
    add(
        "learned_equilibrium_integral_ft",
        equilibrium_ift_pass,
        f"log-IFT CI=[{ci_low:.4f},{ci_high:.4f}], "
        f"ESS={float(equilibrium['exponential_weight_ess']):.1f}, "
        f"streams={equilibrium['unique_streams']}",
    )

    exact_log_ift = float(exact["log_mean_exp_minus_total_entropy"])
    add(
        "exact_gibbs_balance_control",
        abs(exact_log_ift) <= 1.0e-3,
        f"log mean exp(-Delta s_total)={exact_log_ift:.3e}",
    )

    driven_support_pass = (
        float(driven["support_fraction"]) >= 0.95
        and float(driven["exponential_weight_ess"]) >= 1000.0
        and int(driven["unique_streams"]) >= 32
    )
    add(
        "driven_endpoint_support",
        driven_support_pass,
        f"support={float(driven['support_fraction']):.4f}, "
        f"ESS={float(driven['exponential_weight_ess']):.1f}, "
        f"streams={driven['unique_streams']}",
    )

    finite_z = [abs(float(row["z_score"])) for row in stationarity]
    max_z = max(finite_z) if finite_z else float("inf")
    add(
        "endpoint_stationarity",
        max_z <= 3.0,
        f"maximum |stream-level z|={max_z:.3f} over {len(finite_z)} checks",
    )

    interpretation_ready = all(
        item["status"] == "PASS" for item in checks
    )
    overall = "PASS_EXPLORATORY_ENDPOINT" if interpretation_ready else "BLOCKED"
    payload = {
        "overall": overall,
        "selected_sigma_bins": float(selected["sigma_bins"]),
        "checks": checks,
        "interpretation": (
            "This is an n=2 endpoint-density control only; it is not a "
            "long-chain detailed-FT result."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = [
        "# n=2 total-entropy endpoint audit",
        "",
        f"Overall: **{overall}**",
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    lines.extend(
        f"| {item['check']} | {item['status']} | {item['detail']} |"
        for item in checks
    )
    lines.extend([
        "",
        "Passing would validate only the exploratory n=2 endpoint estimator;",
        "it would not establish a long-chain detailed fluctuation theorem.",
        "",
    ])
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output}: {overall}")


if __name__ == "__main__":
    main()
