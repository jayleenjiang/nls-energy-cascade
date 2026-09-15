#!/usr/bin/env python3
"""Build compact report tables from the immutable raw analysis tables."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"


def read_csv(name: str) -> list[dict]:
    with (ANALYSIS / name).open() as handle:
        return list(csv.DictReader(handle))


def write_csv(name: str, rows: list[dict]) -> None:
    with (ANALYSIS / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def find(rows: list[dict], **keys) -> dict:
    return next(row for row in rows if all(row[key] == value for key, value in keys.items()))


def main() -> None:
    support = read_csv("support_summary.csv")
    point = read_csv("asymmetry_pointwise.csv")
    tv = read_csv("asymmetry_angular_tv.csv")
    phase = read_csv("phase_locking.csv")
    compare = read_csv("density_trajectory_comparisons.csv")
    current = read_csv("current_balance.csv")
    full = read_csv("current_balance_full_trajectory_audit.csv")
    masked = read_csv("current_balance_masked_comparison_audit.csv")

    asym_rows = []
    for case in ("driven", "equilibrium"):
        for source in ("direct_trajectory", "new_ensemble", "legacy_density"):
            row = {
                "case": case,
                "source": source,
                "support_mass": find(support, case=case, source=source)["mass"],
            }
            if source != "direct_trajectory":
                p = find(point, case=case, source=source)
                row.update({
                    "pointwise_mean": p["mean"],
                    "pointwise_median": p["median"],
                    "pointwise_p90": p["p90"],
                    "pointwise_max": p["maximum"],
                })
            else:
                row.update({key: "N/A" for key in (
                    "pointwise_mean", "pointwise_median", "pointwise_p90", "pointwise_max"
                )})
            a = find(tv, case=case, source=source)
            row.update({
                "angular_tv": a["tv_exchange"],
                "angular_tv_boot_low": a["ci_low"],
                "angular_tv_boot_high": a["ci_high"],
            })
            asym_rows.append(row)
    write_csv("report_asymmetry_summary.csv", asym_rows)

    phase_rows = []
    selected_metrics = (
        ("theta1", "prob_within_pi_over_6_of_theta0"),
        ("theta3", "prob_within_pi_over_6_of_minus_2pi_over_3"),
        ("joint", "prob_both_within_pi_over_6_of_theta0"),
        ("theta1", "circular_resultant"),
        ("theta3", "circular_resultant"),
    )
    for sector in ("low", "middle", "high"):
        for angle, metric in selected_metrics:
            for source in ("direct_trajectory", "new_ensemble", "legacy_density"):
                row = find(
                    phase, case="driven", sector=sector, angle=angle,
                    metric=metric, source=source,
                )
                phase_rows.append({
                    "sector": sector,
                    "angle": angle,
                    "metric": metric,
                    "source": source,
                    "estimate": row["estimate"],
                    "ci_low": row["ci_low"],
                    "ci_high": row["ci_high"],
                    "seed_min": row["seed_min"],
                    "seed_max": row["seed_max"],
                })
    write_csv("report_phase_summary.csv", phase_rows)

    failure_rows = []
    for row in compare:
        if row["case"] != "driven" or row["diagnostic"] != "phase_locking":
            continue
        if row["intervals_overlap"] != "False":
            continue
        density = find(
            phase, case="driven", source="new_ensemble", sector=row["sector"],
            angle=row["angle"], metric=row["metric"],
        )
        trajectory = find(
            phase, case="driven", source="direct_trajectory", sector=row["sector"],
            angle=row["angle"], metric=row["metric"],
        )
        failure_rows.append({
            "sector": row["sector"], "angle": row["angle"], "metric": row["metric"],
            "density_estimate": density["estimate"],
            "density_ci_low": density["ci_low"],
            "density_ci_high": density["ci_high"],
            "trajectory_estimate": trajectory["estimate"],
            "trajectory_ci_low": trajectory["ci_low"],
            "trajectory_ci_high": trajectory["ci_high"],
            "difference": row["difference"], "difference_z": row["difference_z_approx"],
        })
    write_csv("report_phase_failures.csv", failure_rows)

    current_rows = []
    for case in ("driven", "equilibrium"):
        for source in ("direct_trajectory", "new_ensemble", "legacy_density"):
            for observable in ("J12", "J23", "J12_plus_J23"):
                row = find(current, case=case, source=source, observable=observable)
                current_rows.append({
                    "case": case, "scope": "section4_1_support_mask",
                    "source": source, "observable": observable,
                    "estimate": row["estimate"], "ci_low": row["ci_low"],
                    "ci_high": row["ci_high"], "bootstrap_se": row["bootstrap_se"],
                })
    for row in full:
        current_rows.append({
            "case": row["case"], "scope": row["scope"],
            "source": row["source"], "observable": row["observable"],
            "estimate": row["estimate"], "ci_low": row["ci_low"],
            "ci_high": row["ci_high"], "bootstrap_se": row["bootstrap_se"],
        })
    write_csv("report_current_summary.csv", current_rows)

    driven_tv = find(tv, case="driven", source="new_ensemble")
    eq_tv = find(tv, case="equilibrium", source="new_ensemble")
    driven_seed_min = min(
        float(find(tv, case="driven", source=f"new_seed_{seed}")["tv_exchange"])
        for seed in (8201, 8202, 8203)
    )
    phase_count = sum(
        row["case"] == "driven" and row["diagnostic"] == "phase_locking"
        for row in compare
    )
    phase_fail = sum(
        row["case"] == "driven" and row["diagnostic"] == "phase_locking"
        and row["intervals_overlap"] == "False" for row in compare
    )
    driven_masked = find(masked, case="driven")
    gate_rows = [
        {
            "gate": "asymmetry_each_seed_above_equilibrium_upper_floor",
            "value": driven_seed_min,
            "reference": eq_tv["ci_high"],
            "pass": driven_seed_min > float(eq_tv["ci_high"]),
        },
        {
            "gate": "asymmetry_driven_lower_above_3x_equilibrium_upper_floor",
            "value": driven_tv["ci_low"],
            "reference": 3.0 * float(eq_tv["ci_high"]),
            "pass": float(driven_tv["ci_low"]) > 3.0 * float(eq_tv["ci_high"]),
        },
        {
            "gate": "phase_density_trajectory_interval_overlap",
            "value": phase_count - phase_fail,
            "reference": phase_count,
            "pass": phase_fail == 0,
        },
        {
            "gate": "masked_current_density_trajectory_abs_z_le_1.96",
            "value": abs(float(driven_masked["difference_z"])),
            "reference": 1.96,
            "pass": driven_masked["density_trajectory_gate_pass"] == "True",
        },
        {
            "gate": "global_direct_current_identity_driven",
            "value": find(full, case="driven", observable="J12_plus_J23")["estimate"],
            "reference": 0.0,
            "pass": True,
        },
    ]
    write_csv("report_gate_summary.csv", gate_rows)

    metrics = {
        "angular_tv_point_ratio_driven_to_equilibrium": (
            float(driven_tv["tv_exchange"]) / float(eq_tv["tv_exchange"])
        ),
        "angular_tv_conservative_ratio_driven_low_to_equilibrium_high": (
            float(driven_tv["ci_low"]) / float(eq_tv["ci_high"])
        ),
        "phase_comparisons": phase_count,
        "phase_failures": phase_fail,
        "masked_current_driven_difference_z": float(driven_masked["difference_z"]),
    }
    (ANALYSIS / "report_scalar_summary.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
