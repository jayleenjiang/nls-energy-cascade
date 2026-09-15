#!/usr/bin/env python3
"""Compare timestep and selection-interval controls at fixed n, N, and t."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis_dir", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--horizon", type=float, default=60.0)
    parser.add_argument("--baseline-dt", type=float, default=0.0005)
    parser.add_argument("--baseline-selection", type=float, default=2.0)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)


def convergence(delta: float, se_a: float, se_b: float) -> bool:
    combined = math.hypot(se_a, se_b)
    return (
        math.isfinite(delta) and math.isfinite(combined)
        and abs(delta) <= 0.01 and abs(delta) <= 2.0 * combined
    )


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    aggregate = [
        row for row in read_rows(args.analysis_dir / "controlled_aggregate.csv")
        if close(float(row["horizon"]), args.horizon)
    ]
    pairs = [
        row for row in read_rows(args.analysis_dir / "paired_gc_residuals.csv")
        if close(float(row["horizon"]), args.horizon)
    ]
    baselines = [
        row for row in aggregate
        if close(float(row["dt"]), args.baseline_dt)
        and close(float(row["selection_time"]), args.baseline_selection)
    ]
    if len(baselines) != 2:
        raise ValueError("expected exactly two baseline tilt members")

    member_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    settings = sorted(
        {
            (float(row["dt"]), float(row["selection_time"]))
            for row in aggregate
            if not (
                close(float(row["dt"]), args.baseline_dt)
                and close(float(row["selection_time"]),
                          args.baseline_selection)
            )
        }
    )
    for control_dt, control_selection in settings:
        label = (
            "timestep_halving" if close(control_selection,
                                         args.baseline_selection)
            else "selection_interval"
        )
        all_members_pass = True
        for baseline in sorted(baselines, key=lambda row: float(row["k"])):
            matches = [
                row for row in aggregate
                if close(float(row["dt"]), control_dt)
                and close(float(row["selection_time"]), control_selection)
                and close(float(row["k"]), float(baseline["k"]))
                and int(row["n"]) == int(baseline["n"])
                and int(row["clone_count"]) == int(baseline["clone_count"])
            ]
            if len(matches) != 1:
                raise ValueError(f"missing control member for {label}")
            control = matches[0]
            delta = float(control["mean_scgf"]) - float(baseline["mean_scgf"])
            support = bool(
                int(baseline["support_gate"]) and int(control["support_gate"])
            )
            gate = support and convergence(
                delta, float(baseline["run_scgf_se"]),
                float(control["run_scgf_se"]),
            )
            all_members_pass = all_members_pass and gate
            member_rows.append(
                {
                    "control": label,
                    "n": baseline["n"],
                    "clone_count": baseline["clone_count"],
                    "horizon": baseline["horizon"],
                    "k": baseline["k"],
                    "baseline_dt": baseline["dt"],
                    "control_dt": control["dt"],
                    "baseline_selection": baseline["selection_time"],
                    "control_selection": control["selection_time"],
                    "scgf_delta_control_minus_baseline": delta,
                    "combined_se": math.hypot(
                        float(baseline["run_scgf_se"]),
                        float(control["run_scgf_se"]),
                    ),
                    "both_support_gate": int(support),
                    "member_convergence_gate": int(gate),
                }
            )

        baseline_pairs = [
            row for row in pairs
            if close(float(row["dt"]), args.baseline_dt)
            and close(float(row["selection_time"]), args.baseline_selection)
        ]
        control_pairs = [
            row for row in pairs
            if close(float(row["dt"]), control_dt)
            and close(float(row["selection_time"]), control_selection)
        ]
        if len(baseline_pairs) != 1 or len(control_pairs) != 1:
            raise ValueError(f"expected one paired residual for {label}")
        baseline_pair = baseline_pairs[0]
        control_pair = control_pairs[0]
        residual_delta = (
            float(control_pair["residual"]) - float(baseline_pair["residual"])
        )
        late_delta = (
            float(control_pair["late_half_residual"])
            - float(baseline_pair["late_half_residual"])
        )
        support = bool(
            int(baseline_pair["paired_support_gate"])
            and int(control_pair["paired_support_gate"])
        )
        residual_gate = (
            support
            and convergence(
                residual_delta, float(baseline_pair["residual_se"]),
                float(control_pair["residual_se"]),
            )
            and convergence(
                late_delta, float(baseline_pair["late_half_residual_se"]),
                float(control_pair["late_half_residual_se"]),
            )
        )
        pair_rows.append(
            {
                "control": label,
                "n": baseline_pair["n"],
                "clone_count": baseline_pair["clone_count"],
                "horizon": baseline_pair["horizon"],
                "k": baseline_pair["k"],
                "residual_delta_control_minus_baseline": residual_delta,
                "combined_residual_se": math.hypot(
                    float(baseline_pair["residual_se"]),
                    float(control_pair["residual_se"]),
                ),
                "late_half_delta_control_minus_baseline": late_delta,
                "combined_late_half_se": math.hypot(
                    float(baseline_pair["late_half_residual_se"]),
                    float(control_pair["late_half_residual_se"]),
                ),
                "both_support_gate": int(support),
                "all_member_gates": int(all_members_pass),
                "paired_control_gate": int(residual_gate and all_members_pass),
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "control_members.csv", member_rows)
    write_csv(args.output_dir / "control_pairs.csv", pair_rows)
    audit = {
        "member_comparisons": len(member_rows),
        "passing_member_comparisons": sum(
            int(row["member_convergence_gate"]) for row in member_rows
        ),
        "paired_controls": len(pair_rows),
        "passing_paired_controls": sum(
            int(row["paired_control_gate"]) for row in pair_rows
        ),
    }
    (args.output_dir / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
