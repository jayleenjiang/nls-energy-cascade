#!/usr/bin/env python3
"""Compile endpoint timestep/selection controls into report-ready outputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment_root", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []
    for n in (10, 40):
        base = (
            args.experiment_root
            / f"analysis_numerical_controls_n{n}"
            / "comparisons"
        )
        pairs = read(base / "control_pairs.csv")
        members = read(base / "control_members.csv")
        for pair in sorted(pairs, key=lambda row: row["control"]):
            control = pair["control"]
            matching = [row for row in members if row["control"] == control]
            if len(matching) != 2:
                raise ValueError(f"expected two members for n={n}, {control}")
            rows.append(
                {
                    "n": n,
                    "control": control,
                    "clone_count": int(pair["clone_count"]),
                    "horizon": float(pair["horizon"]),
                    "k": float(pair["k"]),
                    "one_minus_k": 1.0 - float(pair["k"]),
                    "maximum_absolute_member_delta": max(
                        abs(float(row["scgf_delta_control_minus_baseline"]))
                        for row in matching
                    ),
                    "residual_delta": float(
                        pair["residual_delta_control_minus_baseline"]
                    ),
                    "late_half_residual_delta": float(
                        pair["late_half_delta_control_minus_baseline"]
                    ),
                    "support_gate": int(pair["both_support_gate"]),
                    "member_gates": int(pair["all_member_gates"]),
                    "paired_control_gate": int(pair["paired_control_gate"]),
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "final_control_summary.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    audit = {
        "control_rows": len(rows),
        "passing_control_rows": sum(int(row["paired_control_gate"]) for row in rows),
        "all_endpoint_controls_pass": all(
            int(row["paired_control_gate"]) for row in rows
        ),
    }
    (args.output_dir / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    tex_rows = []
    labels = {
        "selection_interval": "selection interval",
        "timestep_halving": "timestep halving",
    }
    for row in rows:
        tex_rows.append(
            f"{int(row['n'])} & {labels.get(str(row['control']), row['control'])} & "
            f"$({float(row['k']):.1f},{float(row['one_minus_k']):.1f})$ & "
            f"{float(row['maximum_absolute_member_delta']):.5f} & "
            f"{float(row['residual_delta']):+.5f} & "
            f"{float(row['late_half_residual_delta']):+.5f} & "
            f"{'pass' if int(row['paired_control_gate']) else 'fail'} \\\\"  # noqa: W605
        )
    (args.output_dir / "final_control_rows.tex").write_text(
        "\n".join(tex_rows) + "\n\\bottomrule\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
