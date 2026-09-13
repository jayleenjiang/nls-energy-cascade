#!/usr/bin/env python3
"""Mechanical integrity audit for the completed open-chain profile outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPECTED_HEADERS = {
    "profile": "j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta,sample_count",
    "quarter_profiles": "quarter,interval_start,interval_end,j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta,sample_count",
    "quarter_differences": "j,mean_I_q3,mean_I_q4,mean_delta_I_q4_minus_q3,se_delta_I,relative_delta_I,z_delta_I,mean_sin_q3,mean_sin_q4,mean_delta_sin_q4_minus_q3,se_delta_sin,z_delta_sin,sample_count",
    "burnin_checkpoints": "checkpoint_time,j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta",
    "checkpoints": "checkpoint_time,j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta",
    "state_checkpoints": "phase,checkpoint,absolute_time,j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta,sample_count",
    "mass_checkpoints": "series,checkpoint,absolute_time,mean_total_action,se_total_action,mean_I_left,se_I_left,mean_I_right,se_I_right,sample_count",
    "trajectory_diagnostics": "trajectory_id,valid,time_averaged_flux,max_I_left,max_I_right",
    "summary": "bc_id,bc_label,n,T1,Tn,gamma,dt_max,burnin,measure,seed,threads,requested_trajectories,valid_trajectories,nonfinite_trajectories,discarded_trajectories,integration_step_samples,projection_count,mean_flux,se_mean_flux,elapsed_seconds",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_raw(path: Path) -> tuple[list[str], str, list[dict[str, str]]]:
    lines = path.read_text().splitlines()
    metadata = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if not line.startswith("#")]
    if not data_lines:
        return metadata, "", []
    return metadata, data_lines[0], list(csv.DictReader(data_lines))


def expected_rows(suffix: str, n: int) -> int:
    return {
        "profile": n,
        "quarter_profiles": 4 * n,
        "quarter_differences": n,
        "burnin_checkpoints": 4 * n,
        "checkpoints": 4 * n,
        "state_checkpoints": 8 * n,
        "mass_checkpoints": 13,
        "trajectory_diagnostics": 256,
        "summary": 1,
    }[suffix]


def main() -> None:
    errors: list[str] = []
    matrix = list(csv.DictReader((ROOT / "RUN_MATRIX.csv").open()))
    expected_paths: list[Path] = []
    total_valid = total_nonfinite = total_discarded = total_projections = 0

    for run in matrix:
        n = int(run["n"])
        prefix = ROOT / "raw" / (
            f"OPEN_{run['temperature_label']}_n{n}_rep{run['replicate']}"
        )
        for suffix, expected_header in EXPECTED_HEADERS.items():
            path = Path(f"{prefix}_{suffix}.csv")
            expected_paths.append(path)
            if not path.is_file():
                errors.append(f"missing:{path}")
                continue
            metadata, header, rows = read_raw(path)
            if header != expected_header:
                errors.append(f"header:{path}")
            if len(rows) != expected_rows(suffix, n):
                errors.append(
                    f"rows:{path}:{len(rows)}!={expected_rows(suffix, n)}"
                )
            if suffix != "summary":
                joined = "\n".join(metadata)
                for token in (
                    "bc_id=4",
                    "right_boundary_bath=none",
                    "right:b_I=0",
                    "right:b_phi=0",
                    "right:sigma_I=0,sigma_phi=0",
                    f"n={n}",
                    f"seed={run['seed']}",
                ):
                    if token not in joined:
                        errors.append(f"metadata:{token}:{path}")
            if suffix in {
                "profile",
                "quarter_profiles",
                "quarter_differences",
                "state_checkpoints",
                "mass_checkpoints",
            }:
                if any(int(row["sample_count"]) != 256 for row in rows):
                    errors.append(f"sample_count:{path}")
            if suffix == "profile":
                for j, row in enumerate(rows, 1):
                    if not all(
                        math.isfinite(float(row[field]))
                        for field in ("mean_I", "se_mean_I")
                    ):
                        errors.append(f"action_nonfinite:{path}:{j}")
                    sine_fields = ("mean_sin_theta", "se_mean_sin_theta")
                    if j < n and not all(
                        row[field] and math.isfinite(float(row[field]))
                        for field in sine_fields
                    ):
                        errors.append(f"bond_sine_nonfinite:{path}:{j}")
                    if j == n and any(row[field] for field in sine_fields):
                        errors.append(f"terminal_nonbond_not_blank:{path}:{j}")
            if suffix == "trajectory_diagnostics":
                if [int(row["trajectory_id"]) for row in rows] != list(range(256)):
                    errors.append(f"trajectory_ids:{path}")
                if any(int(row["valid"]) != 1 for row in rows):
                    errors.append(f"trajectory_valid:{path}")
                for row in rows:
                    if not all(
                        math.isfinite(float(row[field]))
                        for field in ("time_averaged_flux", "max_I_left", "max_I_right")
                    ):
                        errors.append(f"trajectory_nonfinite:{path}")
                        break
            if suffix == "summary" and rows:
                row = rows[0]
                required = {
                    "bc_id": "4",
                    "bc_label": "OPEN_left_BC1_right_free",
                    "Tn": "0",
                    "threads": "2",
                    "requested_trajectories": "256",
                    "valid_trajectories": "256",
                    "nonfinite_trajectories": "0",
                    "discarded_trajectories": "0",
                }
                for field, value in required.items():
                    if row[field] != value:
                        errors.append(f"summary:{field}:{path}")
                total_valid += int(row["valid_trajectories"])
                total_nonfinite += int(row["nonfinite_trajectories"])
                total_discarded += int(row["discarded_trajectories"])
                total_projections += int(row["projection_count"])

    raw_hash_rows = list(csv.DictReader((ROOT / "RAW_HASHES.csv").open()))
    raw_hashes = {Path(row["path"]).resolve(): row["sha256"] for row in raw_hash_rows}
    if len(raw_hash_rows) != 108:
        errors.append(f"raw_hash_count:{len(raw_hash_rows)}")
    for path in expected_paths:
        if raw_hashes.get(path.resolve()) != sha256(path):
            errors.append(f"raw_hash:{path}")

    analysis_hash_rows = list(
        csv.DictReader((ROOT / "analysis/input_hashes.csv").open())
    )
    if len(analysis_hash_rows) != 108:
        errors.append(f"analysis_hash_count:{len(analysis_hash_rows)}")
    for row in analysis_hash_rows:
        path = Path(row["path"])
        if sha256(path) != row["sha256"]:
            errors.append(f"analysis_hash:{path}")

    curated_count = len(list((ROOT / "curated_results").glob("*.csv")))
    if curated_count != 30:
        errors.append(f"curated_count:{curated_count}")

    output = {
        "status": "PASS" if not errors else "FAIL",
        "logical_conditions": 6,
        "replicate_runs": len(matrix),
        "raw_csv_files": len(expected_paths),
        "raw_hash_rows": len(raw_hash_rows),
        "analysis_input_hash_rows": len(analysis_hash_rows),
        "curated_csv_files": curated_count,
        "total_valid_trajectories": total_valid,
        "total_nonfinite_trajectories": total_nonfinite,
        "total_discarded_trajectories": total_discarded,
        "total_projection_events": total_projections,
        "expected_blank_values": (
            "mean_sin_theta and se_mean_sin_theta at j=n, because no bond n exists"
        ),
        "errors": errors,
    }
    destination = ROOT / "analysis" / "VALIDATION_REPORT.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
