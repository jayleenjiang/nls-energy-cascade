#!/usr/bin/env python3
"""Audit float64 direct-state replays against the immutable historical output."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
NEW = HERE / "raw" / "direct_state_f64"
OLD = REPO / "experiments" / "spectral_gap_controlled_2026-09-07" / "raw"
OLD_SPLITS = HERE / "STREAM_SPLITS.csv"
NEW_SPLITS = HERE / "DIRECT_STATE_SPLITS.csv"
OUT = HERE / "DIRECT_STATE_AUDIT.json"
HASHES = HERE / "provenance" / "DIRECT_STATE_SHA256.tsv"

N_STREAMS = 64
N_SNAP = 100_001
N_NEW_COL = 5
N_OLD_COL = 12
EXPECTED_BYTES = N_SNAP * N_NEW_COL * 8
TOL = 5.0e-7
CASES = {
    "driven_dt1e-3": (2.0, 8.0, 1.0e-3, 2026090801),
    "driven_dt2p5e-4": (2.0, 8.0, 2.5e-4, 2026090802),
    "equilibrium_dt1e-3": (5.0, 5.0, 1.0e-3, 2026090803),
    "equilibrium_dt2p5e-4": (5.0, 5.0, 2.5e-4, 2026090804),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_format(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def expected_old_observables(state: np.ndarray) -> np.ndarray:
    i1, i2, i3, theta1, theta3 = state.T
    c1, s1 = np.cos(theta1), np.sin(theta1)
    c3, s3 = np.cos(theta3), np.sin(theta3)
    return np.column_stack((
        c1, s1, c3, s3, np.cos(theta1-theta3), i2, i1+i3, i1-i3,
        c1+c3, c1-c3, s1+s3, s1-s3,
    )).astype("<f4")


def main() -> None:
    old_split_rows = list(csv.DictReader(OLD_SPLITS.open(newline="")))
    split_by_key = {
        (r["branch"], r["case"], int(r["stream_id"])): r["split"]
        for r in old_split_rows
    }
    summaries = []
    manifest = []
    split_rows = []
    overall = True
    for case, (t1, t3, dt, base_seed) in CASES.items():
        root = NEW / case
        fmt = parse_format(root / "FORMAT.txt")
        metadata = list(csv.DictReader((root / "metadata.csv").open(newline="")))
        metadata_by_id = {int(row["stream_id"]): row for row in metadata}
        format_ok = (
            fmt.get("output_precision") == "float64"
            and fmt.get("columns") == "I1,I2,I3,theta1,theta3"
            and fmt.get("shape_per_stream") == "100001x5"
            and fmt.get("projection") == "none"
            and fmt.get("floor") == "none"
            and fmt.get("adaptive_step") == "no"
            and math.isclose(float(fmt.get("T1", "nan")), t1, rel_tol=0, abs_tol=0)
            and math.isclose(float(fmt.get("T3", "nan")), t3, rel_tol=0, abs_tol=0)
            and math.isclose(float(fmt.get("dt", "nan")), dt, rel_tol=0, abs_tol=1e-15)
            and int(fmt.get("base_seed", "-1")) == base_seed
        )
        all_finite = True
        all_positive = True
        diagnostics_zero = True
        minimum_action = math.inf
        max_abs_difference = np.zeros(N_OLD_COL, dtype=np.float64)
        for stream_id in range(N_STREAMS):
            path = root / f"stream_{stream_id:03d}.f64"
            if path.stat().st_size != EXPECTED_BYTES:
                raise RuntimeError(f"wrong byte count: {path}")
            state = np.fromfile(path, dtype="<f8").reshape(N_SNAP, N_NEW_COL)
            old = np.fromfile(
                OLD / case / f"stream_{stream_id:03d}.f32", dtype="<f4"
            ).reshape(N_SNAP, N_OLD_COL)
            all_finite &= bool(np.isfinite(state).all())
            all_positive &= bool((state[:, :3] > 0).all())
            minimum_action = min(minimum_action, float(state[:, :3].min()))
            expected = expected_old_observables(state)
            max_abs_difference = np.maximum(
                max_abs_difference,
                np.max(np.abs(old.astype(np.float64)-expected.astype(np.float64)), axis=0),
            )
            meta = metadata_by_id.get(stream_id, {})
            diagnostics_zero &= bool(meta) and all(
                int(meta[key]) == 0 for key in (
                    "projection_count", "floor_count", "zero_radius_events",
                    "midpoint_failure_count", "nonfinite",
                )
            )
            digest = sha256(path)
            rel = str(path.relative_to(REPO))
            manifest.append({"sha256": digest, "bytes": path.stat().st_size, "path": rel})
            for branch in ("density", "modes"):
                split_rows.append({
                    "branch": branch, "case": case, "dt": dt,
                    "stream_id": stream_id, "base_seed": base_seed,
                    "stream_seed": meta["stream_seed"],
                    "split": split_by_key[(branch, case, stream_id)],
                    "relative_path": rel, "sha256": digest,
                })
        cancellation_rows_positive = True
        cancellation_values = []
        if case == "equilibrium_dt2p5e-4":
            for stream_id, row_index, action_index, name in (
                (2, 40629, 2, "I3"), (34, 46751, 0, "I1")
            ):
                path = root / f"stream_{stream_id:03d}.f64"
                state = np.fromfile(path, dtype="<f8").reshape(N_SNAP, N_NEW_COL)
                value = float(state[row_index, action_index])
                cancellation_rows_positive &= value > 0
                cancellation_values.append({
                    "stream_id": stream_id, "row_index": row_index,
                    "action": name, "direct_value": value,
                })
        case_pass = bool(
            format_ok and len(metadata_by_id) == N_STREAMS and all_finite
            and all_positive and diagnostics_zero
            and float(max_abs_difference.max()) <= TOL
            and cancellation_rows_positive
        )
        overall &= case_pass
        summaries.append({
            "case": case, "pass": case_pass, "format_ok": format_ok,
            "all_finite": all_finite, "all_actions_strictly_positive": all_positive,
            "diagnostics_zero": diagnostics_zero,
            "minimum_action": minimum_action,
            "max_abs_difference_by_old_column": max_abs_difference.tolist(),
            "max_abs_difference": float(max_abs_difference.max()),
            "tolerance": TOL,
            "cancellation_rows_positive": cancellation_rows_positive,
            "cancellation_values": cancellation_values,
        })

    with HASHES.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("sha256", "bytes", "path"), delimiter="\t")
        writer.writeheader(); writer.writerows(manifest)
    with NEW_SPLITS.open("w", newline="") as handle:
        fields = ("branch", "case", "dt", "stream_id", "base_seed", "stream_seed",
                  "split", "relative_path", "sha256")
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(split_rows)
    result = {
        "protocol_version": "section4-v1-amendment-001",
        "old_input_audit_preserved": str(HERE / "INPUT_AUDIT.json"),
        "cases": summaries, "overall_pass": overall,
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

