#!/usr/bin/env python3
"""Audit frozen controlled trajectories and materialize stream-level splits."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
RAW = REPO / "experiments" / "spectral_gap_controlled_2026-09-07" / "raw"
OUT_AUDIT = HERE / "INPUT_AUDIT.json"
OUT_SPLITS = HERE / "STREAM_SPLITS.csv"
OUT_HASHES = HERE / "provenance" / "INPUT_SHA256.tsv"

N_STREAMS = 64
N_SNAP = 100_001
N_COL = 12
ROW_BYTES = N_COL * np.dtype("<f4").itemsize
EXPECTED_BYTES = N_SNAP * ROW_BYTES

CASES = {
    "driven_dt1e-3": {"T1": 2.0, "T3": 8.0, "dt": 1.0e-3, "base_seed": 2026090801},
    "driven_dt2p5e-4": {"T1": 2.0, "T3": 8.0, "dt": 2.5e-4, "base_seed": 2026090802},
    "equilibrium_dt1e-3": {"T1": 5.0, "T3": 5.0, "dt": 1.0e-3, "base_seed": 2026090803},
    "equilibrium_dt2p5e-4": {"T1": 5.0, "T3": 5.0, "dt": 2.5e-4, "base_seed": 2026090804},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_format(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def split_for(branch: str, case: str, dt: float) -> dict[int, str]:
    ranked = []
    for stream_id in range(N_STREAMS):
        token = f"section4-v1|{branch}|{case}|{dt:.12g}|{stream_id}".encode()
        ranked.append((hashlib.sha256(token).hexdigest(), stream_id))
    ranked.sort()
    result = {}
    for rank, (_, stream_id) in enumerate(ranked):
        result[stream_id] = "train" if rank < 32 else "validation" if rank < 48 else "test"
    return result


def finite_max(current: float, values: np.ndarray) -> float:
    return max(current, float(np.max(np.abs(values))))


def audit_case(case: str, expected: dict[str, float | int]) -> tuple[dict, list[dict]]:
    root = RAW / case
    format_path = root / "FORMAT.txt"
    metadata_path = root / "metadata.csv"
    if not format_path.exists() or not metadata_path.exists():
        raise FileNotFoundError(f"missing FORMAT/metadata for {case}")

    fmt = parse_format(format_path)
    required = {
        "state_precision": "float64",
        "output_precision": "float32",
        "projection": "none",
        "floor": "none",
        "adaptive_step": "no",
        "shape_per_stream": "100001x12",
        "streams": "64",
    }
    format_failures = {k: {"expected": v, "actual": fmt.get(k)} for k, v in required.items()
                       if fmt.get(k) != v}
    numeric_checks = {
        "T1": math.isclose(float(fmt["T1"]), float(expected["T1"]), abs_tol=0.0),
        "T3": math.isclose(float(fmt["T3"]), float(expected["T3"]), abs_tol=0.0),
        "dt": math.isclose(float(fmt["dt"]), float(expected["dt"]), rel_tol=0.0, abs_tol=1e-15),
        "base_seed": int(fmt["base_seed"]) == int(expected["base_seed"]),
    }

    with metadata_path.open(newline="") as handle:
        metadata = list(csv.DictReader(handle))
    meta_by_id = {int(row["stream_id"]): row for row in metadata}
    metadata_ids_ok = sorted(meta_by_id) == list(range(N_STREAMS))

    max_unit_circle_error = 0.0
    max_cos_difference_error = 0.0
    max_redundancy_error = 0.0
    minimum_action = math.inf
    maximum_action = -math.inf
    all_finite = True
    positive_actions = True
    stream_rows: list[dict] = []

    for stream_id in range(N_STREAMS):
        path = root / f"stream_{stream_id:03d}.f32"
        if not path.exists():
            raise FileNotFoundError(path)
        size = path.stat().st_size
        if size != EXPECTED_BYTES:
            raise RuntimeError(f"wrong byte count for {path}: {size} != {EXPECTED_BYTES}")
        data = np.fromfile(path, dtype="<f4").reshape(N_SNAP, N_COL)
        finite = bool(np.isfinite(data).all())
        all_finite &= finite
        if not finite:
            raise RuntimeError(f"nonfinite input in {path}")

        c1, s1, c3, s3 = (data[:, i].astype(np.float64) for i in range(4))
        i2 = data[:, 5].astype(np.float64)
        isum = data[:, 6].astype(np.float64)
        idiff = data[:, 7].astype(np.float64)
        i1 = 0.5 * (isum + idiff)
        i3 = 0.5 * (isum - idiff)
        minimum_action = min(minimum_action, float(i1.min()), float(i2.min()), float(i3.min()))
        maximum_action = max(maximum_action, float(i1.max()), float(i2.max()), float(i3.max()))
        positive_actions &= bool((i1 > 0).all() and (i2 > 0).all() and (i3 > 0).all())

        max_unit_circle_error = finite_max(max_unit_circle_error, c1*c1 + s1*s1 - 1.0)
        max_unit_circle_error = finite_max(max_unit_circle_error, c3*c3 + s3*s3 - 1.0)
        max_cos_difference_error = finite_max(
            max_cos_difference_error, data[:, 4].astype(np.float64) - (c1*c3 + s1*s3)
        )
        redundancies = (
            data[:, 8].astype(np.float64) - (c1 + c3),
            data[:, 9].astype(np.float64) - (c1 - c3),
            data[:, 10].astype(np.float64) - (s1 + s3),
            data[:, 11].astype(np.float64) - (s1 - s3),
        )
        for residual in redundancies:
            max_redundancy_error = finite_max(max_redundancy_error, residual)

        meta = meta_by_id[stream_id]
        if any(int(meta[key]) != 0 for key in
               ("projection_count", "floor_count", "zero_radius_events",
                "midpoint_failure_count", "nonfinite")):
            raise RuntimeError(f"controlled-integrator diagnostic failure in {case} stream {stream_id}")
        stream_rows.append({
            "case": case,
            "dt": expected["dt"],
            "stream_id": stream_id,
            "base_seed": expected["base_seed"],
            "stream_seed": meta["stream_seed"],
            "relative_path": str(path.relative_to(REPO)),
            "bytes": size,
            "sha256": sha256(path),
        })

    passed = (
        not format_failures and all(numeric_checks.values()) and metadata_ids_ok and
        all_finite and positive_actions and max_unit_circle_error <= 5e-7 and
        max_cos_difference_error <= 5e-7 and max_redundancy_error <= 5e-7
    )
    summary = {
        "case": case,
        "format_failures": format_failures,
        "numeric_checks": numeric_checks,
        "metadata_ids_ok": metadata_ids_ok,
        "streams": N_STREAMS,
        "snapshots_per_stream": N_SNAP,
        "raw_snapshots": N_STREAMS * N_SNAP,
        "all_finite": all_finite,
        "positive_reconstructed_actions": positive_actions,
        "minimum_action": minimum_action,
        "maximum_action": maximum_action,
        "max_unit_circle_error": max_unit_circle_error,
        "max_cos_difference_error": max_cos_difference_error,
        "max_redundancy_error": max_redundancy_error,
        "pass": passed,
    }
    return summary, stream_rows


def main() -> None:
    summaries = []
    rows = []
    for case, expected in CASES.items():
        summary, case_rows = audit_case(case, expected)
        summaries.append(summary)
        rows.extend(case_rows)

    split_maps = {
        (branch, case): split_for(branch, case, float(expected["dt"]))
        for branch in ("density", "modes") for case, expected in CASES.items()
    }

    with OUT_SPLITS.open("w", newline="") as handle:
        fields = ["branch", "case", "dt", "stream_id", "base_seed", "stream_seed",
                  "split", "relative_path", "sha256"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for branch in ("density", "modes"):
            for row in rows:
                out = {k: row[k] for k in fields if k not in ("branch", "split")}
                out["branch"] = branch
                out["split"] = split_maps[(branch, row["case"])][int(row["stream_id"])]
                writer.writerow({key: out[key] for key in fields})

    OUT_HASHES.parent.mkdir(parents=True, exist_ok=True)
    with OUT_HASHES.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sha256", "bytes", "path"])
        for row in rows:
            writer.writerow([row["sha256"], row["bytes"], row["relative_path"]])

    split_counts = {}
    for branch in ("density", "modes"):
        for case in CASES:
            counts = {name: 0 for name in ("train", "validation", "test")}
            for value in split_maps[(branch, case)].values():
                counts[value] += 1
            split_counts[f"{branch}:{case}"] = counts

    result = {
        "protocol_version": "section4-v1",
        "input_root": str(RAW),
        "expected_shape": [N_SNAP, N_COL],
        "case_summaries": summaries,
        "split_counts": split_counts,
        "overall_pass": all(item["pass"] for item in summaries),
    }
    OUT_AUDIT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
