#!/usr/bin/env python3
"""Independent integrity audit for the frozen multi-affinity timestep study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


EXPECTED_HEADER = (
    "stream_id,block_id,q_left,q_right,delta_energy,entropy_medium,"
    "entropy_rate,action_current,energy_balance_error"
)
EXPECTED_ROWS = 4_000_256
BLOCKS_PER_STREAM = 31_252
EXPECTED_STREAMS = 128
SOURCE_SHA = "98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3"
BINARY_SHA = "4c4880d721733897d200f2601690da873f5b226aaa1013df23df2351c1dfe7d1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_new_raw(path: Path, expected_hash: str) -> dict:
    digest = hashlib.sha256()
    rows = 0
    min_value = math.inf
    max_value = -math.inf
    with path.open("rb") as handle:
        first = handle.readline()
        digest.update(first)
        header = first.rstrip(b"\r\n").decode("ascii")
        if header != EXPECTED_HEADER:
            raise RuntimeError(f"header mismatch: {path}")
        for raw_line in handle:
            digest.update(raw_line)
            fields = raw_line.rstrip(b"\r\n").split(b",")
            if len(fields) != 9:
                raise RuntimeError(f"field count mismatch at data row {rows}: {path}")
            stream_id = int(fields[0])
            block_id = int(fields[1])
            if stream_id != rows // BLOCKS_PER_STREAM:
                raise RuntimeError(f"stream id mismatch at data row {rows}: {path}")
            if block_id != rows % BLOCKS_PER_STREAM:
                raise RuntimeError(f"block id mismatch at data row {rows}: {path}")
            for item in fields[2:]:
                value = float(item)
                if not math.isfinite(value):
                    raise RuntimeError(f"non-finite value at data row {rows}: {path}")
                min_value = min(min_value, value)
                max_value = max(max_value, value)
            rows += 1
    actual_hash = digest.hexdigest()
    if rows != EXPECTED_ROWS:
        raise RuntimeError(f"row count {rows} != {EXPECTED_ROWS}: {path}")
    if rows // BLOCKS_PER_STREAM != EXPECTED_STREAMS:
        raise RuntimeError(f"stream count mismatch: {path}")
    if actual_hash != expected_hash:
        raise RuntimeError(f"SHA-256 mismatch: {path}")
    return {
        "path": str(path),
        "rows": rows,
        "streams": EXPECTED_STREAMS,
        "finite": True,
        "ids_contiguous": True,
        "numeric_min": min_value,
        "numeric_max": max_value,
        "sha256": actual_hash,
        "status": "PASS",
    }


def csv_rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_dir.resolve()

    raw_manifest = {}
    with (root / "RAW_DATA_MANIFEST.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw_manifest[row["case"]] = row

    raw_results = []
    for case, row in raw_manifest.items():
        path = root / "raw" / f"{case}_n10_tau5_blocks.csv"
        raw_results.append(audit_new_raw(path, row["sha256"]))

    input_results = []
    for row in csv_rows(root / "analysis" / "input_hashes.csv"):
        path = Path(row["path"])
        actual = sha256(path)
        if actual != row["sha256"]:
            raise RuntimeError(f"analysis input hash mismatch: {path}")
        input_results.append({
            "case": row["case"],
            "path": str(path),
            "rows_declared": int(row["rows"]),
            "sha256": actual,
            "status": "PASS",
        })

    expected_analysis_rows = {
        "per_tau_results.csv": 72,
        "slope_bootstrap.csv": 8_000,
        "extrapolation_summary.csv": 16,
        "extrapolation_bootstrap.csv": 16_000,
        "dt_zero_extrapolation.csv": 4,
        "dt_zero_bootstrap.csv": 4_000,
        "three_affinity_finest_summary.csv": 6,
        "equilibrium_weighted_baseline.csv": 2,
        "equilibrium_baseline_difference.csv": 1,
        "input_hashes.csv": 8,
    }
    analysis_results = []
    for name, expected in expected_analysis_rows.items():
        path = root / "analysis" / name
        rows = csv_rows(path)
        if len(rows) != expected:
            raise RuntimeError(f"{name}: {len(rows)} rows != {expected}")
        analysis_results.append({"file": name, "rows": len(rows), "status": "PASS"})

    per_tau = csv_rows(root / "analysis" / "per_tau_results.csv")
    resolved_per_tau = [
        row for row in per_tau if int(row["full_sample_resolved"]) == 1
    ]
    if min(int(row["bootstrap_accepted"]) for row in resolved_per_tau) < 800:
        raise RuntimeError("a per-duration fit has fewer than 800 accepted bootstraps")
    extrapolation = csv_rows(root / "analysis" / "extrapolation_summary.csv")
    if min(int(row["bootstrap_accepted"]) for row in extrapolation) < 800:
        raise RuntimeError("an extrapolation has fewer than 800 accepted bootstraps")

    source = root.parent / "entropy_ft_window_length_2026-09-10/source/NLS_entropy_ft_1905cf.cpp"
    binary = root.parent / "entropy_ft_window_length_2026-09-10/bin/entropy_ft_window"
    source_hash = sha256(source)
    binary_hash = sha256(binary)
    if source_hash != SOURCE_SHA or binary_hash != BINARY_SHA:
        raise RuntimeError("source or binary hash mismatch")

    report = root / "report/multiaffinity_timestep_report.pdf"
    if not report.read_bytes().startswith(b"%PDF"):
        raise RuntimeError("report is not a PDF")
    expected_figures = [
        root / "figures" / f"{stem}.{suffix}"
        for stem in (
            "timestep_summary",
            "finite_time_weak",
            "finite_time_moderate",
            "equilibrium_estimator_floor",
        )
        for suffix in ("pdf", "png")
    ]
    if any(not path.is_file() or path.stat().st_size == 0 for path in expected_figures):
        raise RuntimeError("one or more expected figures are missing or empty")

    output = {
        "status": "PASS",
        "new_raw_files": raw_results,
        "all_analysis_inputs": input_results,
        "analysis_row_counts": analysis_results,
        "bootstrap_minimums": {
            "resolved_per_tau_min_accepted": min(
                int(row["bootstrap_accepted"]) for row in resolved_per_tau
            ),
            "unresolved_per_tau_count": sum(
                int(row["full_sample_resolved"]) == 0 for row in per_tau
            ),
            "extrapolation_min_accepted": min(
                int(row["bootstrap_accepted"]) for row in extrapolation
            ),
        },
        "source_sha256": source_hash,
        "binary_sha256": binary_hash,
        "report_pdf": str(report),
        "report_sha256": sha256(report),
        "figure_count": len(expected_figures),
    }
    (root / "VALIDATION_AUDIT.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
