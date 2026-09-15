#!/usr/bin/env python3
"""Audit the frozen Part-1 slow-band resolution-scaling outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "part1_analysis"
CONTROL = HERE.parent / "spectral_gap_controlled_2026-09-07"
PREVIOUS = HERE.parent / "slow_band_multiexponential_2026-09-07" / "run_slow_band.py"

CASES = ("driven_dt2p5e-4", "equilibrium_dt2p5e-4")
STREAM_COUNTS = (64, 128, 256, 512, 1024, 2048)
SEPARATIONS = (0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60)
WINDOWS = ("fixed", "extended")
CONSTRAINTS = ("signed", "nonnegative")
REPLICATES = 500


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cell_key(row: dict[str, str]) -> tuple[str, int, float, str, str]:
    return (
        row["case"], int(row["N_streams"]), float(row["separation"]),
        row["window"], row["constraint"],
    )


def interval(values: list[float]) -> tuple[float, float]:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return math.nan, math.nan
    return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    manifest = json.loads((ANALYSIS / "run_manifest.json").read_text())
    gate = json.loads((HERE / "PART1_GATE.json").read_text())
    summary = read_csv(ANALYSIS / "synthetic_scaling_summary.csv")
    raw = read_csv(ANALYSIS / "synthetic_scaling_raw.csv")
    curves = read_csv(ANALYSIS / "resolution_curve.csv")
    extrap = read_csv(ANALYSIS / "delta020_extrapolation.csv")

    expected_keys = {
        (case, n, delta, window, constraint)
        for case in CASES for n in STREAM_COUNTS for delta in SEPARATIONS
        for window in WINDOWS for constraint in CONSTRAINTS
    }
    summary_keys = {cell_key(row) for row in summary}
    check("manifest complete", manifest.get("status") == "complete",
          f"status={manifest.get('status')}")
    check("summary cell count and grid", len(summary) == 336 and summary_keys == expected_keys,
          f"rows={len(summary)}, unique_cells={len(summary_keys)}, expected=336")
    check("raw row count", len(raw) == 168_000,
          f"rows={len(raw)}, expected=168000")

    by_cell: dict[tuple[str, int, float, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        by_cell[cell_key(row)].append(row)
    replicate_ok = (
        set(by_cell) == expected_keys
        and all({int(row["replicate"]) for row in rows} == set(range(REPLICATES))
                and len(rows) == REPLICATES for rows in by_cell.values())
    )
    check("raw replicate coverage", replicate_ok,
          f"cells={len(by_cell)}, replicates_per_cell={REPLICATES}")

    summary_by_key = {cell_key(row): row for row in summary}
    aggregate_ok = True
    aggregate_max_error = 0.0
    for key, rows in by_cell.items():
        target = summary_by_key[key]
        good = [row for row in rows if row["success"] == "1"]
        orders = [int(row["selected_m"]) for row in good]
        fraction = sum(order >= 2 for order in orders) / len(orders) if orders else 0.0
        matched = [row for row in good if int(row["selected_m"]) >= 2]
        slow_errors = [float(row["abs_error_slow"]) for row in matched]
        fast_errors = [float(row["abs_error_fast"]) for row in matched]
        slow_med = float(np.median(slow_errors)) if slow_errors else math.nan
        fast_med = float(np.median(fast_errors)) if fast_errors else math.nan
        tolerance = max(0.05, key[2] / 4.0)
        passed = fraction >= 0.90 and slow_med <= tolerance and fast_med <= tolerance
        comparisons = [
            (float(target["fraction_select_m_ge_2"]), fraction),
            (float(target["median_abs_error_slow"]), slow_med),
            (float(target["median_abs_error_fast"]), fast_med),
        ]
        for stored, rebuilt in comparisons:
            if math.isfinite(stored) and math.isfinite(rebuilt):
                aggregate_max_error = max(aggregate_max_error, abs(stored - rebuilt))
            elif not (math.isnan(stored) and math.isnan(rebuilt)):
                aggregate_ok = False
        if int(target["successful_fits"]) != len(good):
            aggregate_ok = False
        if int(target["resolution_pass"]) != int(passed):
            aggregate_ok = False
    check("summary reconstructed from raw rows",
          aggregate_ok and aggregate_max_error < 1e-12,
          f"maximum_absolute_reconstruction_difference={aggregate_max_error:.3e}")

    finite_summary = all(
        math.isfinite(float(row[column]))
        for row in summary
        for column in (
            "fraction_select_m_ge_2", "median_abs_error_slow",
            "median_abs_error_fast", "tolerance",
            "fraction_numerical_identifiability_pass",
        )
    )
    check("finite summary statistics", finite_summary,
          f"rows={len(summary)}")

    expected_curve_rows = len(CASES) * len(STREAM_COUNTS) * len(WINDOWS) * len(CONSTRAINTS)
    check("resolution curve row count", len(curves) == expected_curve_rows,
          f"rows={len(curves)}, expected={expected_curve_rows}")
    check("delta=0.20 extrapolation row count", len(extrap) == 8,
          f"rows={len(extrap)}, expected=8")

    delta020 = [row for row in summary if abs(float(row["separation"]) - 0.20) < 1e-12]
    delta020_passes = sum(row["resolution_pass"] == "1" for row in delta020)
    check("frozen Part-2 gate reproduced",
          not gate["part2_allowed"] and delta020_passes == 0
          and manifest["part2_gate"] == gate,
          f"delta020_passing_cells={delta020_passes}; part2_allowed={gate['part2_allowed']}")

    source_hash = sha256(HERE / "run_scaling.py")
    previous_hash = sha256(PREVIOUS)
    input_hashes = {
        case: sha256(CONTROL / "analysis" / f"{case}_trajectory_correlations.npz")
        for case in CASES
    }
    check("analysis source hash", source_hash == manifest["analysis_source_sha256"],
          f"sha256={source_hash}")
    check("previous pipeline hash", previous_hash == manifest["previous_pipeline_sha256"],
          f"sha256={previous_hash}")
    check("input NPZ hashes", input_hashes == manifest["input_npz_sha256"],
          json.dumps(input_hashes, sort_keys=True))

    task_seeds = sorted({
        (row["case"], int(row["N_streams"]), float(row["separation"]), int(row["task_seed"]))
        for row in raw
    })
    seed_consistency = len(task_seeds) == 84
    check("deterministic task-seed coverage", seed_consistency,
          f"unique_case_N_separation_seeds={len(task_seeds)}, expected=84")
    with (ANALYSIS / "task_seeds.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("case", "N_streams", "separation", "task_seed"))
        writer.writerows(task_seeds)

    pass_rows = [row for row in summary if row["resolution_pass"] == "1"]
    with (ANALYSIS / "observed_resolution_passes.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(pass_rows)
    with (ANALYSIS / "delta020_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(delta020)

    output_files = sorted(path for path in ANALYSIS.iterdir() if path.is_file())
    hashes = [
        {"path": str(path.relative_to(HERE)), "sha256": sha256(path), "bytes": path.stat().st_size}
        for path in output_files
    ]
    result = {
        "overall_pass": all(bool(item["passed"]) for item in checks),
        "checks": checks,
        "counts": {
            "cells": len(summary),
            "raw_replicates": len(raw),
            "observed_resolution_pass_cells": len(pass_rows),
            "delta020_pass_cells": delta020_passes,
            "task_seeds": len(task_seeds),
        },
        "part2_gate": gate,
        "artifact_hashes": hashes,
    }
    provenance = HERE / "provenance"
    provenance.mkdir(exist_ok=True)
    (provenance / "integrity_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    with (provenance / "artifact_sha256.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=("path", "sha256", "bytes"))
        writer.writeheader()
        writer.writerows(hashes)
    print(json.dumps(result, indent=2))
    if not result["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
