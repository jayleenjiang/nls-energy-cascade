#!/usr/bin/env python3
"""Integrity, provenance, and finite-value audit for controlled n=3 data."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CASES = (
    "driven_dt1e-3",
    "driven_dt2p5e-4",
    "equilibrium_dt1e-3",
    "equilibrium_dt2p5e-4",
)
NSTREAM = 64
NSNAP = 100001
NOBS = 12
EXPECTED_BYTES = NSNAP * NOBS * 4


def hash_file(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    results = []
    manifest = []
    all_pass = True
    for case in CASES:
        raw = ROOT / "raw" / case
        files = sorted(raw.glob("stream_*.f32"))
        size_ok = len(files) == NSTREAM and all(p.stat().st_size == EXPECTED_BYTES for p in files)
        finite = True
        combo_max = np.zeros(4)
        for path in files:
            x = np.memmap(path, mode="r", dtype="<f4", shape=(NSNAP, NOBS))
            finite &= bool(np.isfinite(x).all())
            combo_max = np.maximum(combo_max, [
                np.max(np.abs(x[:, 8].astype(float) - x[:, 0] - x[:, 2])),
                np.max(np.abs(x[:, 9].astype(float) - x[:, 0] + x[:, 2])),
                np.max(np.abs(x[:, 10].astype(float) - x[:, 1] - x[:, 3])),
                np.max(np.abs(x[:, 11].astype(float) - x[:, 1] + x[:, 3])),
            ])
            manifest.append((hash_file(path), str(path.relative_to(ROOT))))
        with (raw / "metadata.csv").open(newline="") as f:
            metadata = list(csv.DictReader(f))
        ids_ok = [int(r["stream_id"]) for r in metadata] == list(range(NSTREAM))
        counts = {
            key: sum(int(r[key]) for r in metadata)
            for key in (
                "projection_count", "floor_count", "zero_radius_events",
                "midpoint_failure_count", "nonfinite",
            )
        }
        min_action = min(float(r["min_sampled_action"]) for r in metadata)
        case_pass = (
            size_ok and finite and ids_ok and len(metadata) == NSTREAM
            and all(v == 0 for v in counts.values())
        )
        all_pass &= case_pass
        results.append({
            "case": case,
            "stream_count": len(files),
            "expected_stream_count": NSTREAM,
            "bytes_per_stream": EXPECTED_BYTES,
            "size_pass": size_ok,
            "finite_pass": finite,
            "metadata_rows": len(metadata),
            "stream_ids_contiguous": ids_ok,
            **counts,
            "minimum_sampled_action": min_action,
            "max_float32_combination_roundoff": combo_max.tolist(),
            "pass": case_pass,
        })
        for extra in ("metadata.csv", "FORMAT.txt", "run.log"):
            p = raw / extra
            manifest.append((hash_file(p), str(p.relative_to(ROOT))))

    for path in (
        ROOT / "PROTOCOL.md",
        ROOT / "source" / "NLS_stationary_autocorr_controlled.cpp",
        ROOT / "bin" / "NLS_stationary_autocorr_controlled",
        ROOT / "analyze_controlled.py",
        ROOT / "analysis" / "controlled_summary.json",
        ROOT / "analysis" / "early_odd_damped_fits.csv",
        ROOT / "analysis" / "previous_vs_controlled_rates.csv",
        ROOT / "make_figures.py",
        ROOT / "make_report_tables.py",
        ROOT / "audit_results.py",
        ROOT / "figures" / "controlled_local_rates.pdf",
        ROOT / "figures" / "controlled_local_rates.png",
        ROOT / "figures" / "controlled_odd_damped_fits.pdf",
        ROOT / "figures" / "controlled_odd_damped_fits.png",
        ROOT / "report" / "controlled_relaxation_report.tex",
        ROOT / "report" / "controlled_relaxation_report.pdf",
        ROOT / "report" / "generated_tables.tex",
    ):
        manifest.append((hash_file(path), str(path.relative_to(ROOT))))

    (ROOT / "provenance").mkdir(exist_ok=True)
    with (ROOT / "provenance" / "raw_artifact_manifest.sha256").open("w") as f:
        for digest, path in sorted(manifest, key=lambda x: x[1]):
            f.write(f"{digest}  {path}\n")
    with (ROOT / "provenance" / "integrity_audit.json").open("w") as f:
        json.dump({"status": "PASS" if all_pass else "FAIL", "cases": results}, f, indent=2)
    print(json.dumps({"status": "PASS" if all_pass else "FAIL", "cases": results}, indent=2))
    raise SystemExit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
