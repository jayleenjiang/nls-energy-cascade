#!/usr/bin/env python3
"""Integrity audit for the frozen n=3 autocorrelation production."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CASES = {
    "driven_T2_T8": 2026090701,
    "equilibrium_T5_T5": 2026090702,
}
N_STREAM = 64
N_SNAP = 100001
N_OBS = 7
EXPECTED_BYTES = N_SNAP * N_OBS * 4


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    result = {"status": "PASS", "cases": {}, "analysis": {}}
    errors = []
    for case, base_seed in CASES.items():
        d = ROOT / "raw" / case
        files = sorted(d.glob("stream_*.f32"))
        rows = list(csv.DictReader((d / "metadata.csv").open()))
        finite = True
        sizes = []
        for p in files:
            sizes.append(p.stat().st_size)
            x = np.memmap(p, dtype="<f4", mode="r")
            if not np.all(np.isfinite(x)):
                finite = False
        seeds = [int(r["seed"]) for r in rows]
        projections = sum(int(r["projection_count"]) for r in rows)
        nonfinite_flags = sum(int(r["nonfinite"]) for r in rows)
        checks = {
            "stream_file_count": len(files),
            "metadata_row_count": len(rows),
            "all_exact_size": bool(sizes and all(x == EXPECTED_BYTES for x in sizes)),
            "expected_bytes_per_stream": EXPECTED_BYTES,
            "all_values_finite": bool(finite),
            "seeds_contiguous": seeds == list(range(base_seed, base_seed + N_STREAM)),
            "nonfinite_flags": nonfinite_flags,
            "projection_count": projections,
            "metadata_sha256": sha256(d / "metadata.csv"),
        }
        result["cases"][case] = checks
        if not (
            len(files) == N_STREAM and len(rows) == N_STREAM and checks["all_exact_size"]
            and finite and checks["seeds_contiguous"] and nonfinite_flags == 0
        ):
            errors.append(case)

        plateau_rows = list(csv.DictReader((ROOT / "analysis" / f"{case}_plateaus.csv").open()))
        result["analysis"][case] = {
            "plateau_rows": len(plateau_rows),
            "all_stationarity_gates_pass": all(r["stationarity_pass"] == "1" for r in plateau_rows),
            "accepted_plateau_observables": [r["observable"] for r in plateau_rows if r["lambda_R"]],
            "autocorrelation_csv_rows": sum(
                1 for _ in (ROOT / "analysis" / f"{case}_autocorrelations.csv").open()
            ) - 1,
        }
        if result["analysis"][case]["plateau_rows"] != 7:
            errors.append(f"{case}: plateau rows")
        if result["analysis"][case]["autocorrelation_csv_rows"] != 623:
            errors.append(f"{case}: autocorrelation rows")
        if not result["analysis"][case]["all_stationarity_gates_pass"]:
            errors.append(f"{case}: stationarity")

    if errors:
        result["status"] = "FAIL"
        result["errors"] = errors
    with (ROOT / "provenance" / "integrity_audit.json").open("w") as f:
        json.dump(result, f, indent=2)
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

