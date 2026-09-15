#!/usr/bin/env python3
"""Integrity and provenance audit for the analysis-only EDMD extension."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
AN = HERE / "analysis"
INPUT = HERE.parent / "spectral_gap_controlled_2026-09-07" / "raw"
PROV = HERE / "provenance"
PROV.mkdir(exist_ok=True)
CASES = ("driven_dt1e-3", "driven_dt2p5e-4",
         "equilibrium_dt1e-3", "equilibrium_dt2p5e-4")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def csv_rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


audit = {"status": "PASS", "checks": [], "cases": {}}
input_hashes = []
for case in CASES:
    files = sorted((INPUT / case).glob("stream_*.f32"))
    size_ok = len(files) == 64 and all(p.stat().st_size == 100001 * 12 * 4 for p in files)
    for p in files:
        input_hashes.append((sha256(p), str(p.resolve())))
    osc = csv_rows(AN / f"{case}_bootstrap_oscillatory_raw.csv")
    real = csv_rows(AN / f"{case}_bootstrap_real_raw.csv")
    spectra = csv_rows(AN / f"{case}_spectra.csv")
    gram = csv_rows(AN / f"{case}_gram_diagnostics.csv")
    z = np.load(AN / f"{case}_aggregate_matrices.npz")
    matrix_ok = (z["A"].shape == (519, 519) and z["B"].shape == (6, 519, 519)
                 and z["b"].shape == (519, 4) and z["origins"].shape == (2048,))
    osc_ids = sorted(int(r["replicate"]) for r in osc)
    osc_ok = len(osc) == 500 and osc_ids == list(range(500))
    real_modes = sorted(set(int(r["reference_mode"]) for r in real))
    real_ok = all(len([r for r in real if int(r["reference_mode"]) == m]) == 500
                  for m in real_modes)
    constant_ok = all(float(r["trivial_error"]) <= 1e-8 for r in gram)
    finite_spectra = all(np.isfinite(float(r["lambda_real"]))
                         and np.isfinite(float(r["lambda_imag"])) for r in spectra)
    ok = size_ok and matrix_ok and osc_ok and real_ok and constant_ok and finite_spectra
    audit["cases"][case] = {
        "streams": len(files), "stream_size_ok": size_ok,
        "osc_bootstrap_rows": len(osc), "osc_bootstrap_ok": osc_ok,
        "real_bootstrap_rows": len(real), "real_modes": real_modes,
        "real_bootstrap_ok": real_ok, "aggregate_matrix_shapes_ok": matrix_ok,
        "spectrum_rows": len(spectra), "finite_spectra": finite_spectra,
        "gram_rows": len(gram), "constant_mode_ok": constant_ok,
        "max_trivial_error": max(float(r["trivial_error"]) for r in gram),
        "overall": ok,
    }
    if not ok:
        audit["status"] = "FAIL"

with (PROV / "input_stream_hashes.sha256").open("w") as f:
    for digest, path in input_hashes:
        f.write(f"{digest}  {path}\n")

output_hashes = []
for p in sorted(AN.glob("*")):
    if p.is_file():
        output_hashes.append((sha256(p), str(p.resolve())))
with (PROV / "analysis_output_hashes.sha256").open("w") as f:
    for digest, path in output_hashes:
        f.write(f"{digest}  {path}\n")

audit["input_hash_count"] = len(input_hashes)
audit["analysis_output_hash_count"] = len(output_hashes)
audit["analysis_source_sha256"] = sha256(HERE / "run_extended_edmd.py")
audit["summary_source_sha256"] = sha256(HERE / "summarize_extended_edmd.py")
audit["figure_source_sha256"] = sha256(HERE / "make_figures.py")
with (HERE / "AUDIT.json").open("w") as f:
    json.dump(audit, f, indent=2)
print(json.dumps(audit, indent=2))
