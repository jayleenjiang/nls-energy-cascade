#!/usr/bin/env python3
"""Integrity audit for the analysis-only sector follow-up."""

import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EXPECTED = {
    ROOT / "analysis/driven_T2_T8_autocorrelations.csv": "a7c17121d6d0cc77a89162e31c5caa361e14f12a4d8698163f7e4e8ae1164c5e",
    ROOT / "analysis/equilibrium_T5_T5_autocorrelations.csv": "c36dee1b2da89d8a193cce009ca9ca5775a357ef555c5d7fb3cba9854eb60223",
    ROOT / "analysis/driven_T2_T8_trajectory_correlations.npz": "9bb7e2c8f9a40a373867ee0bb18d7480178ae884c94d21438a0b43b35c7fff81",
    ROOT / "analysis/equilibrium_T5_T5_trajectory_correlations.npz": "367cc1e99a2f22a353aebd49e28f212bbdcdebc534e630713a0cc50c079a31f8",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


input_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in EXPECTED}
fits = read_csv(HERE / "results/sign_changing_early_damped_fits.csv")
curves = read_csv(HERE / "results/sign_changing_fit_curves.csv")
even = json.loads((HERE / "results/even_candidate_common_plateaus.json").read_text())

finite_fields = (
    "lambda_R", "lambda_R_ci_low", "lambda_R_ci_high",
    "lambda_I", "lambda_I_ci_low", "lambda_I_ci_high", "period",
    "delta_AIC_damped_minus_exponential", "mean_correlation_first_zero",
    "fitted_first_zero",
)

checks = {
    "input_hashes_match": all(input_hashes[str(p.relative_to(ROOT))] == h for p, h in EXPECTED.items()),
    "four_fit_rows": len(fits) == 4,
    "384_curve_rows": len(curves) == 384,
    "all_fit_windows_fixed": all(float(r["fit_start"]) == 0.05 and float(r["fit_end"]) == 1.0 for r in fits),
    "all_bootstraps_accepted": all(int(r["bootstrap_accept"]) == 2000 for r in fits),
    "all_reported_fit_values_finite": all(math.isfinite(float(r[k])) for r in fits for k in finite_fields),
    "two_even_group_rows": len(even) == 2,
    "no_simulator_reference_in_analysis": "NLS_stationary_autocorr" not in (HERE / "analyze_sectors.py").read_text(),
}

audit = {
    "status": "PASS" if all(checks.values()) else "FAIL",
    "checks": checks,
    "input_sha256": input_hashes,
    "output_sha256": {
        str(p.relative_to(HERE)): sha256(p)
        for p in sorted((HERE / "results").glob("*")) if p.is_file()
    },
}
(HERE / "integrity_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print(json.dumps(audit, indent=2))
