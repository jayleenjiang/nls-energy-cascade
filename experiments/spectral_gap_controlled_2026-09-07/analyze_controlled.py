#!/usr/bin/env python3
"""Frozen controlled-integrator autocorrelation and oscillation analysis."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
LEGACY = ROOT.parent / "spectral_gap_autocorr_2026-09-07" / "analyze_autocorr.py"
LEGACY_SHA256 = "0a379e0e9e0356cf3240c353b14df699727e22453d7ca1908fd9683865b5a5b9"
OBS = [
    "cos_theta1",
    "sin_theta1",
    "cos_theta3",
    "sin_theta3",
    "cos_theta1_minus_theta3",
    "I2",
    "I1_plus_I3",
    "I1_minus_I3",
    "cos_even_sum",
    "cos_odd_diff",
    "sin_even_sum",
    "sin_odd_diff",
]
CASES = [
    ("driven_dt1e-3", "driven", 1.0e-3),
    ("driven_dt2p5e-4", "driven", 2.5e-4),
    ("equilibrium_dt1e-3", "equilibrium", 1.0e-3),
    ("equilibrium_dt2p5e-4", "equilibrium", 2.5e-4),
]
EVEN = [
    "cos_theta1_minus_theta3",
    "I2",
    "I1_plus_I3",
    "cos_even_sum",
    "sin_even_sum",
]
ODD = ["I1_minus_I3", "cos_odd_diff", "sin_odd_diff"]
EARLY_ODD = ODD
FIT_START = 0.05
FIT_END = 1.00
EARLY_BOOTSTRAPS = 2000
EARLY_BOOT_SEED = 2026090813
MIN_BOOT_ACCEPT_FRAC = 0.95
DELTA_AIC_GATE = -10.0
HALF_CYCLE_THRESHOLD = math.pi / (FIT_END - FIT_START)
BOOT_SEED = 2026090811


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_legacy():
    actual = sha256(LEGACY)
    if actual != LEGACY_SHA256:
        raise RuntimeError(f"frozen legacy analyzer hash mismatch: {actual}")
    spec = importlib.util.spec_from_file_location("frozen_legacy_autocorr", LEGACY)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    module.OBS = OBS
    module.NOBS = len(OBS)
    module.BOOT_SEED = BOOT_SEED
    # The production question uses the separately frozen early-window fit below;
    # suppress the legacy late-support damped bootstrap, not the plateau bootstrap.
    module.DAMPED_BOOTSTRAPS = 0
    return module


def read_rows(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_rows(path: Path, rows):
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def record_from_row(row):
    if row["plateau_start"] == "":
        return None
    return {
        "start": float(row["plateau_start"]),
        "end": float(row["plateau_end"]),
        "lambda": float(row["lambda_R"]),
        "ci_low": float(row["ci_low"]),
        "ci_high": float(row["ci_high"]),
    }


def first_zero(t, y):
    for i in range(len(t) - 1):
        if y[i] == 0.0:
            return float(t[i])
        if y[i] * y[i + 1] < 0.0:
            return float(t[i] - y[i] * (t[i + 1] - t[i]) / (y[i + 1] - y[i]))
    return math.nan


def damped_prediction(p, t):
    lr, li, a, b = p
    return np.exp(lr * t) * (a * np.cos(li * t) + b * np.sin(li * t))


def early_odd_fits(legacy, case, case_index, analysis_dir):
    data = np.load(analysis_dir / f"{case}_trajectory_correlations.npz")
    names = [str(x) for x in data["observables"]]
    times = np.asarray(data["times"], dtype=float)
    corr = np.asarray(data["C_by_trajectory"], dtype=float)
    ids = np.where((times >= FIT_START - 1e-12) & (times <= FIT_END + 1e-12))[0]
    t = times[ids]
    output = []
    curves = []
    for obs_index, obs in enumerate(EARLY_ODD):
        q = names.index(obs)
        per_stream = corr[:, q, ids]
        y = np.mean(per_stream, axis=0)
        se = np.std(per_stream, axis=0, ddof=1) / math.sqrt(per_stream.shape[0])
        p = legacy.damped_fit(t, y, se)
        null = legacy.exponential_fit(t, y, se)
        if p is None or null is None:
            output.append({"case": case, "observable": obs, "status": "FIT_FAILED"})
            continue
        rss = float(np.sum(legacy.damped_residual(p, t, y, se) ** 2))
        null_rss = float(null[0])
        n = len(t)
        delta_aic = (
            n * math.log(max(rss / n, 1e-300)) + 2 * 4
            - (n * math.log(max(null_rss / n, 1e-300)) + 2 * 2)
        )
        dense_t = np.linspace(FIT_START, FIT_END, 20001)
        rng = np.random.default_rng(EARLY_BOOT_SEED + 100 * case_index + obs_index)
        vals = []
        zeros = []
        for _ in range(EARLY_BOOTSTRAPS):
            sample = rng.integers(0, per_stream.shape[0], size=per_stream.shape[0])
            yb = np.mean(per_stream[sample], axis=0)
            pb = legacy.damped_fit(t, yb, se, initial=p)
            if pb is not None and np.all(np.isfinite(pb)):
                vals.append(pb[:2])
                zeros.append(first_zero(dense_t, damped_prediction(pb, dense_t)))
        vals = np.asarray(vals, dtype=float)
        zeros = np.asarray(zeros, dtype=float)
        if len(vals):
            lr_lo, li_lo = np.percentile(vals, 2.5, axis=0)
            lr_hi, li_hi = np.percentile(vals, 97.5, axis=0)
        else:
            lr_lo = lr_hi = li_lo = li_hi = math.nan
        finite_zeros = zeros[np.isfinite(zeros)]
        zero_ci = (
            np.percentile(finite_zeros, [2.5, 97.5])
            if len(finite_zeros)
            else (math.nan, math.nan)
        )
        resolved = (
            len(vals) >= MIN_BOOT_ACCEPT_FRAC * EARLY_BOOTSTRAPS
            and li_lo > 0.0
            and delta_aic <= DELTA_AIC_GATE
            and p[1] >= HALF_CYCLE_THRESHOLD
            and p[1] < 20.0 - 1e-8
        )
        output.append({
            "case": case,
            "observable": obs,
            "status": "FIT",
            "fit_start": FIT_START,
            "fit_end": FIT_END,
            "n_time_points": len(t),
            "lambda_R": p[0],
            "lambda_R_ci_low": lr_lo,
            "lambda_R_ci_high": lr_hi,
            "lambda_I": p[1],
            "lambda_I_ci_low": li_lo,
            "lambda_I_ci_high": li_hi,
            "period": 2.0 * math.pi / p[1] if p[1] > 0 else math.inf,
            "period_ci_low": 2.0 * math.pi / li_hi if li_hi > 0 else math.nan,
            "period_ci_high": 2.0 * math.pi / li_lo if li_lo > 0 else math.nan,
            "delta_AIC_damped_minus_exponential": delta_aic,
            "bootstrap_accept": len(vals),
            "half_cycle_threshold": HALF_CYCLE_THRESHOLD,
            "oscillation_resolved": int(resolved),
            "mean_correlation_first_zero": first_zero(t, y),
            "fitted_first_zero": first_zero(dense_t, damped_prediction(p, dense_t)),
            "fitted_first_zero_ci_low": zero_ci[0],
            "fitted_first_zero_ci_high": zero_ci[1],
            "historical_reference_lambda_I": 5.35,
            "historical_reference_period": 2.0 * math.pi / 5.35,
        })
        pred = damped_prediction(p, t)
        null_pred = null[1][1] * np.exp(null[1][0] * t)
        for ti, yi, sei, pi, ni in zip(t, y, se, pred, null_pred):
            curves.append({
                "case": case, "observable": obs, "time": ti,
                "C": yi, "SE_C": sei, "damped_fit": pi,
                "real_exponential_fit": ni,
            })
    return output, curves


def compare_timesteps(plateaus_by_case):
    rows = []
    old = {
        "driven": {
            "cos_theta1": None,
            "sin_theta1": None,
            "cos_theta3": None,
            "sin_theta3": "NOT_SAVED",
            "cos_theta1_minus_theta3": None,
            "I2": (-0.9925655214, -1.0517372600, -0.9336587890),
            "I1_plus_I3": (-1.1088171245, -1.2013330743, -1.0228371013),
            "I1_minus_I3": None,
        },
        "equilibrium": {
            "cos_theta1": (-1.1878330875, -1.3607165857, -1.0159711099),
            "sin_theta1": None,
            "cos_theta3": (-1.1333583808, -1.3365572930, -0.9464618584),
            "sin_theta3": "NOT_SAVED",
            "cos_theta1_minus_theta3": None,
            "I2": (-0.9065189091, -0.9664823188, -0.8441868234),
            "I1_plus_I3": None,
            "I1_minus_I3": None,
        },
    }
    primary = OBS[:8]
    for bath in ("driven", "equilibrium"):
        coarse = {r["observable"]: r for r in plateaus_by_case[f"{bath}_dt1e-3"]}
        fine = {r["observable"]: r for r in plateaus_by_case[f"{bath}_dt2p5e-4"]}
        for obs in primary:
            a = coarse[obs]
            b = fine[obs]
            resolved_a = a["lambda_R"] != ""
            resolved_b = b["lambda_R"] != ""
            ci_overlap = False
            rel_difference = math.nan
            timestep_consistent = False
            if resolved_a and resolved_b:
                la, lb = float(a["lambda_R"]), float(b["lambda_R"])
                ci_overlap = max(float(a["ci_low"]), float(b["ci_low"])) <= min(
                    float(a["ci_high"]), float(b["ci_high"])
                )
                rel_difference = abs(la - lb) / abs(np.median([la, lb]))
                timestep_consistent = ci_overlap and rel_difference <= 0.20
            prior = old[bath].get(obs)
            rows.append({
                "bath_case": bath,
                "observable": obs,
                "previous_lambda_R": "" if prior is None or isinstance(prior, str) else prior[0],
                "previous_ci_low": "" if prior is None or isinstance(prior, str) else prior[1],
                "previous_ci_high": "" if prior is None or isinstance(prior, str) else prior[2],
                "previous_status": (
                    prior if isinstance(prior, str) else ("UNRESOLVED" if prior is None else "RESOLVED")
                ),
                "new_dt1e-3_lambda_R": a["lambda_R"],
                "new_dt1e-3_ci_low": a["ci_low"],
                "new_dt1e-3_ci_high": a["ci_high"],
                "new_dt2p5e-4_lambda_R": b["lambda_R"],
                "new_dt2p5e-4_ci_low": b["ci_low"],
                "new_dt2p5e-4_ci_high": b["ci_high"],
                "ci_overlap_between_timesteps": int(ci_overlap),
                "relative_timestep_difference": rel_difference,
                "timestep_consistent": int(timestep_consistent),
            })
    return rows


def main():
    legacy = load_legacy()
    analysis_dir = ROOT / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    case_summaries = []
    plateaus_by_case = {}
    sector_results = []
    early_rows = []
    early_curves = []
    for case_index, (case, bath, dt) in enumerate(CASES):
        summary = legacy.analyze_case(case, ROOT / "raw" / case, analysis_dir)
        summary.update({"bath_case": bath, "dt": dt})
        case_summaries.append(summary)
        rows = read_rows(analysis_dir / f"{case}_plateaus.csv")
        plateaus_by_case[case] = rows
        by_name = {r["observable"]: record_from_row(r) for r in rows}
        for label, names in (("exchange_even", EVEN), ("exchange_odd", ODD)):
            sector_results.append(
                legacy.common_plateau(case + ":" + label, names, [by_name[n] for n in names])
            )
        fits, curves = early_odd_fits(legacy, case, case_index, analysis_dir)
        early_rows.extend(fits)
        early_curves.extend(curves)

    comparison = compare_timesteps(plateaus_by_case)
    write_rows(analysis_dir / "early_odd_damped_fits.csv", early_rows)
    write_rows(analysis_dir / "early_odd_fit_curves.csv", early_curves)
    write_rows(analysis_dir / "previous_vs_controlled_rates.csv", comparison)
    with (analysis_dir / "controlled_summary.json").open("w") as f:
        json.dump({
            "legacy_analysis_sha256": LEGACY_SHA256,
            "plateau_bootstraps": legacy.BOOTSTRAPS,
            "plateau_bootstrap_seed": BOOT_SEED,
            "early_fit_window": [FIT_START, FIT_END],
            "early_fit_bootstraps": EARLY_BOOTSTRAPS,
            "early_fit_bootstrap_seed": EARLY_BOOT_SEED,
            "case_summaries": case_summaries,
            "sector_common_plateaus": sector_results,
        }, f, indent=2)


if __name__ == "__main__":
    main()
