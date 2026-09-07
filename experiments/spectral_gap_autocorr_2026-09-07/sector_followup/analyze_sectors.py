#!/usr/bin/env python3
"""Analysis-only symmetry/sign-changing follow-up for saved n=3 correlations."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "results"
CASES = ("driven_T2_T8", "equilibrium_T5_T5")
EVEN_CANDIDATES = (
    "cos_theta1",
    "cos_theta3",
    "cos_theta1_minus_theta3",
    "I2",
    "I1_plus_I3",
)
SIGN_CHANGING = ("sin_theta1", "I1_minus_I3")
FIT_START = 0.05
FIT_END = 1.00
BOOTSTRAPS = 2000
BOOT_SEED = 2026090713
MIN_BOOT_ACCEPT_FRAC = 0.95
DELTA_AIC_GATE = -10.0
HALF_CYCLE_THRESHOLD = math.pi / (FIT_END - FIT_START)


def read_dicts(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_dicts(path: Path, rows):
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def damped_prediction(p, t):
    lr, li, a, b = p
    return np.exp(lr * t) * (a * np.cos(li * t) + b * np.sin(li * t))


def damped_residual(p, t, y, se):
    return (damped_prediction(p, t) - y) / np.maximum(se, 1.0e-12)


def exponential_residual(p, t, y, se):
    lr, a = p
    return (a * np.exp(lr * t) - y) / np.maximum(se, 1.0e-12)


def fit_damped(t, y, se, initial=None):
    amp = max(abs(float(y[0])), 1.0e-8)
    if initial is None:
        starts = [
            np.array([lr, li, math.copysign(amp, y[0]), 0.0])
            for lr in (-0.5, -1.0, -2.0, -4.0, -8.0)
            for li in (0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0)
        ]
    else:
        starts = [np.asarray(initial, dtype=float)]
    lower = np.array([-20.0, 0.0, -1000.0, -1000.0])
    upper = np.array([0.0, 20.0, 1000.0, 1000.0])
    best = None
    for start in starts:
        try:
            fit = least_squares(
                damped_residual,
                np.clip(start, lower + 1.0e-10, upper - 1.0e-10),
                args=(t, y, se),
                bounds=(lower, upper),
                max_nfev=5000,
            )
        except Exception:
            continue
        rss = float(np.sum(fit.fun**2))
        if fit.success and np.all(np.isfinite(fit.x)) and (best is None or rss < best[0]):
            best = (rss, fit.x)
    return best


def fit_exponential(t, y, se):
    starts = [(lr, y[0]) for lr in (-0.5, -1.0, -2.0, -4.0, -8.0)]
    best = None
    for start in starts:
        try:
            fit = least_squares(
                exponential_residual,
                start,
                args=(t, y, se),
                bounds=([-20.0, -1000.0], [0.0, 1000.0]),
                max_nfev=5000,
            )
        except Exception:
            continue
        rss = float(np.sum(fit.fun**2))
        if fit.success and np.all(np.isfinite(fit.x)) and (best is None or rss < best[0]):
            best = (rss, fit.x)
    return best


def first_zero(t, y):
    for i in range(len(t) - 1):
        if y[i] == 0.0:
            return float(t[i])
        if y[i] * y[i + 1] < 0.0:
            return float(t[i] - y[i] * (t[i + 1] - t[i]) / (y[i + 1] - y[i]))
    return math.nan


def even_common(case):
    rows = read_dicts(ROOT / "analysis" / f"{case}_plateaus.csv")
    selected = [r for r in rows if r["observable"] in EVEN_CANDIDATES]
    accepted = [r for r in selected if r["plateau_start"] != ""]
    result = {
        "case": case,
        "group": "even-candidate",
        "requested_observables": list(EVEN_CANDIDATES),
        "accepted_observables": [r["observable"] for r in accepted],
        "unresolved_observables": [r["observable"] for r in selected if r["plateau_start"] == ""],
        "status": "NO COMMON PLATEAU",
    }
    if accepted:
        rates = np.asarray([float(r["lambda_R"]) for r in accepted])
        med = float(np.median(rates))
        result.update({
            "rate_min": float(np.min(rates)),
            "rate_max": float(np.max(rates)),
            "absolute_spread": float(np.ptp(rates)),
            "relative_spread_over_abs_median": float(np.ptp(rates) / abs(med)),
        })
    if len(accepted) < 3:
        result["failure_reason"] = "fewer than three individually resolved observables"
        return result
    start = max(float(r["plateau_start"]) for r in accepted)
    end = min(float(r["plateau_end"]) for r in accepted)
    ci_low = max(float(r["ci_low"]) for r in accepted)
    ci_high = min(float(r["ci_high"]) for r in accepted)
    time_ok = start <= end
    spread_ok = result["relative_spread_over_abs_median"] <= 0.20
    ci_ok = ci_low <= ci_high
    result.update({
        "common_time_start": start,
        "common_time_end": end,
        "time_overlap": time_ok,
        "ci_intersection_low": ci_low,
        "ci_intersection_high": ci_high,
        "ci_overlap": ci_ok,
        "rate_spread_ok": spread_ok,
    })
    if time_ok and spread_ok and ci_ok:
        result["status"] = "COMMON PLATEAU"
        ses = np.asarray([
            (float(r["ci_high"]) - float(r["ci_low"])) / (2.0 * 1.96)
            for r in accepted
        ])
        weights = 1.0 / np.maximum(ses, 1.0e-12) ** 2
        result["reported_rate"] = float(np.sum(weights * np.asarray([
            float(r["lambda_R"]) for r in accepted
        ])) / np.sum(weights))
        result["reported_ci_intersection"] = [ci_low, ci_high]
    else:
        failures = []
        if not time_ok:
            failures.append("no common time interval")
        if not spread_ok:
            failures.append("rate spread exceeds 20%")
        if not ci_ok:
            failures.append("95% confidence intervals do not intersect")
        result["failure_reason"] = "; ".join(failures)
    return result


def fit_sign_changing(case, case_index):
    data = np.load(ROOT / "analysis" / f"{case}_trajectory_correlations.npz")
    names = [str(x) for x in data["observables"]]
    times = np.asarray(data["times"], dtype=float)
    corr = np.asarray(data["C_by_trajectory"], dtype=float)
    ids = np.where((times >= FIT_START - 1.0e-12) & (times <= FIT_END + 1.0e-12))[0]
    t = times[ids]
    diag = {r["observable"]: r for r in read_dicts(ROOT / "analysis" / f"{case}_diagnostics.csv")}
    rows = []
    curves = []
    bootstrap_parameters = {}
    for obs_index, obs in enumerate(SIGN_CHANGING):
        q = names.index(obs)
        per_traj = corr[:, q, ids]
        y = np.mean(per_traj, axis=0)
        se = np.std(per_traj, axis=0, ddof=1) / math.sqrt(per_traj.shape[0])
        damped = fit_damped(t, y, se)
        null = fit_exponential(t, y, se)
        if damped is None or null is None:
            rows.append({"case": case, "observable": obs, "status": "FIT_FAILED"})
            continue
        rss, p = damped
        null_rss, null_p = null
        n = len(t)
        aic = n * math.log(max(rss / n, 1.0e-300)) + 2 * 4
        null_aic = n * math.log(max(null_rss / n, 1.0e-300)) + 2 * 2
        delta_aic = aic - null_aic
        observed_zero = first_zero(t, y)
        dense_t = np.linspace(FIT_START, FIT_END, 20001)
        fitted_zero = first_zero(dense_t, damped_prediction(p, dense_t))

        rng = np.random.default_rng(BOOT_SEED + 100 * case_index + obs_index)
        vals = []
        zeros = []
        for _ in range(BOOTSTRAPS):
            sample = rng.integers(0, per_traj.shape[0], size=per_traj.shape[0])
            yb = np.mean(per_traj[sample], axis=0)
            result = fit_damped(t, yb, se, initial=p)
            if result is None:
                continue
            _, pb = result
            vals.append(pb[:2])
            zeros.append(first_zero(dense_t, damped_prediction(pb, dense_t)))
        vals = np.asarray(vals, dtype=float)
        zeros = np.asarray(zeros, dtype=float)
        bootstrap_parameters[obs] = vals
        lr_low, li_low = np.percentile(vals, 2.5, axis=0)
        lr_high, li_high = np.percentile(vals, 97.5, axis=0)
        finite_zeros = zeros[np.isfinite(zeros)]
        zero_low, zero_high = (
            np.percentile(finite_zeros, [2.5, 97.5])
            if finite_zeros.size else (math.nan, math.nan)
        )
        resolved = (
            len(vals) >= MIN_BOOT_ACCEPT_FRAC * BOOTSTRAPS
            and li_low > 0.0
            and delta_aic <= DELTA_AIC_GATE
            and p[1] >= HALF_CYCLE_THRESHOLD
            and p[1] < 20.0 - 1.0e-8
        )
        rows.append({
            "case": case,
            "observable": obs,
            "status": "FIT",
            "fit_start": FIT_START,
            "fit_end": FIT_END,
            "n_time_points": n,
            "lambda_R": p[0],
            "lambda_R_ci_low": lr_low,
            "lambda_R_ci_high": lr_high,
            "lambda_I": p[1],
            "lambda_I_ci_low": li_low,
            "lambda_I_ci_high": li_high,
            "period": 2.0 * math.pi / p[1] if p[1] > 0.0 else math.inf,
            "period_ci_low": 2.0 * math.pi / li_high if li_high > 0.0 else math.nan,
            "period_ci_high": 2.0 * math.pi / li_low if li_low > 0.0 else math.nan,
            "delta_AIC_damped_minus_exponential": delta_aic,
            "bootstrap_accept": len(vals),
            "half_cycle_threshold": HALF_CYCLE_THRESHOLD,
            "oscillation_resolved": int(resolved),
            "mean_correlation_first_zero": observed_zero,
            "fitted_first_zero": fitted_zero,
            "fitted_first_zero_ci_low": zero_low,
            "fitted_first_zero_ci_high": zero_high,
            "first_supported_negative_time": diag[obs]["first_supported_negative_time"],
        })
        pred = damped_prediction(p, t)
        null_pred = null_p[1] * np.exp(null_p[0] * t)
        for ti, yi, sei, pi, ni in zip(t, y, se, pred, null_pred):
            curves.append({
                "case": case,
                "observable": obs,
                "time": ti,
                "C": yi,
                "SE_C": sei,
                "damped_fit": pi,
                "real_exponential_fit": ni,
            })
    return rows, curves, bootstrap_parameters


def intervals_overlap(a_low, a_high, b_low, b_high):
    return max(float(a_low), float(b_low)) <= min(float(a_high), float(b_high))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    even = [even_common(case) for case in CASES]
    odd_rows = []
    curve_rows = []
    agreement = []
    for case_index, case in enumerate(CASES):
        rows, curves, _ = fit_sign_changing(case, case_index)
        odd_rows.extend(rows)
        curve_rows.extend(curves)
        good = {r["observable"]: r for r in rows if r["status"] == "FIT"}
        if len(good) == 2:
            a, b = (good[x] for x in SIGN_CHANGING)
            agreement.append({
                "case": case,
                "lambda_R_ci_overlap": intervals_overlap(
                    a["lambda_R_ci_low"], a["lambda_R_ci_high"],
                    b["lambda_R_ci_low"], b["lambda_R_ci_high"],
                ),
                "lambda_I_ci_overlap": intervals_overlap(
                    a["lambda_I_ci_low"], a["lambda_I_ci_high"],
                    b["lambda_I_ci_low"], b["lambda_I_ci_high"],
                ),
            })
    write_dicts(OUT / "sign_changing_early_damped_fits.csv", odd_rows)
    write_dicts(OUT / "sign_changing_fit_curves.csv", curve_rows)
    with (OUT / "even_candidate_common_plateaus.json").open("w") as f:
        json.dump(even, f, indent=2)
    with (OUT / "sector_summary.json").open("w") as f:
        json.dump({
            "fit_window": [FIT_START, FIT_END],
            "bootstrap_replicates": BOOTSTRAPS,
            "half_cycle_threshold": HALF_CYCLE_THRESHOLD,
            "even_candidate_results": even,
            "sign_changing_pair_agreement": agreement,
            "symmetry_label_limitation": (
                "sin(theta1) is not exchange odd by itself; the saved observable set "
                "does not span the exact exchange sectors"
            ),
        }, f, indent=2)
    print(json.dumps({"even": even, "odd": odd_rows, "agreement": agreement}, indent=2))


if __name__ == "__main__":
    main()
