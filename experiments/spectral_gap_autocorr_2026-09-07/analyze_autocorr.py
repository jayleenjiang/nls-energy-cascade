#!/usr/bin/env python3
"""Frozen stationary-autocorrelation analysis for the n=3 NLS experiment."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


OBS = [
    "cos_theta1",
    "sin_theta1",
    "cos_theta3",
    "cos_theta1_minus_theta3",
    "I2",
    "I1_plus_I3",
    "I1_minus_I3",
]
NOBS = len(OBS)
SNAPSHOT = 0.01
MEASURE = 1000.0
MAX_T = 10.0
MAX_LAG = int(round(MAX_T / SNAPSHOT))
BOOTSTRAPS = 2000
DAMPED_BOOTSTRAPS = 500
BOOT_SEED = 2026090711
DAMPED_BOOT_SEED = 2026090712
T_MIN = 1.0
SNR_MIN = 3.0
VARIATION_FRAC = 0.20
MIN_POINTS = 8
MIN_SPAN_RATIO = 1.5


def next_pow_two(x: int) -> int:
    return 1 << (x - 1).bit_length()


def autocov_fft(x: np.ndarray, mean: float, max_lag: int) -> np.ndarray:
    y = np.asarray(x, dtype=np.float64) - mean
    n = y.size
    nfft = next_pow_two(2 * n)
    z = np.fft.rfft(y, nfft)
    ac = np.fft.irfft(z * np.conjugate(z), nfft)[: max_lag + 1]
    return ac / np.arange(n, n - max_lag - 1, -1, dtype=np.float64)


def lag_grid() -> np.ndarray:
    positive = np.unique(np.rint(np.geomspace(1, MAX_LAG, 121)).astype(int))
    return np.concatenate(([0], positive))


def local_rate(c: np.ndarray, times: np.ndarray) -> np.ndarray:
    out = np.full(c.shape, np.nan, dtype=np.float64)
    for i in range(2, len(c) - 2):
        if c[i - 2] > 0.0 and c[i + 2] > 0.0:
            out[i] = (math.log(c[i + 2]) - math.log(c[i - 2])) / (
                times[i + 2] - times[i - 2]
            )
    return out


def select_plateau(times, c, se, rates):
    valid = (
        (times >= T_MIN)
        & np.isfinite(rates)
        & (rates < 0.0)
        & (c > 0.0)
        & (c / np.maximum(se, np.finfo(float).tiny) >= SNR_MIN)
    )
    candidates = []
    n = len(times)
    for i in range(n):
        if not valid[i]:
            continue
        for j in range(i + MIN_POINTS - 1, n):
            if not np.all(valid[i : j + 1]):
                break
            if times[j] / times[i] < MIN_SPAN_RATIO:
                continue
            rr = rates[i : j + 1]
            med = float(np.median(rr))
            if med >= 0.0:
                continue
            if float(np.max(rr) - np.min(rr)) <= VARIATION_FRAC * abs(med):
                candidates.append((math.log(times[j] / times[i]), times[i], i, j))
    if not candidates:
        return None
    # Maximum log span; on a tie, choose the later interval.
    _, _, i, j = max(candidates, key=lambda x: (x[0], x[1]))
    return i, j


def weighted_log_slope(times, c, se, i, j):
    t = times[i : j + 1]
    y = np.log(c[i : j + 1])
    sy = np.maximum(se[i : j + 1] / c[i : j + 1], 1.0e-12)
    w = 1.0 / (sy * sy)
    X = np.column_stack((np.ones_like(t), t))
    beta = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * y))
    return float(beta[1]), float(beta[0])


def ips_tau(c: np.ndarray) -> tuple[float, bool]:
    rho = c / c[0]
    last = 0
    for k in range(1, len(rho) - 1, 2):
        if not np.isfinite(rho[k] + rho[k + 1]) or rho[k] + rho[k + 1] <= 0.0:
            break
        last = k + 1
    truncated = last == len(rho) - 1 or last == len(rho) - 2
    tau = SNAPSHOT * (0.5 + float(np.sum(rho[1 : last + 1])))
    return max(tau, 0.5 * SNAPSHOT), truncated


def damped_residual(params, t, y, se):
    lr, li, a, b = params
    pred = np.exp(lr * t) * (a * np.cos(li * t) + b * np.sin(li * t))
    return (pred - y) / np.maximum(se, 1.0e-12)


def exponential_residual(params, t, y, se):
    lr, a = params
    pred = a * np.exp(lr * t)
    return (pred - y) / np.maximum(se, 1.0e-12)


def exponential_fit(t, y, se):
    starts = [(-0.25, y[0]), (-0.5, y[0]), (-1.0, y[0]), (-2.0, y[0]), (-4.0, y[0])]
    best = None
    for start in starts:
        try:
            fit = least_squares(
                exponential_residual,
                start,
                args=(t, y, se),
                bounds=([-20.0, -1000.0], [0.0, 1000.0]),
                max_nfev=3000,
            )
        except Exception:
            continue
        rss = float(np.sum(fit.fun * fit.fun))
        if fit.success and (best is None or rss < best[0]):
            best = (rss, fit.x)
    return best


def damped_fit(t, y, se, initial=None):
    if initial is not None:
        starts = [initial]
    else:
        starts = []
        amp = max(abs(y[0]), 1.0e-6)
        for lr in (-0.25, -0.5, -1.0, -2.0, -4.0):
            for li in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
                starts.append(np.array([lr, li, math.copysign(amp, y[0]), 0.0]))
    best = None
    lower = np.array([-20.0, 0.0, -1000.0, -1000.0])
    upper = np.array([0.0, 20.0, 1000.0, 1000.0])
    for start in starts:
        try:
            fit = least_squares(
                damped_residual,
                np.clip(start, lower + 1e-10, upper - 1e-10),
                args=(t, y, se),
                bounds=(lower, upper),
                max_nfev=3000,
            )
        except Exception:
            continue
        rss = float(np.sum(fit.fun * fit.fun))
        if fit.success and (best is None or rss < best[0]):
            best = (rss, fit.x)
    return None if best is None else best[1]


def load_case(case_dir: Path):
    files = sorted(case_dir.glob("stream_*.f32"))
    if len(files) != 64:
        raise RuntimeError(f"{case_dir}: expected 64 streams, found {len(files)}")
    sizes = {p.stat().st_size for p in files}
    if len(sizes) != 1:
        raise RuntimeError(f"{case_dir}: unequal stream file sizes")
    n_snap = next(iter(sizes)) // (4 * NOBS)
    expected = int(round(MEASURE / SNAPSHOT)) + 1
    if n_snap != expected or next(iter(sizes)) != n_snap * NOBS * 4:
        raise RuntimeError(f"{case_dir}: invalid binary shape {n_snap}")
    arrays = [np.memmap(p, mode="r", dtype="<f4", shape=(n_snap, NOBS)) for p in files]
    return files, arrays, n_snap


def analyze_case(case_name: str, case_dir: Path, outdir: Path):
    files, arrays, n_snap = load_case(case_dir)
    means = np.zeros(NOBS)
    for x in arrays:
        means += np.sum(x, axis=0, dtype=np.float64)
    means /= len(arrays) * n_snap

    corr = np.empty((len(arrays), NOBS, MAX_LAG + 1), dtype=np.float64)
    corr_tail = np.empty_like(corr)
    tail_start = n_snap // 4
    for r, x in enumerate(arrays):
        for q in range(NOBS):
            corr[r, q] = autocov_fft(x[:, q], means[q], MAX_LAG)
            corr_tail[r, q] = autocov_fft(x[tail_start:, q], means[q], MAX_LAG)

    grid = lag_grid()
    times = grid * SNAPSHOT
    c_grid = corr[:, :, grid]
    ct_grid = corr_tail[:, :, grid]
    mean_c = np.mean(c_grid, axis=0)
    se_c = np.std(c_grid, axis=0, ddof=1) / math.sqrt(len(arrays))
    mean_ct = np.mean(ct_grid, axis=0)
    rates = np.stack([local_rate(mean_c[q], times) for q in range(NOBS)])

    rng = np.random.default_rng(BOOT_SEED)
    boot_indices = rng.integers(0, len(arrays), size=(BOOTSTRAPS, len(arrays)))
    boot_c = np.empty((BOOTSTRAPS, NOBS, len(grid)), dtype=np.float64)
    boot_rates = np.empty_like(boot_c)
    for b, idx in enumerate(boot_indices):
        boot_c[b] = np.mean(c_grid[idx], axis=0)
        for q in range(NOBS):
            boot_rates[b, q] = local_rate(boot_c[b, q], times)
    rate_lo = np.nanpercentile(boot_rates, 2.5, axis=0)
    rate_hi = np.nanpercentile(boot_rates, 97.5, axis=0)

    dense_mean = np.mean(corr, axis=0)
    summaries = []
    plateau_rows = []
    damped_rows = []
    plateau_records = []

    for q, obs in enumerate(OBS):
        c0 = float(mean_c[q, 0])
        stationarity_rms = float(
            np.sqrt(np.mean((mean_c[q, times <= 5] - mean_ct[q, times <= 5]) ** 2))
            / abs(c0)
        )
        tau, tau_truncated = ips_tau(dense_mean[q])
        n_eff = len(arrays) * MEASURE / (2.0 * tau)

        rate_support = (
            (times >= T_MIN)
            & np.isfinite(rates[q])
            & (rates[q] < 0.0)
            & (mean_c[q] > 0.0)
            & (mean_c[q] / np.maximum(se_c[q], np.finfo(float).tiny) >= SNR_MIN)
        )
        support_ids = np.where(rate_support)[0]
        if support_ids.size:
            terminal = [int(support_ids[-1])]
            for idx in support_ids[-2::-1]:
                if idx == terminal[-1] - 1:
                    terminal.append(int(idx))
                else:
                    break
            terminal = np.asarray(terminal[::-1][-5:], dtype=int)
            last_supported_time = float(times[support_ids[-1]])
            last_supported_rate = float(rates[q, support_ids[-1]])
            late_drift = (
                float(np.polyfit(times[terminal], rates[q, terminal], 1)[0])
                if terminal.size >= 3 else np.nan
            )
        else:
            last_supported_time = last_supported_rate = late_drift = np.nan

        interval = select_plateau(times, mean_c[q], se_c[q], rates[q])
        if interval is None:
            plateau_records.append(None)
            plateau_rows.append({
                "case": case_name, "observable": obs, "C0": c0,
                "plateau_start": "", "plateau_end": "", "lambda_R": "",
                "ci_low": "", "ci_high": "", "bootstrap_accept": 0,
                "tau_int": tau, "tau_truncated_at_t10": int(tau_truncated),
                "effective_count": n_eff, "stationarity_rms_over_C0": stationarity_rms,
                "stationarity_pass": int(stationarity_rms <= 0.05),
                "last_supported_rate_time": last_supported_time,
                "last_supported_lambda_eff": last_supported_rate,
                "late_dlambda_dt_last5": late_drift,
            })
        else:
            i, j = interval
            slope, intercept = weighted_log_slope(times, mean_c[q], se_c[q], i, j)
            slopes = []
            for b in range(BOOTSTRAPS):
                cb = boot_c[b, q]
                if np.all(cb[i : j + 1] > 0.0):
                    try:
                        sb, _ = weighted_log_slope(times, cb, se_c[q], i, j)
                        if np.isfinite(sb): slopes.append(sb)
                    except np.linalg.LinAlgError:
                        pass
            lo, hi = (np.percentile(slopes, [2.5, 97.5]) if slopes else (np.nan, np.nan))
            rec = {
                "i": i, "j": j, "start": float(times[i]), "end": float(times[j]),
                "lambda": slope, "ci_low": float(lo), "ci_high": float(hi),
            }
            plateau_records.append(rec)
            plateau_rows.append({
                "case": case_name, "observable": obs, "C0": c0,
                "plateau_start": times[i], "plateau_end": times[j], "lambda_R": slope,
                "ci_low": lo, "ci_high": hi, "bootstrap_accept": len(slopes),
                "tau_int": tau, "tau_truncated_at_t10": int(tau_truncated),
                "effective_count": n_eff, "stationarity_rms_over_C0": stationarity_rms,
                "stationarity_pass": int(stationarity_rms <= 0.05),
                "last_supported_rate_time": last_supported_time,
                "last_supported_lambda_eff": last_supported_rate,
                "late_dlambda_dt_last5": late_drift,
            })

        snr_abs = np.abs(mean_c[q]) / np.maximum(se_c[q], np.finfo(float).tiny)
        eligible = np.where((times >= T_MIN) & (snr_abs >= SNR_MIN))[0]
        if eligible.size == 0:
            damped_rows.append({"case": case_name, "observable": obs, "status": "UNRESOLVED"})
        else:
            last = int(eligible[-1])
            ids = np.where((times >= T_MIN) & (np.arange(len(times)) <= last))[0]
            if ids.size < 12:
                damped_rows.append({"case": case_name, "observable": obs, "status": "UNRESOLVED"})
            else:
                fit = damped_fit(times[ids], mean_c[q, ids], se_c[q, ids])
                if fit is None:
                    damped_rows.append({"case": case_name, "observable": obs, "status": "FIT_FAILED"})
                else:
                    full_rss = float(np.sum(damped_residual(fit, times[ids], mean_c[q, ids], se_c[q, ids]) ** 2))
                    null = exponential_fit(times[ids], mean_c[q, ids], se_c[q, ids])
                    null_rss = np.nan if null is None else float(null[0])
                    n_fit = len(ids)
                    aic_full = n_fit * math.log(max(full_rss / n_fit, 1e-300)) + 2 * 4
                    aic_null = (
                        np.nan if not np.isfinite(null_rss)
                        else n_fit * math.log(max(null_rss / n_fit, 1e-300)) + 2 * 2
                    )
                    delta_aic = aic_full - aic_null
                    rayleigh_half_cycle = math.pi / (times[ids[-1]] - times[ids[0]])
                    drng = np.random.default_rng(DAMPED_BOOT_SEED + q)
                    vals = []
                    for _ in range(DAMPED_BOOTSTRAPS):
                        idx = drng.integers(0, len(arrays), size=len(arrays))
                        cb = np.mean(c_grid[idx, q], axis=0)
                        fb = damped_fit(times[ids], cb[ids], se_c[q, ids], initial=fit)
                        if fb is not None and np.all(np.isfinite(fb)):
                            vals.append(fb[:2])
                    vals = np.asarray(vals)
                    if len(vals):
                        lrlo, lilo = np.percentile(vals, 2.5, axis=0)
                        lrhi, lihi = np.percentile(vals, 97.5, axis=0)
                    else:
                        lrlo = lrhi = lilo = lihi = np.nan
                    damped_rows.append({
                        "case": case_name, "observable": obs, "status": "FIT",
                        "fit_start": times[ids[0]], "fit_end": times[ids[-1]],
                        "lambda_R": fit[0], "lambda_R_ci_low": lrlo,
                        "lambda_R_ci_high": lrhi, "lambda_I": fit[1],
                        "lambda_I_ci_low": lilo, "lambda_I_ci_high": lihi,
                        "bootstrap_accept": len(vals),
                        "frequency_half_cycle_threshold": rayleigh_half_cycle,
                        "delta_AIC_damped_minus_exponential": delta_aic,
                        "lambda_I_resolved": int(
                            lilo > 0
                            and fit[1] >= rayleigh_half_cycle
                            and delta_aic <= -10.0
                            and fit[1] < 20 - 1e-8
                        ),
                    })

        dense_se = np.std(corr[:, q], axis=0, ddof=1) / math.sqrt(len(arrays))
        diff = corr[:, q, 1:] - corr[:, q, :-1]
        dmean = np.mean(diff, axis=0)
        dse = np.std(diff, axis=0, ddof=1) / math.sqrt(len(arrays))
        significant_increases = int(np.sum(dmean > 1.96 * dse))
        supported_negative = np.where(dense_mean[q] < -3.0 * dense_se)[0]
        summaries.append({
            "case": case_name, "observable": obs, "mean_g": means[q], "C0": c0,
            "significant_monotonicity_violations": significant_increases,
            "supported_negative_lag_count": len(supported_negative),
            "first_supported_negative_time": (
                float(supported_negative[0] * SNAPSHOT) if len(supported_negative) else ""
            ),
        })

    common = common_plateau(case_name, OBS, plateau_records)

    outdir.mkdir(parents=True, exist_ok=True)
    write_dict_csv(outdir / f"{case_name}_plateaus.csv", plateau_rows)
    write_dict_csv(outdir / f"{case_name}_damped_fits.csv", damped_rows)
    write_dict_csv(outdir / f"{case_name}_diagnostics.csv", summaries)
    with (outdir / f"{case_name}_common_plateau.json").open("w") as f:
        json.dump(common, f, indent=2)

    with (outdir / f"{case_name}_autocorrelations.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "observable", "lag_index", "time", "C", "SE_C", "SNR_C",
            "C_discard_first_quarter", "lambda_eff", "lambda_ci_low", "lambda_ci_high",
        ])
        for q, obs in enumerate(OBS):
            for m, lag in enumerate(grid):
                w.writerow([
                    obs, lag, times[m], mean_c[q, m], se_c[q, m],
                    mean_c[q, m] / max(se_c[q, m], np.finfo(float).tiny),
                    mean_ct[q, m], rates[q, m], rate_lo[q, m], rate_hi[q, m],
                ])

    np.savez_compressed(
        outdir / f"{case_name}_trajectory_correlations.npz",
        observables=np.asarray(OBS), lags=np.arange(MAX_LAG + 1),
        times=np.arange(MAX_LAG + 1) * SNAPSHOT, C_by_trajectory=corr,
        C_discard_first_quarter_by_trajectory=corr_tail,
    )
    return {
        "case": case_name, "streams": len(arrays), "snapshots_per_stream": n_snap,
        "means": dict(zip(OBS, means.tolist())), "common_plateau": common,
    }


def common_plateau(case, names, records):
    valid = [(name, r) for name, r in zip(names, records) if r is not None]
    result = {"case": case, "status": "NO COMMON PLATEAU", "eligible": [x[0] for x in valid]}
    if len(valid) < 3:
        return result
    start = max(r["start"] for _, r in valid)
    end = min(r["end"] for _, r in valid)
    rates = np.asarray([r["lambda"] for _, r in valid])
    med = float(np.median(rates))
    rate_spread_ok = float(np.max(rates) - np.min(rates)) <= 0.20 * abs(med)
    ci_low = max(r["ci_low"] for _, r in valid)
    ci_high = min(r["ci_high"] for _, r in valid)
    ci_overlap = ci_low <= ci_high
    time_overlap = start <= end
    result.update({
        "common_time_start": start, "common_time_end": end,
        "rate_min": float(np.min(rates)), "rate_max": float(np.max(rates)),
        "rate_spread_ok": bool(rate_spread_ok), "ci_intersection_low": ci_low,
        "ci_intersection_high": ci_high, "ci_overlap": bool(ci_overlap),
        "time_overlap": bool(time_overlap),
    })
    if rate_spread_ok and ci_overlap and time_overlap:
        result["status"] = "COMMON PLATEAU"
        result["reported_rate"] = med
        result["reported_ci_intersection"] = [ci_low, ci_high]
    return result


def write_dict_csv(path: Path, rows):
    keys = []
    for row in rows:
        for key in row:
            if key not in keys: keys.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args()
    exp = args.experiment.resolve()
    out = exp / "analysis"
    results = []
    for case in ("driven_T2_T8", "equilibrium_T5_T5"):
        results.append(analyze_case(case, exp / "raw" / case, out))
    with (out / "analysis_summary.json").open("w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
