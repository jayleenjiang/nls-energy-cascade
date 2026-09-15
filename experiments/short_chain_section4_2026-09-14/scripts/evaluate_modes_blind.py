#!/usr/bin/env python3
"""One-shot blind-test evaluation of the frozen Section-4.3 EDMD modes."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from mode_common import (
    ROOT, LAGS, NROWS, SNAPSHOT, feature_matrix, load_state, origin_indices,
    split_paths,
)


OUT = ROOT/"modes"
SELECTION = OUT/"MODE_SELECTION.json"
BOOTSTRAPS = 500
BOOTSTRAP_SEED = 2026091404
FIT_START = 0.05
FIT_END = 1.00
INDEPENDENT_IMAG = {
    "driven_dt1e-3": 5.401,
    "driven_dt2p5e-4": 5.500,
    "equilibrium_dt1e-3": 5.081,
    "equilibrium_dt2p5e-4": 5.264,
}

LEGACY_SOURCE = (
    ROOT.parent/"spectral_gap_autocorr_2026-09-07"/"analyze_autocorr.py"
)
SPEC = importlib.util.spec_from_file_location("section4_frozen_autocorr", LEGACY_SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {LEGACY_SOURCE}")
legacy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(legacy)


def evaluate_q(path: Path, coefficient: np.ndarray, indices: np.ndarray,
               log_mean: np.ndarray, log_sd: np.ndarray,
               chunk: int = 4096) -> np.ndarray:
    state = load_state(path)
    output = np.empty(NROWS, dtype=np.complex128)
    for start in range(0, NROWS, chunk):
        features = feature_matrix(state[start:start+chunk], log_mean, log_sd)
        output[start:start+chunk] = features[:, indices]@coefficient
    return output


def phase_align(coefficient: np.ndarray) -> np.ndarray:
    result = coefficient.copy()
    pivot = int(np.argmax(np.abs(result)))
    result *= np.exp(-1j*np.angle(result[pivot]))
    return result


def test_residual(q: np.ndarray, mu: complex, lag_steps: int) -> tuple[float, float]:
    origins = origin_indices()
    delta = q[:, origins+lag_steps]-mu*q[:, origins]
    return float(np.sum(np.abs(delta)**2)), float(np.sum(np.abs(q[:, origins])**2))


def correlation_by_stream(signals: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    maximum_lag = legacy.MAX_LAG
    mean = float(signals.mean())
    corr = np.asarray([
        legacy.autocov_fft(signal, mean, maximum_lag) for signal in signals
    ])
    grid = legacy.lag_grid()
    return corr, grid, corr[:, grid]


def bootstrap_residual(numerators, denominators, rng):
    values = []
    n = len(numerators)
    for _ in range(BOOTSTRAPS):
        index = rng.integers(0, n, size=n)
        values.append(math.sqrt(numerators[index].sum()/denominators[index].sum()))
    return np.asarray(values)


def fit_real_mode(case: str, cgrid: np.ndarray, grid: np.ndarray, rng):
    times = grid*SNAPSHOT
    mean = cgrid.mean(axis=0)
    se = cgrid.std(axis=0, ddof=1)/math.sqrt(len(cgrid))
    rates = legacy.local_rate(mean, times)
    interval = legacy.select_plateau(times, mean, se, rates)
    result = {
        "case": case, "kind": "real", "plateau_resolved": interval is not None,
    }
    if interval is None:
        return result, mean, se, rates
    left, right = interval
    slope, intercept = legacy.weighted_log_slope(times, mean, se, left, right)
    values = []
    for _ in range(BOOTSTRAPS):
        sample = cgrid[rng.integers(0, len(cgrid), size=len(cgrid))].mean(axis=0)
        if np.all(sample[left:right+1] > 0):
            try:
                value, _ = legacy.weighted_log_slope(times, sample, se, left, right)
                if np.isfinite(value):
                    values.append(value)
            except np.linalg.LinAlgError:
                pass
    lo, hi = np.percentile(values, (2.5, 97.5)) if values else (np.nan, np.nan)
    result.update({
        "plateau_start": float(times[left]), "plateau_end": float(times[right]),
        "lambda_R": slope, "lambda_R_ci": [float(lo), float(hi)],
        "intercept": intercept, "bootstrap_accept": len(values),
    })
    return result, mean, se, rates


def fit_complex_mode(case: str, cgrid: np.ndarray, grid: np.ndarray, rng):
    times = grid*SNAPSHOT
    selected = (times >= FIT_START-1e-12)&(times <= FIT_END+1e-12)
    t = times[selected]
    per_stream = cgrid[:, selected]
    mean = per_stream.mean(axis=0)
    se = per_stream.std(axis=0, ddof=1)/math.sqrt(len(per_stream))
    damped = legacy.damped_fit(t, mean, se)
    exponential = legacy.exponential_fit(t, mean, se)
    result = {"case": case, "kind": "complex", "fit_start": FIT_START,
              "fit_end": FIT_END, "fit_resolved": False}
    if damped is None or exponential is None:
        return result, mean, se, np.full_like(mean, np.nan)
    rss = float(np.sum(legacy.damped_residual(damped, t, mean, se)**2))
    null_rss = float(exponential[0])
    n = len(t)
    delta_aic = n*math.log(max(rss/n, 1e-300))+8 - (
        n*math.log(max(null_rss/n, 1e-300))+4
    )
    values = []
    for _ in range(BOOTSTRAPS):
        sample = per_stream[rng.integers(0, len(per_stream), size=len(per_stream))].mean(axis=0)
        fitted = legacy.damped_fit(t, sample, se, initial=damped)
        if fitted is not None and np.all(np.isfinite(fitted)):
            values.append(fitted[:2])
    values = np.asarray(values)
    if len(values):
        low = np.percentile(values, 2.5, axis=0)
        high = np.percentile(values, 97.5, axis=0)
    else:
        low = high = np.asarray([np.nan, np.nan])
    result.update({
        "fit_resolved": bool(len(values) >= 450 and low[1] > 0 and delta_aic <= -10),
        "lambda_R": float(damped[0]), "lambda_R_ci": [float(low[0]), float(high[0])],
        "lambda_I": float(damped[1]), "lambda_I_ci": [float(low[1]), float(high[1])],
        "period": float(2*np.pi/damped[1]),
        "delta_AIC_damped_minus_exponential": delta_aic,
        "bootstrap_accept": len(values),
        "independent_lambda_I": INDEPENDENT_IMAG[case],
        "independent_difference": float(damped[1]-INDEPENDENT_IMAG[case]),
    })
    prediction = np.exp(damped[0]*t)*(
        damped[2]*np.cos(damped[1]*t)+damped[3]*np.sin(damped[1]*t)
    )
    return result, mean, se, prediction


def intervals_overlap(left, right):
    return max(left[0], right[0]) <= min(left[1], right[1])


def main() -> None:
    selection = json.loads(SELECTION.read_text())
    selection_sha256 = hashlib.sha256(SELECTION.read_bytes()).hexdigest()
    marker = OUT/"BLIND_TEST_READ_MARKER.json"
    if marker.exists():
        raise RuntimeError("blind test has already been marked as read")
    # Complete all path/model preflight checks before opening any blind data.
    for key, chosen in sorted(selection["selections"].items()):
        case, _ = key.split(":")
        fit_path = ROOT/chosen["frozen_fit"]
        if not fit_path.is_file():
            raise RuntimeError(f"missing frozen fit: {fit_path}")
        paths = split_paths(case, "test")
        if len(paths) != 16 or any(not path.is_file() for path in paths):
            raise RuntimeError(f"bad blind-test path inventory for {case}")
    marker.write_text(json.dumps({
        "selection_file": str(SELECTION.relative_to(ROOT)),
        "selection_sha256": selection_sha256,
        "test_streams_read_once": True,
        "status": "started",
    }, indent=2, sort_keys=True)+"\n")
    summaries = {}
    correlation_rows = []
    residual_rows = []
    for number, (key, chosen) in enumerate(sorted(selection["selections"].items())):
        case, kind = key.split(":")
        fit_data = np.load(ROOT/chosen["frozen_fit"])
        coefficient = phase_align(fit_data["coefficient"])
        indices = fit_data["dictionary_indices"]
        mu = complex(*fit_data["mu"])
        lag_steps = int(round(chosen["tau"]/SNAPSHOT))
        q_values = []
        numerators, denominators = [], []
        paths = split_paths(case, "test")
        if len(paths) != 16:
            raise RuntimeError(f"bad blind-test stream count for {case}")
        stream_ids = []
        for stream_number, path in enumerate(paths, start=1):
            stream_id = int(path.stem.rsplit("_", 1)[-1])
            stream_ids.append(stream_id)
            q = evaluate_q(path, coefficient, indices,
                           fit_data["log_mean"], fit_data["log_sd"])
            numerator, denominator = test_residual(q[None, :], mu, lag_steps)
            numerators.append(numerator); denominators.append(denominator)
            q_values.append(q.real)
            print(f"blind {key}: stream {stream_number}/16", flush=True)
        signals = np.asarray(q_values)
        corr, grid, cgrid = correlation_by_stream(signals)
        rng = np.random.default_rng(BOOTSTRAP_SEED+number)
        numerators = np.asarray(numerators); denominators = np.asarray(denominators)
        residual = math.sqrt(numerators.sum()/denominators.sum())
        residual_boot = bootstrap_residual(numerators, denominators, rng)
        if kind == "real":
            result, mean, se, curve = fit_real_mode(case, cgrid, grid, rng)
            times = grid*SNAPSHOT
        else:
            result, mean, se, curve = fit_complex_mode(case, cgrid, grid, rng)
            times = grid[(grid*SNAPSHOT >= FIT_START-1e-12)&
                         (grid*SNAPSHOT <= FIT_END+1e-12)]*SNAPSHOT
        result.update({
            "selection_local_admissible": chosen["local_admissible"],
            "selection_timestep_pair_admissible": chosen["timestep_pair_admissible"],
            "validation_residual": chosen["validation_residual"],
            "test_residual": residual,
            "test_residual_ci": np.percentile(residual_boot, (2.5, 97.5)).tolist(),
            "residual_growth_ratio": residual/chosen["validation_residual"],
            "residual_gate_pass": bool(residual <= 1.25*chosen["validation_residual"]),
            "test_streams": 16,
        })
        summaries[key] = result
        for stream_order, (stream_id, num, den) in enumerate(
                zip(stream_ids, numerators, denominators), start=1):
            residual_rows.append({"selection": key,
                                  "stream_order": stream_order,
                                  "stream_id": stream_id,
                                  "squared_residual_sum": num, "q_squared_sum": den})
        for t, value, error, fitted in zip(times, mean, se, curve):
            correlation_rows.append({"selection": key, "time": t,
                                     "C": value, "SE_C": error,
                                     "fitted_or_lambda_eff": fitted})
        np.savez_compressed(
            OUT/f"{case}_{kind}_blind_correlations.npz",
            lags=np.arange(legacy.MAX_LAG+1),
            times=np.arange(legacy.MAX_LAG+1)*SNAPSHOT,
            C_by_stream=corr,
        )

    for bath in ("driven", "equilibrium"):
        for kind in ("real", "complex"):
            coarse = summaries[f"{bath}_dt1e-3:{kind}"]
            fine = summaries[f"{bath}_dt2p5e-4:{kind}"]
            if kind == "real":
                overlap = (
                    coarse.get("plateau_resolved", False)
                    and fine.get("plateau_resolved", False)
                    and intervals_overlap(coarse["lambda_R_ci"], fine["lambda_R_ci"])
                )
                for result in (coarse, fine):
                    result["timestep_blind_ci_overlap"] = bool(overlap)
                    result["reportable_edmd_visible_mode"] = bool(
                        result["selection_local_admissible"]
                        and result["selection_timestep_pair_admissible"]
                        and result["residual_gate_pass"]
                        and result.get("plateau_resolved", False) and overlap
                    )
            else:
                overlap = (
                    coarse.get("fit_resolved", False) and fine.get("fit_resolved", False)
                    and intervals_overlap(coarse["lambda_R_ci"], fine["lambda_R_ci"])
                    and intervals_overlap(coarse["lambda_I_ci"], fine["lambda_I_ci"])
                )
                for result in (coarse, fine):
                    result["timestep_blind_ci_overlap"] = bool(overlap)
                    result["reportable_oscillatory_mode"] = bool(
                        result["selection_local_admissible"]
                        and result["selection_timestep_pair_admissible"]
                        and result["residual_gate_pass"]
                        and result.get("fit_resolved", False) and overlap
                    )

    with (OUT/"blind_mode_residual_components.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(residual_rows[0]))
        writer.writeheader(); writer.writerows(residual_rows)
    with (OUT/"blind_mode_correlations.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(correlation_rows[0]))
        writer.writeheader(); writer.writerows(correlation_rows)
    result = {
        "protocol_version": "section4-v1", "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_seed": BOOTSTRAP_SEED, "summaries": summaries,
        "unique_spectral_gap_claim_allowed": False,
        "unique_gap_blocker": (
            "the pre-existing common controlled-autocorrelation plateau gate fails"
        ),
    }
    (OUT/"BLIND_MODE_RESULTS.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    marker.write_text(json.dumps({
        "selection_file": str(SELECTION.relative_to(ROOT)),
        "selection_sha256": selection_sha256,
        "test_streams_read_once": True,
        "status": "complete",
        "blind_results": "modes/BLIND_MODE_RESULTS.json",
    }, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
