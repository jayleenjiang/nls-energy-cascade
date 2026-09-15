#!/usr/bin/env python3
"""Targeted EDMD dictionary extension on existing controlled n=3 trajectories.

Scientific choices are frozen in PROTOCOL.md.  This script is analysis-only:
it reads the saved streams and never advances the stochastic dynamics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy import linalg


HERE = Path(__file__).resolve().parent
BASE_SCRIPT = HERE.parent / "edmd_three_mode_2026-09-07" / "run_edmd.py"
SPEC = importlib.util.spec_from_file_location("base_edmd", BASE_SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {BASE_SCRIPT}")
base = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = base
SPEC.loader.exec_module(base)

NROWS = base.NROWS
NCOLS = base.NCOLS
SNAPSHOT = base.SNAPSHOT
N_ORIGINS = base.N_ORIGINS
CASES = base.CASES
CUTOFF = 1.0e-10
CUTOFF_SENSITIVITY = (1.0e-8, 1.0e-10, 1.0e-12)
LAGS = (0.02, 0.05, 0.10, 0.25, 0.50, 1.00)
LAG_STEPS = tuple(int(round(t / SNAPSHOT)) for t in LAGS)
DICTIONARIES = ("E1", "E2", "E3")
BOOTSTRAPS = 500
BOOTSTRAP_SEED = 2026090702
OBSERVABLES = ("I2", "I1_plus_I3", "cos_theta3", "cos_even_sum")

EXTRA_LABELS = (
    "I1", "I2", "I3",
    "I1_sq", "I1_I2", "I1_I3", "I2_sq", "I2_I3", "I3_sq",
    "dI", "dI_cos_theta1", "dI_sin_theta1", "dI_cos_theta3",
    "dI_sin_theta3", "dI_cos_theta1_minus_theta3",
    "dI_sin_theta1_minus_theta3", "dI_sq", "dI_sumI", "dI_I2",
)
EXTRA_INDICES = np.arange(500, 519, dtype=np.int64)
DICT_INDICES = {
    "E1": np.concatenate((base.DICT_INDICES["D1"], EXTRA_INDICES)),
    "E2": np.concatenate((base.DICT_INDICES["D2"], EXTRA_INDICES)),
    "E3": np.concatenate((base.DICT_INDICES["D3"], EXTRA_INDICES)),
}
assert [len(DICT_INDICES[x]) for x in DICTIONARIES] == [109, 269, 519]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def extended_features(raw: np.ndarray, log_mean: np.ndarray, log_sd: np.ndarray) -> np.ndarray:
    old = base.features(raw, log_mean, log_sd)
    i1, i2, i3, theta1, theta3 = base.reconstruct(raw)
    d = i1 - i3
    extras = np.column_stack((
        i1, i2, i3,
        i1 * i1, i1 * i2, i1 * i3, i2 * i2, i2 * i3, i3 * i3,
        d,
        d * np.cos(theta1), d * np.sin(theta1),
        d * np.cos(theta3), d * np.sin(theta3),
        d * np.cos(theta1 - theta3), d * np.sin(theta1 - theta3),
        d * d, d * (i1 + i3), d * i2,
    ))
    return np.column_stack((old, extras))


def observable_values(raw: np.ndarray) -> np.ndarray:
    return np.column_stack((raw[:, 5], raw[:, 6], raw[:, 2], raw[:, 8])).astype(np.float64)


def matrix_pass(case_dir: Path, log_mean: np.ndarray, log_sd: np.ndarray):
    files = sorted(case_dir.glob("stream_*.f32"))
    if len(files) != 64:
        raise RuntimeError(f"{case_dir}: expected 64 streams, found {len(files)}")
    max_lag = max(LAG_STEPS)
    origins = np.rint(np.linspace(0, NROWS - 1 - max_lag, N_ORIGINS)).astype(np.int64)
    if len(np.unique(origins)) != N_ORIGINS:
        raise RuntimeError("origin quadrature contains duplicates")

    K = 519
    As = np.empty((64, K, K), dtype=np.float32)
    Bs = np.empty((len(LAGS), 64, K, K), dtype=np.float32)
    bs = np.empty((64, K, len(OBSERVABLES)), dtype=np.float32)
    gsums = np.empty((64, len(OBSERVABLES)), dtype=np.float64)
    for s, path in enumerate(files):
        raw = np.memmap(path, dtype="<f4", mode="r", shape=(NROWS, NCOLS))
        xraw = np.asarray(raw[origins])
        X = extended_features(xraw, log_mean, log_sd)
        g = observable_values(xraw)
        As[s] = (X.T @ X / N_ORIGINS).astype(np.float32)
        bs[s] = (X.T @ g / N_ORIGINS).astype(np.float32)
        gsums[s] = g.mean(axis=0)
        for ell, step in enumerate(LAG_STEPS):
            yraw = np.asarray(raw[origins + step])
            Y = extended_features(yraw, log_mean, log_sd)
            Bs[ell, s] = (X.T @ Y / N_ORIGINS).astype(np.float32)
        print(f"  matrices {case_dir.name}: stream {s + 1}/64", flush=True)
    return As, Bs, bs, gsums, origins


def subset_matrix(M: np.ndarray, idx: np.ndarray):
    return M[np.ix_(idx, idx)].astype(np.float64, copy=False)


def spectrum_rows(case, dictionary, tau, cutoff, fit):
    rows = []
    for rank, (j, mu, lam) in enumerate(base.canonical_modes(fit), start=0):
        rows.append((case, dictionary, tau, math.pi / tau, cutoff, rank, j,
                     mu.real, mu.imag, abs(mu), lam.real, lam.imag,
                     int(abs(lam.imag) >= 0.8 * math.pi / tau),
                     fit.kept, fit.condition, fit.trivial_error))
    return rows


def select_oscillatory(fit):
    candidates = []
    for j, mu, lam in base.canonical_modes(fit):
        if 4.5 <= lam.imag <= 6.0 and -20.0 <= lam.real <= -0.05:
            candidates.append((j, mu, lam))
    if not candidates:
        return None
    return max(candidates, key=lambda x: x[2].real)


def match(reference: complex, fit):
    candidates = fit.lam.copy()
    candidates[candidates.imag < -1.0e-8] = np.nan + 1j * np.nan
    return base.match_rate(reference, candidates)


def bootstrap_oscillatory(case, As, Bs, reference, output: Path):
    rng = np.random.default_rng(BOOTSTRAP_SEED + CASES.index(case))
    rows = []
    Bstream = Bs[LAGS.index(0.05)]
    for b in range(BOOTSTRAPS):
        picked = rng.integers(0, 64, size=64)
        counts = np.bincount(picked, minlength=64).astype(np.float64) / 64.0
        A = np.tensordot(counts, As, axes=(0, 0))
        B = np.tensordot(counts, Bstream, axes=(0, 0))
        try:
            fit = base.solve_edmd(A, B, 0.05, CUTOFF, vectors=False)
            j, distance = match(reference, fit)
            if j is None:
                rows.append((case, b, reference.real, reference.imag, 0, distance,
                             math.nan, math.nan, math.nan, math.nan,
                             fit.kept, fit.trivial_error))
            else:
                rows.append((case, b, reference.real, reference.imag, 1, distance,
                             fit.mu[j].real, fit.mu[j].imag,
                             fit.lam[j].real, fit.lam[j].imag,
                             fit.kept, fit.trivial_error))
        except Exception:
            rows.append((case, b, reference.real, reference.imag, 0, math.nan,
                         math.nan, math.nan, math.nan, math.nan, 0, math.nan))
        if (b + 1) % 25 == 0:
            print(f"  oscillatory bootstrap {case}: {b + 1}/{BOOTSTRAPS}", flush=True)
    write_csv(output / f"{case}_bootstrap_oscillatory_raw.csv",
              ("case", "replicate", "reference_real", "reference_imag", "matched",
               "distance", "mu_real", "mu_imag", "lambda_real", "lambda_imag",
               "gram_kept", "trivial_error"), rows)


def bootstrap_real(case, As, Bs, reference_fit, output: Path):
    refs = base.reference_modes(reference_fit)
    rng = np.random.default_rng(BOOTSTRAP_SEED + 100 + CASES.index(case))
    rows = []
    Bstream = Bs[LAGS.index(0.50)]
    for b in range(BOOTSTRAPS):
        picked = rng.integers(0, 64, size=64)
        counts = np.bincount(picked, minlength=64).astype(np.float64) / 64.0
        A = np.tensordot(counts, As, axes=(0, 0))
        B = np.tensordot(counts, Bstream, axes=(0, 0))
        try:
            fit = base.solve_edmd(A, B, 0.50, CUTOFF, vectors=False)
            for r, (_, _, lam_ref) in enumerate(refs):
                j, distance = match(lam_ref, fit)
                if j is None:
                    rows.append((case, b, r, lam_ref.real, lam_ref.imag, 0,
                                 distance, math.nan, math.nan, math.nan, math.nan,
                                 fit.kept, fit.trivial_error))
                else:
                    rows.append((case, b, r, lam_ref.real, lam_ref.imag, 1,
                                 distance, fit.mu[j].real, fit.mu[j].imag,
                                 fit.lam[j].real, fit.lam[j].imag,
                                 fit.kept, fit.trivial_error))
        except Exception:
            for r, (_, _, lam_ref) in enumerate(refs):
                rows.append((case, b, r, lam_ref.real, lam_ref.imag, 0, math.nan,
                             math.nan, math.nan, math.nan, math.nan, 0, math.nan))
        if (b + 1) % 25 == 0:
            print(f"  real bootstrap {case}: {b + 1}/{BOOTSTRAPS}", flush=True)
    write_csv(output / f"{case}_bootstrap_real_raw.csv",
              ("case", "replicate", "reference_mode", "reference_real",
               "reference_imag", "matched", "distance", "mu_real", "mu_imag",
               "lambda_real", "lambda_imag", "gram_kept", "trivial_error"), rows)
    return refs


def modal_weights(case, A, b, gmean, fit, refs, output: Path):
    constant = 0
    Bc = b - A[:, [constant]] * gmean[None, :]
    d = fit.transform.T @ Bc
    V = fit.reduced_vectors
    Vinv = linalg.inv(V, check_finite=False)
    amplitudes = Vinv @ d
    inner = Bc.T @ (fit.transform @ V)
    raw_contrib = inner.T * amplitudes
    reconstructed_var = np.sum(raw_contrib, axis=0)
    rows = []
    for r, (j, _, lam) in enumerate(refs):
        if lam.imag > 1.0e-8:
            partner = int(np.argmin(np.abs(fit.lam - np.conj(lam))))
            contribution = raw_contrib[j] + raw_contrib[partner]
        else:
            partner = -1
            contribution = raw_contrib[j]
        for q, obs in enumerate(OBSERVABLES):
            denom = reconstructed_var[q]
            weight = contribution[q] / denom if abs(denom) > 0 else np.nan
            rows.append((case, r, j, partner, lam.real, lam.imag, obs,
                         contribution[q].real, contribution[q].imag,
                         denom.real, denom.imag, weight.real, weight.imag, abs(weight)))
    write_csv(output / f"{case}_observable_modal_weights.csv",
              ("case", "reference_mode", "eigen_index", "conjugate_index",
               "lambda_real", "lambda_imag", "observable", "contribution_real",
               "contribution_imag", "reconstructed_variance_real",
               "reconstructed_variance_imag", "weight_real", "weight_imag",
               "weight_abs"), rows)


def save_slices(case, fit_osc, osc, fit_real, real, log_mean, log_sd, output: Path):
    grid = np.linspace(-math.pi, math.pi, 81)
    rows = []
    modes = (("oscillatory", fit_osc, osc), ("leading_real", fit_real, real))
    for label, fit, selected in modes:
        if selected is None:
            continue
        j, _, lam = selected
        coeff = fit.coeff[:, j].copy()
        pivot = int(np.argmax(np.abs(coeff)))
        if abs(coeff[pivot]) > 0:
            coeff *= np.exp(-1j * np.angle(coeff[pivot]))
        for action in (1.0, 2.0, 4.0):
            for th1 in grid:
                raw = np.zeros((len(grid), NCOLS), dtype=np.float64)
                raw[:, 0] = np.cos(th1)
                raw[:, 1] = np.sin(th1)
                raw[:, 2] = np.cos(grid)
                raw[:, 3] = np.sin(grid)
                raw[:, 5] = action
                raw[:, 6] = 2.0 * action
                raw[:, 7] = 0.0
                raw[:, 8] = raw[:, 0] + raw[:, 2]
                phi = extended_features(raw, log_mean, log_sd) @ coeff
                for th3, value in zip(grid, phi):
                    rows.append((case, label, lam.real, lam.imag, action,
                                 th1, th3, value.real, value.imag))
    write_csv(output / f"{case}_eigenfunction_slices.csv",
              ("case", "mode", "lambda_real", "lambda_imag", "I_equal",
               "theta1", "theta3", "phi_real", "phi_imag"), rows)


def run_case(case: str, input_root: Path, output: Path):
    started = time.time()
    case_dir = input_root / case
    files = sorted(case_dir.glob("stream_*.f32"))
    if len(files) != 64:
        raise RuntimeError(f"{case}: found {len(files)} stream files")
    expected = NROWS * NCOLS * 4
    for p in files:
        if p.stat().st_size != expected:
            raise RuntimeError(f"bad file size {p}: {p.stat().st_size} != {expected}")

    mean, sd, mins, maxs, count = base.full_log_stats(files)
    write_csv(output / f"{case}_standardization.csv",
              ("case", "variable", "log_mean", "log_sd", "min_action",
               "max_action", "snapshot_count"),
              [(case, f"I{j+1}", mean[j], sd[j], mins[j], maxs[j], count)
               for j in range(3)])
    As, Bs, bs, gsums, origins = matrix_pass(case_dir, mean, sd)
    Amax = As.mean(axis=0, dtype=np.float64)
    Bmax = Bs.mean(axis=1, dtype=np.float64)
    bmax = bs.mean(axis=0, dtype=np.float64)
    gmean = gsums.mean(axis=0)
    np.savez_compressed(output / f"{case}_aggregate_matrices.npz",
                        A=Amax, B=Bmax, b=bmax, gmean=gmean, origins=origins,
                        lags=np.asarray(LAGS))

    spectrum, diagnostics, fits = [], [], {}
    for label in DICTIONARIES:
        idx = DICT_INDICES[label]
        A = subset_matrix(Amax, idx)
        for ell, tau in enumerate(LAGS):
            B = subset_matrix(Bmax[ell], idx)
            fit = base.solve_edmd(A, B, tau, CUTOFF, vectors=True)
            fits[(label, tau, CUTOFF)] = fit
            spectrum.extend(spectrum_rows(case, label, tau, CUTOFF, fit))
            diagnostics.append((case, label, len(idx), tau, math.pi / tau,
                                CUTOFF, fit.kept, fit.condition,
                                fit.trivial_error, fit.mu[fit.trivial_index].real,
                                fit.mu[fit.trivial_index].imag))
    for tau in (0.05, 0.50):
        for cutoff in CUTOFF_SENSITIVITY:
            if cutoff == CUTOFF:
                continue
            ell = LAGS.index(tau)
            fit = base.solve_edmd(Amax, Bmax[ell], tau, cutoff, vectors=True)
            fits[("E3", tau, cutoff)] = fit
            spectrum.extend(spectrum_rows(case, "E3", tau, cutoff, fit))
            diagnostics.append((case, "E3", 519, tau, math.pi / tau, cutoff,
                                fit.kept, fit.condition, fit.trivial_error,
                                fit.mu[fit.trivial_index].real,
                                fit.mu[fit.trivial_index].imag))

    write_csv(output / f"{case}_spectra.csv",
              ("case", "dictionary", "tau", "nyquist", "cutoff", "rank",
               "eigen_index", "mu_real", "mu_imag", "mu_abs", "lambda_real",
               "lambda_imag", "aliased_20pct", "gram_kept", "gram_condition",
               "trivial_error"), spectrum)
    write_csv(output / f"{case}_gram_diagnostics.csv",
              ("case", "dictionary", "nominal_K", "tau", "nyquist", "cutoff",
               "kept", "condition", "trivial_error", "trivial_mu_real",
               "trivial_mu_imag"), diagnostics)

    fit_osc = fits[("E3", 0.05, CUTOFF)]
    osc = select_oscillatory(fit_osc)
    if osc is None:
        osc_ref = complex(math.nan, math.nan)
    else:
        osc_ref = osc[2]
    bootstrap_oscillatory(case, As, Bs, osc_ref, output)

    fit_real = fits[("E3", 0.50, CUTOFF)]
    refs = bootstrap_real(case, As, Bs, fit_real, output)
    modal_weights(case, Amax, bmax, gmean, fit_real, refs, output)
    real = refs[0] if refs else None
    save_slices(case, fit_osc, osc, fit_real, real, mean, sd, output)

    del As, Bs, bs
    return {
        "case": case,
        "elapsed_seconds": time.time() - started,
        "input_streams": 64,
        "snapshot_count": count,
        "origin_count_per_stream": N_ORIGINS,
        "pair_count": N_ORIGINS * 64,
        "oscillatory_reference": None if osc is None else {
            "lambda_real": float(osc[2].real), "lambda_imag": float(osc[2].imag)},
        "real_references": [
            {"mode": r, "lambda_real": float(lam.real),
             "lambda_imag": float(lam.imag)}
            for r, (_, _, lam) in enumerate(refs)
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--case", choices=CASES, action="append")
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    selected = tuple(args.case) if args.case else CASES
    started = time.time()
    manifest = {
        "started_unix": started,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "input": str(args.input.resolve()),
        "output": str(args.output.resolve()),
        "cases": list(selected),
        "lags": list(LAGS),
        "nyquist": {str(t): math.pi / t for t in LAGS},
        "origins_per_stream": N_ORIGINS,
        "bootstraps": BOOTSTRAPS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "cutoff": CUTOFF,
        "cutoff_sensitivity": list(CUTOFF_SENSITIVITY),
        "dictionaries": [{"label": x, "nominal_K": len(DICT_INDICES[x])}
                         for x in DICTIONARIES],
        "extra_features": list(EXTRA_LABELS),
        "base_analysis_source": str(BASE_SCRIPT),
        "base_analysis_sha256": sha256(BASE_SCRIPT),
        "results": [],
    }
    for case in selected:
        manifest["results"].append(run_case(case, args.input, args.output))
    manifest["elapsed_seconds"] = time.time() - started
    manifest["analysis_source_sha256"] = sha256(Path(__file__).resolve())
    with (args.output / "run_manifest.json").open("w") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps({"status": "complete",
                      "elapsed_seconds": manifest["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
