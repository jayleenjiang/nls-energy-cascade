#!/usr/bin/env python3
"""Controlled EDMD analysis for the saved n=3 trajectories.

All scientific choices are frozen in PROTOCOL.md.  This script does not alter
or regenerate trajectory data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import linalg


NROWS = 100001
NCOLS = 12
SNAPSHOT = 0.01
N_ORIGINS = 2048
LAGS = (0.10, 0.25, 0.50, 1.00)
LAG_STEPS = tuple(int(round(t / SNAPSHOT)) for t in LAGS)
CUTOFF = 1.0e-10
CUTOFF_SENSITIVITY = (1.0e-8, 1.0e-10, 1.0e-12)
DICTIONARIES = (('D1', 2, 1), ('D2', 2, 2), ('D3', 3, 2))
BOOTSTRAPS = 500
BOOTSTRAP_SEED = 2026090701
MAX_REFERENCE_MODES = 12
OBSERVABLES = ('I2', 'I1_plus_I3', 'cos_theta3', 'cos_even_sum', 'I1_minus_I3')
CASES = (
    'driven_dt1e-3',
    'driven_dt2p5e-4',
    'equilibrium_dt1e-3',
    'equilibrium_dt2p5e-4',
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def action_exponents(dmax: int = 3):
    out = []
    for degree in range(dmax + 1):
        for a in range(degree + 1):
            for b in range(degree - a + 1):
                c = degree - a - b
                out.append((a, b, c))
    return out


def angle_terms(mmax: int = 2):
    out = [('const', 0, 0)]
    for m1 in range(-mmax, mmax + 1):
        for m3 in range(-mmax, mmax + 1):
            if m1 == 0 and m3 == 0:
                continue
            # One representative of each m/-m pair.
            if not (m1 > 0 or (m1 == 0 and m3 > 0)):
                continue
            out.append(('cos', m1, m3))
            out.append(('sin', m1, m3))
    return out


ACTION_EXPS = action_exponents(3)
ANGLE_TERMS = angle_terms(2)
assert len(ACTION_EXPS) == 20
assert len(ANGLE_TERMS) == 25


def feature_labels():
    labels = []
    for exp in ACTION_EXPS:
        for kind, m1, m3 in ANGLE_TERMS:
            labels.append((exp, kind, m1, m3))
    return labels


FEATURE_LABELS = feature_labels()
CONSTANT_INDEX = FEATURE_LABELS.index(((0, 0, 0), 'const', 0, 0))


def dictionary_indices(d: int, m: int) -> np.ndarray:
    idx = []
    for i, (exp, kind, m1, m3) in enumerate(FEATURE_LABELS):
        if sum(exp) <= d and abs(m1) <= m and abs(m3) <= m:
            idx.append(i)
    return np.asarray(idx, dtype=np.int64)


DICT_INDICES = {label: dictionary_indices(d, m) for label, d, m in DICTIONARIES}
assert [len(DICT_INDICES[x]) for x in ('D1', 'D2', 'D3')] == [90, 250, 500]


def reconstruct(raw: np.ndarray):
    i2 = raw[:, 5].astype(np.float64)
    isum = raw[:, 6].astype(np.float64)
    idiff = raw[:, 7].astype(np.float64)
    i1 = 0.5 * (isum + idiff)
    i3 = 0.5 * (isum - idiff)
    theta1 = np.arctan2(raw[:, 1].astype(np.float64), raw[:, 0].astype(np.float64))
    theta3 = np.arctan2(raw[:, 3].astype(np.float64), raw[:, 2].astype(np.float64))
    return i1, i2, i3, theta1, theta3


def features(raw: np.ndarray, log_mean: np.ndarray, log_sd: np.ndarray) -> np.ndarray:
    i1, i2, i3, theta1, theta3 = reconstruct(raw)
    actions = np.column_stack((i1, i2, i3))
    if np.any(~np.isfinite(actions)) or np.any(actions < 0):
        raise ValueError('nonfinite or negative reconstructed action')
    z = (np.log(np.maximum(actions, 1.0e-12)) - log_mean) / log_sd

    amon = np.empty((len(raw), len(ACTION_EXPS)), dtype=np.float64)
    for j, (a, b, c) in enumerate(ACTION_EXPS):
        amon[:, j] = (z[:, 0] ** a) * (z[:, 1] ** b) * (z[:, 2] ** c)

    ang = np.empty((len(raw), len(ANGLE_TERMS)), dtype=np.float64)
    for j, (kind, m1, m3) in enumerate(ANGLE_TERMS):
        if kind == 'const':
            ang[:, j] = 1.0
        else:
            phase = m1 * theta1 + m3 * theta3
            ang[:, j] = np.cos(phase) if kind == 'cos' else np.sin(phase)
    return (amon[:, :, None] * ang[:, None, :]).reshape(len(raw), -1)


def observable_values(raw: np.ndarray) -> np.ndarray:
    return np.column_stack((raw[:, 5], raw[:, 6], raw[:, 2], raw[:, 8], raw[:, 7])).astype(np.float64)


@dataclass
class EdmdResult:
    mu: np.ndarray
    lam: np.ndarray
    coeff: np.ndarray
    reduced_vectors: np.ndarray
    transform: np.ndarray
    gram_eigenvalues: np.ndarray
    kept: int
    condition: float
    trivial_index: int
    trivial_error: float


def solve_edmd(A: np.ndarray, B: np.ndarray, tau: float, cutoff: float, vectors: bool = True) -> EdmdResult:
    A = 0.5 * (A + A.T)
    vals, vecs = linalg.eigh(A, check_finite=False, overwrite_a=False)
    vmax = float(vals[-1])
    keep = vals > cutoff * vmax
    if np.count_nonzero(keep) < 2:
        raise RuntimeError('Gram truncation retained fewer than two directions')
    kept_vals = vals[keep]
    U = vecs[:, keep]
    T = U / np.sqrt(kept_vals)[None, :]
    Kwhite = T.T @ B @ T
    if vectors:
        mu, V = linalg.eig(Kwhite, check_finite=False)
        coeff = T @ V
    else:
        mu = linalg.eigvals(Kwhite, check_finite=False)
        V = np.empty((Kwhite.shape[0], 0), dtype=np.complex128)
        coeff = np.empty((A.shape[0], 0), dtype=np.complex128)
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = np.log(mu.astype(np.complex128)) / tau
    trivial = int(np.nanargmin(np.abs(mu - 1.0)))
    return EdmdResult(
        mu=mu,
        lam=lam,
        coeff=coeff,
        reduced_vectors=V,
        transform=T,
        gram_eigenvalues=vals,
        kept=int(np.count_nonzero(keep)),
        condition=float(vmax / kept_vals[0]),
        trivial_index=trivial,
        trivial_error=float(abs(mu[trivial] - 1.0)),
    )


def canonical_modes(result: EdmdResult, limit: int | None = None):
    rows = []
    for j, (mu, lam) in enumerate(zip(result.mu, result.lam)):
        if not (np.isfinite(mu.real) and np.isfinite(mu.imag) and np.isfinite(lam.real) and np.isfinite(lam.imag)):
            continue
        if lam.imag < -1.0e-8:
            continue
        rows.append((j, mu, lam))
    rows.sort(key=lambda x: (-x[2].real, x[2].imag))
    return rows if limit is None else rows[:limit]


def reference_modes(result: EdmdResult):
    out = []
    for j, mu, lam in canonical_modes(result):
        if -3.0 <= lam.real <= -0.05 and -1.0e-8 <= lam.imag <= 12.0 and 0 < abs(mu) <= 1.05:
            out.append((j, mu, lam))
        if len(out) >= MAX_REFERENCE_MODES:
            break
    return out


def match_rate(reference: complex, candidates: np.ndarray):
    valid = np.isfinite(candidates.real) & np.isfinite(candidates.imag) & (candidates.imag >= -1.0e-8)
    idx = np.flatnonzero(valid)
    if idx.size == 0:
        return None, math.inf
    distances = np.abs(candidates[idx] - reference)
    k = int(np.argmin(distances))
    j = int(idx[k])
    distance = float(distances[k])
    threshold = max(0.50, 0.35 * abs(reference))
    return (j if distance <= threshold else None), distance


def write_csv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def save_matrix_csv_gz(path: Path, matrix: np.ndarray):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', newline='') as f:
        np.savetxt(f, matrix, delimiter=',', fmt='%.12e')


def full_log_stats(stream_files):
    count = 0
    s1 = np.zeros(3)
    s2 = np.zeros(3)
    mins = np.full(3, np.inf)
    maxs = np.full(3, -np.inf)
    for path in stream_files:
        raw = np.memmap(path, dtype='<f4', mode='r', shape=(NROWS, NCOLS))
        i1, i2, i3, _, _ = reconstruct(raw)
        actions = np.column_stack((i1, i2, i3))
        logs = np.log(np.maximum(actions, 1.0e-12))
        count += len(logs)
        s1 += logs.sum(axis=0)
        s2 += np.square(logs).sum(axis=0)
        mins = np.minimum(mins, actions.min(axis=0))
        maxs = np.maximum(maxs, actions.max(axis=0))
    mean = s1 / count
    var = (s2 - count * mean * mean) / (count - 1)
    return mean, np.sqrt(var), mins, maxs, count


def matrix_pass(case_dir: Path, log_mean, log_sd):
    files = sorted(case_dir.glob('stream_*.f32'))
    if len(files) != 64:
        raise RuntimeError(f'{case_dir}: expected 64 streams, found {len(files)}')
    max_lag = max(LAG_STEPS)
    origins = np.rint(np.linspace(0, NROWS - 1 - max_lag, N_ORIGINS)).astype(np.int64)
    if len(np.unique(origins)) != N_ORIGINS:
        raise RuntimeError('origin quadrature contains duplicates')

    K = 500
    As = np.empty((64, K, K), dtype=np.float32)
    Bs = np.empty((len(LAGS), 64, K, K), dtype=np.float32)
    bs = np.empty((64, K, len(OBSERVABLES)), dtype=np.float32)
    gsums = np.empty((64, len(OBSERVABLES)), dtype=np.float64)
    for s, path in enumerate(files):
        raw = np.memmap(path, dtype='<f4', mode='r', shape=(NROWS, NCOLS))
        xraw = np.asarray(raw[origins])
        X = features(xraw, log_mean, log_sd)
        g = observable_values(xraw)
        As[s] = (X.T @ X / N_ORIGINS).astype(np.float32)
        bs[s] = (X.T @ g / N_ORIGINS).astype(np.float32)
        gsums[s] = g.mean(axis=0)
        for ell, step in enumerate(LAG_STEPS):
            yraw = np.asarray(raw[origins + step])
            Y = features(yraw, log_mean, log_sd)
            Bs[ell, s] = (X.T @ Y / N_ORIGINS).astype(np.float32)
        print(f'  matrices {case_dir.name}: stream {s + 1}/64', flush=True)
    return As, Bs, bs, gsums, origins


def subset_matrix(M: np.ndarray, idx: np.ndarray):
    return M[np.ix_(idx, idx)].astype(np.float64, copy=False)


def spectrum_rows(case, dictionary, tau, cutoff, result):
    rows = []
    for rank, (j, mu, lam) in enumerate(canonical_modes(result), start=0):
        rows.append((case, dictionary, tau, cutoff, rank, j, mu.real, mu.imag,
                     abs(mu), lam.real, lam.imag, result.kept, result.condition,
                     result.trivial_error))
    return rows


def bootstrap_case(case: str, As: np.ndarray, Bs: np.ndarray, reference: EdmdResult, output: Path):
    refs = reference_modes(reference)
    rng = np.random.default_rng(BOOTSTRAP_SEED + CASES.index(case))
    rows = []
    Bstream = Bs[LAGS.index(0.50)]
    for b in range(BOOTSTRAPS):
        picked = rng.integers(0, 64, size=64)
        counts = np.bincount(picked, minlength=64).astype(np.float64) / 64.0
        A = np.tensordot(counts, As, axes=(0, 0))
        B = np.tensordot(counts, Bstream, axes=(0, 0))
        try:
            fit = solve_edmd(A, B, 0.50, CUTOFF, vectors=False)
            candidates = fit.lam.copy()
            candidates[candidates.imag < 0] = np.nan + 1j * np.nan
            for r, (_, _, lam_ref) in enumerate(refs):
                j, distance = match_rate(lam_ref, candidates)
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
                rows.append((case, b, r, lam_ref.real, lam_ref.imag, 0,
                             math.nan, math.nan, math.nan, math.nan, math.nan,
                             0, math.nan))
        if (b + 1) % 25 == 0:
            print(f'  bootstrap {case}: {b + 1}/{BOOTSTRAPS}', flush=True)
    write_csv(output / f'{case}_bootstrap_raw.csv',
              ('case', 'replicate', 'reference_mode', 'reference_real', 'reference_imag',
               'matched', 'distance', 'mu_real', 'mu_imag', 'lambda_real',
               'lambda_imag', 'gram_kept', 'trivial_error'), rows)
    return refs


def modal_weights(case, A, b, gmean, fit: EdmdResult, refs, output: Path):
    # Project centered observables onto the retained dictionary.
    econst = np.zeros(A.shape[0])
    econst[CONSTANT_INDEX] = 1.0
    Bc = b - A[:, [CONSTANT_INDEX]] * gmean[None, :]
    d = fit.transform.T @ Bc
    V = fit.reduced_vectors
    Vinv = linalg.inv(V, check_finite=False)
    amplitudes = Vinv @ d
    inner = Bc.T @ (fit.transform @ V)
    raw_contrib = inner.T * amplitudes
    reconstructed_var = np.sum(raw_contrib, axis=0)

    rows = []
    for r, (j, mu, lam) in enumerate(refs):
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
                         denom.real, denom.imag, weight.real, weight.imag,
                         abs(weight)))
    write_csv(output / f'{case}_observable_modal_weights.csv',
              ('case', 'reference_mode', 'eigen_index', 'conjugate_index',
               'lambda_real', 'lambda_imag', 'observable', 'contribution_real',
               'contribution_imag', 'reconstructed_variance_real',
               'reconstructed_variance_imag', 'weight_real', 'weight_imag',
               'weight_abs'), rows)


def save_slices(case, fit: EdmdResult, refs, log_mean, log_sd, output: Path):
    # Save the first four diagnostic/reference modes.  Validation labels are
    # applied later after all convergence checks are available.
    grid = np.linspace(-math.pi, math.pi, 81)
    rows = []
    for r, (j, _, lam) in enumerate(refs[:4]):
        coeff = fit.coeff[:, j].copy()
        pivot = int(np.argmax(np.abs(coeff)))
        if abs(coeff[pivot]) > 0:
            coeff *= np.exp(-1j * np.angle(coeff[pivot]))
        for action in (1.0, 2.0, 4.0):
            for th1 in grid:
                raw = np.empty((len(grid), NCOLS), dtype=np.float64)
                raw.fill(0.0)
                raw[:, 0] = np.cos(th1)
                raw[:, 1] = np.sin(th1)
                raw[:, 2] = np.cos(grid)
                raw[:, 3] = np.sin(grid)
                raw[:, 5] = action
                raw[:, 6] = 2.0 * action
                raw[:, 7] = 0.0
                raw[:, 8] = raw[:, 0] + raw[:, 2]
                Phi = features(raw, log_mean, log_sd)
                value = Phi @ coeff
                for th3, val in zip(grid, value):
                    rows.append((case, r, lam.real, lam.imag, action, th1, th3,
                                 val.real, val.imag))
    write_csv(output / f'{case}_eigenfunction_slices.csv',
              ('case', 'reference_mode', 'lambda_real', 'lambda_imag', 'I_equal',
               'theta1', 'theta3', 'phi_real', 'phi_imag'), rows)


def run_case(case: str, input_root: Path, output: Path):
    t0 = time.time()
    case_dir = input_root / case
    files = sorted(case_dir.glob('stream_*.f32'))
    for p in files:
        expected = NROWS * NCOLS * 4
        if p.stat().st_size != expected:
            raise RuntimeError(f'bad file size {p}: {p.stat().st_size} != {expected}')
    mean, sd, mins, maxs, count = full_log_stats(files)
    write_csv(output / f'{case}_standardization.csv',
              ('case', 'variable', 'log_mean', 'log_sd', 'min_action', 'max_action', 'snapshot_count'),
              [(case, f'I{j+1}', mean[j], sd[j], mins[j], maxs[j], count) for j in range(3)])
    print(f'{case}: log-action mean={mean}, sd={sd}', flush=True)
    As, Bs, bs, gsums, origins = matrix_pass(case_dir, mean, sd)
    Amax = As.mean(axis=0, dtype=np.float64)
    Bmax = Bs.mean(axis=1, dtype=np.float64)
    bmax = bs.mean(axis=0, dtype=np.float64)
    gmean = gsums.mean(axis=0)

    matrix_dir = output / 'matrices'
    save_matrix_csv_gz(matrix_dir / f'{case}_D3_A.csv.gz', Amax)
    for ell, tau in enumerate(LAGS):
        save_matrix_csv_gz(matrix_dir / f'{case}_D3_B_tau{tau:.2f}.csv.gz', Bmax[ell])

    spectrum = []
    fits = {}
    diagnostics = []
    for label, d, m in DICTIONARIES:
        idx = DICT_INDICES[label]
        A = subset_matrix(Amax, idx)
        for ell, tau in enumerate(LAGS):
            B = subset_matrix(Bmax[ell], idx)
            fit = solve_edmd(A, B, tau, CUTOFF, vectors=True)
            fits[(label, tau, CUTOFF)] = fit
            spectrum.extend(spectrum_rows(case, label, tau, CUTOFF, fit))
            diagnostics.append((case, label, d, m, len(idx), tau, CUTOFF,
                                fit.kept, fit.condition, fit.trivial_error,
                                fit.mu[fit.trivial_index].real,
                                fit.mu[fit.trivial_index].imag))
    # Threshold sensitivity is evaluated at the primary D3, tau=0.5 fit.
    idx = DICT_INDICES['D3']
    for cutoff in CUTOFF_SENSITIVITY:
        if cutoff == CUTOFF:
            continue
        fit = solve_edmd(Amax, Bmax[LAGS.index(0.50)], 0.50, cutoff, vectors=True)
        fits[('D3', 0.50, cutoff)] = fit
        spectrum.extend(spectrum_rows(case, 'D3', 0.50, cutoff, fit))
        diagnostics.append((case, 'D3', 3, 2, 500, 0.50, cutoff, fit.kept,
                            fit.condition, fit.trivial_error,
                            fit.mu[fit.trivial_index].real,
                            fit.mu[fit.trivial_index].imag))

    write_csv(output / f'{case}_spectra.csv',
              ('case', 'dictionary', 'tau', 'cutoff', 'rank', 'eigen_index',
               'mu_real', 'mu_imag', 'mu_abs', 'lambda_real', 'lambda_imag',
               'gram_kept', 'gram_condition', 'trivial_error'), spectrum)
    write_csv(output / f'{case}_gram_diagnostics.csv',
              ('case', 'dictionary', 'd', 'M', 'K', 'tau', 'cutoff', 'kept',
               'condition', 'trivial_error', 'trivial_mu_real', 'trivial_mu_imag'),
              diagnostics)

    primary = fits[('D3', 0.50, CUTOFF)]
    refs = bootstrap_case(case, As, Bs, primary, output)
    modal_weights(case, Amax, bmax, gmean, primary, refs, output)
    save_slices(case, primary, refs, mean, sd, output)

    del As, Bs, bs
    return {
        'case': case,
        'elapsed_seconds': time.time() - t0,
        'log_mean': mean.tolist(),
        'log_sd': sd.tolist(),
        'action_min': mins.tolist(),
        'action_max': maxs.tolist(),
        'snapshot_count': count,
        'origin_count_per_stream': N_ORIGINS,
        'pair_count': N_ORIGINS * 64,
        'reference_modes': [
            {'mode': r, 'lambda_real': float(lam.real), 'lambda_imag': float(lam.imag),
             'mu_real': float(mu.real), 'mu_imag': float(mu.imag)}
            for r, (_, mu, lam) in enumerate(refs)
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--case', choices=CASES, action='append')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    selected = tuple(args.case) if args.case else CASES
    started = time.time()
    manifest = {
        'started_unix': started,
        'python': sys.version,
        'platform': platform.platform(),
        'numpy': np.__version__,
        'scipy': __import__('scipy').__version__,
        'input': str(args.input.resolve()),
        'output': str(args.output.resolve()),
        'cases': list(selected),
        'lags': list(LAGS),
        'origins_per_stream': N_ORIGINS,
        'bootstraps': BOOTSTRAPS,
        'bootstrap_seed': BOOTSTRAP_SEED,
        'cutoff': CUTOFF,
        'cutoff_sensitivity': list(CUTOFF_SENSITIVITY),
        'dictionaries': [{'label': a, 'd': b, 'M': c, 'K': len(DICT_INDICES[a])}
                         for a, b, c in DICTIONARIES],
        'feature_labels': [
            {'index': i, 'action_exp': list(exp), 'kind': kind, 'm1': m1, 'm3': m3}
            for i, (exp, kind, m1, m3) in enumerate(FEATURE_LABELS)
        ],
        'results': [],
    }
    for case in selected:
        manifest['results'].append(run_case(case, args.input, args.output))
    manifest['elapsed_seconds'] = time.time() - started
    script = Path(__file__).resolve()
    manifest['analysis_source_sha256'] = sha256(script)
    with (args.output / 'run_manifest.json').open('w') as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps({'status': 'complete', 'elapsed_seconds': manifest['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()

