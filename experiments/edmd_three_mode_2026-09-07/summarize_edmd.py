#!/usr/bin/env python3
"""Apply the frozen convergence gates to run_edmd.py outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


CASES = (
    'driven_dt1e-3', 'driven_dt2p5e-4',
    'equilibrium_dt1e-3', 'equilibrium_dt2p5e-4')
LAGS = (0.10, 0.25, 0.50, 1.00)
DICTS = ('D1', 'D2', 'D3')
CUTOFF = 1e-10


def read_csv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def write_csv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def cnum(row, prefix='lambda'):
    return complex(float(row[prefix + '_real']), float(row[prefix + '_imag']))


def select(rows, dictionary, tau, cutoff=CUTOFF):
    return [r for r in rows
            if r['dictionary'] == dictionary
            and abs(float(r['tau']) - tau) < 1e-12
            and abs(float(r['cutoff']) - cutoff) <= 1e-18]


def reference(rows):
    out = []
    for r in select(rows, 'D3', 0.5):
        lam, mu = cnum(r), cnum(r, 'mu')
        if -3 <= lam.real <= -0.05 and -1e-8 <= lam.imag <= 12 and 0 < abs(mu) <= 1.05:
            out.append(r)
        if len(out) == 12:
            break
    return out


def global_match(ref_rows, candidate_rows):
    refs = np.asarray([cnum(r) for r in ref_rows])
    cands = np.asarray([cnum(r) for r in candidate_rows])
    if len(cands) == 0:
        return [None] * len(refs)
    cost = np.abs(refs[:, None] - cands[None, :])
    matched = [None] * len(refs)
    # The frozen protocol defines each tracked mode by its own nearest
    # complex-rate match.  Do not use a global assignment: when a smaller
    # dictionary merges two nearby modes, a Hungarian assignment can steal the
    # physically nearest candidate from the leading reference mode.
    for a in range(len(refs)):
        b = int(np.argmin(cost[a]))
        threshold = max(0.50, 0.35 * abs(refs[a]))
        if cost[a, b] <= threshold:
            matched[a] = candidate_rows[b]
    return matched


def stability(values):
    if any(v is None for v in values):
        return False, math.nan, math.nan, math.nan, math.nan
    z = np.asarray(values, dtype=np.complex128)
    real_range = float(np.ptp(z.real))
    imag_range = float(np.ptp(np.abs(z.imag)))
    real_tol = max(0.15, 0.20 * abs(float(np.median(z.real))))
    imag_tol = max(0.50, 0.20 * float(np.median(np.abs(z.imag))))
    return real_range <= real_tol and imag_range <= imag_tol, real_range, real_tol, imag_range, imag_tol


def fmt_values(values):
    return ';'.join('NA' if v is None else f'{v.real:.9g}{v.imag:+.9g}i' for v in values)


def bootstrap_summary(rows, mode):
    r = [x for x in rows if int(x['reference_mode']) == mode and int(x['matched']) == 1]
    re = np.asarray([float(x['lambda_real']) for x in r])
    im = np.asarray([float(x['lambda_imag']) for x in r])
    if len(re) == 0:
        return 0, (math.nan,) * 6
    return len(re), (float(np.mean(re)), *np.quantile(re, [0.025, 0.975]),
                     float(np.mean(im)), *np.quantile(im, [0.025, 0.975]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--analysis', type=Path, required=True)
    args = p.parse_args()
    root = args.analysis
    spectra = {case: read_csv(root / f'{case}_spectra.csv') for case in CASES}
    boots = {case: read_csv(root / f'{case}_bootstrap_raw.csv') for case in CASES}
    refs = {case: reference(spectra[case]) for case in CASES}

    other_dt = {
        'driven_dt1e-3': 'driven_dt2p5e-4', 'driven_dt2p5e-4': 'driven_dt1e-3',
        'equilibrium_dt1e-3': 'equilibrium_dt2p5e-4', 'equilibrium_dt2p5e-4': 'equilibrium_dt1e-3'}

    records = []
    dictionary_rows, lag_rows, timestep_rows, cutoff_rows = [], [], [], []
    for case in CASES:
        rr = refs[case]
        dict_matches = {d: global_match(rr, select(spectra[case], d, 0.5)) for d in DICTS}
        lag_matches = {t: global_match(rr, select(spectra[case], 'D3', t)) for t in LAGS}
        dt_case = other_dt[case]
        dt_match = global_match(rr, select(spectra[dt_case], 'D3', 0.5))
        cutoff_matches = {c: global_match(rr, select(spectra[case], 'D3', 0.5, c))
                          for c in (1e-8, 1e-10, 1e-12)}
        for mode, refrow in enumerate(rr):
            lam = cnum(refrow)
            dv = [None if dict_matches[d][mode] is None else cnum(dict_matches[d][mode]) for d in DICTS]
            lv = [None if lag_matches[t][mode] is None else cnum(lag_matches[t][mode]) for t in LAGS]
            tv = [lam, None if dt_match[mode] is None else cnum(dt_match[mode])]
            cv = [None if cutoff_matches[c][mode] is None else cnum(cutoff_matches[c][mode])
                  for c in (1e-8, 1e-10, 1e-12)]
            dpass, drange, dtol, dirange, ditol = stability(dv)
            lpass, lrange, ltol, lirange, litol = stability(lv)
            tpass, trange, ttol, tirange, titol = stability(tv)
            cpass, crange, ctol, cirange, citol = stability(cv)
            nboot, bs = bootstrap_summary(boots[case], mode)
            preliminary = dpass and lpass and tpass and cpass and nboot >= 450
            records.append({
                'case': case, 'mode': mode, 'lambda': lam,
                'dictionary_values': dv, 'lag_values': lv, 'timestep_values': tv,
                'cutoff_values': cv, 'dictionary_pass': dpass, 'lag_pass': lpass,
                'timestep_pass': tpass, 'cutoff_pass': cpass,
                'bootstrap_n': nboot, 'bootstrap': bs, 'preliminary': preliminary,
                'dict_stats': (drange, dtol, dirange, ditol),
                'lag_stats': (lrange, ltol, lirange, litol),
                'time_stats': (trange, ttol, tirange, titol),
                'cut_stats': (crange, ctol, cirange, citol),
            })
            for d, v in zip(DICTS, dv):
                dictionary_rows.append((case, mode, lam.real, lam.imag, d,
                                        math.nan if v is None else v.real,
                                        math.nan if v is None else v.imag,
                                        dpass, drange, dtol, dirange, ditol))
            for t, v in zip(LAGS, lv):
                lag_rows.append((case, mode, lam.real, lam.imag, t,
                                 math.nan if v is None else v.real,
                                 math.nan if v is None else v.imag,
                                 lpass, lrange, ltol, lirange, litol))
            timestep_rows.append((case, mode, lam.real, lam.imag, dt_case,
                                  math.nan if tv[1] is None else tv[1].real,
                                  math.nan if tv[1] is None else tv[1].imag,
                                  tpass, trange, ttol, tirange, titol))
            for c, v in zip((1e-8, 1e-10, 1e-12), cv):
                cutoff_rows.append((case, mode, lam.real, lam.imag, c,
                                    math.nan if v is None else v.real,
                                    math.nan if v is None else v.imag,
                                    cpass, crange, ctol, cirange, citol))

    # Neighbour separation is evaluated among candidates passing every other gate.
    for case in CASES:
        group = sorted([r for r in records if r['case'] == case and r['preliminary']],
                       key=lambda x: x['lambda'].real, reverse=True)
        for r in group:
            lo, hi = r['bootstrap'][1], r['bootstrap'][2]
            overlap = False
            for q in group:
                if q is r:
                    continue
                qlo, qhi = q['bootstrap'][1], q['bootstrap'][2]
                if max(lo, qlo) <= min(hi, qhi):
                    overlap = True
            r['neighbour_separated'] = not overlap
            r['local_reportable'] = not overlap
        for r in records:
            if r['case'] == case and 'local_reportable' not in r:
                r['neighbour_separated'] = False
                r['local_reportable'] = False

    # Timestep convergence is a two-sided requirement: the matched mode must
    # independently pass dictionary, lag, cutoff, bootstrap, and separation
    # gates in both timestep datasets.  A mode stable in only one discretization
    # is not promoted by a merely nearby point estimate in the other.
    for r in records:
        peers = [q for q in records if q['case'] == other_dt[r['case']] and q['local_reportable']]
        if peers:
            q = min(peers, key=lambda x: abs(x['lambda'] - r['lambda']))
            threshold = max(0.50, 0.35 * abs(r['lambda']))
            reverse_peers = [p for p in records if p['case'] == r['case'] and p['local_reportable']]
            reverse = min(reverse_peers, key=lambda x: abs(x['lambda'] - q['lambda']))
            reciprocal = reverse is r
            r['mutual_timestep_mode'] = q['mode'] if reciprocal and abs(q['lambda'] - r['lambda']) <= threshold else -1
        else:
            r['mutual_timestep_mode'] = -1
        r['reportable'] = r['local_reportable'] and r['mutual_timestep_mode'] >= 0

    write_csv(root / 'dictionary_convergence.csv',
              ('case','mode','reference_real','reference_imag','dictionary','lambda_real','lambda_imag',
               'pass','real_range','real_tolerance','imag_range','imag_tolerance'), dictionary_rows)
    write_csv(root / 'lag_convergence.csv',
              ('case','mode','reference_real','reference_imag','tau','lambda_real','lambda_imag',
               'pass','real_range','real_tolerance','imag_range','imag_tolerance'), lag_rows)
    write_csv(root / 'timestep_convergence.csv',
              ('case','mode','reference_real','reference_imag','comparison_case','lambda_real','lambda_imag',
               'pass','real_range','real_tolerance','imag_range','imag_tolerance'), timestep_rows)
    write_csv(root / 'cutoff_sensitivity.csv',
              ('case','mode','reference_real','reference_imag','cutoff','lambda_real','lambda_imag',
               'pass','real_range','real_tolerance','imag_range','imag_tolerance'), cutoff_rows)

    leading = []
    for r in records:
        bm, blo, bhi, im, ilo, ihi = r['bootstrap']
        leading.append((r['case'], r['mode'], r['lambda'].real, r['lambda'].imag,
                        r['bootstrap_n'], bm, blo, bhi, im, ilo, ihi,
                        r['dictionary_pass'], r['lag_pass'], r['timestep_pass'],
                        r['cutoff_pass'], r['neighbour_separated'], r['local_reportable'],
                        r['mutual_timestep_mode'], r['reportable'],
                        fmt_values(r['dictionary_values']), fmt_values(r['lag_values']),
                        fmt_values(r['timestep_values']), fmt_values(r['cutoff_values'])))
    write_csv(root / 'leading_modes_summary.csv',
              ('case','mode','lambda_real','lambda_imag','bootstrap_matches','bootstrap_mean_real',
               'ci95_real_low','ci95_real_high','bootstrap_mean_imag','ci95_imag_low','ci95_imag_high',
               'dictionary_pass','lag_pass','timestep_pass','cutoff_pass','neighbour_separated',
               'local_reportable','mutual_timestep_mode','reportable',
               'dictionary_values_D1_D2_D3','lag_values_0p1_0p25_0p5_1p0',
               'timestep_values_this_other','cutoff_values_1e8_1e10_1e12'), leading)

    # Known 5.3 oscillation diagnostic, deliberately without post-hoc rescue.
    osc_rows = []
    for case in CASES:
        base = select(spectra[case], 'D3', 0.5)
        complex_base = [r for r in base if cnum(r).imag > 1 and -20 < cnum(r).real < -0.05]
        if not complex_base:
            continue
        refrow = min(complex_base, key=lambda r: abs(cnum(r).imag - 5.3))
        refone = [refrow]
        dv = []
        for d in DICTS:
            m = global_match(refone, select(spectra[case], d, .5))[0]
            dv.append(None if m is None else cnum(m))
        lv = []
        for t in LAGS:
            m = global_match(refone, select(spectra[case], 'D3', t))[0]
            lv.append(None if m is None else cnum(m))
        dt_case = other_dt[case]
        m = global_match(refone, select(spectra[dt_case], 'D3', .5))[0]
        tv = [cnum(refrow), None if m is None else cnum(m)]
        dp = stability(dv)[0]; lp = stability(lv)[0]; tp = stability(tv)[0]
        target = abs(cnum(refrow).imag - 5.3) <= .8
        osc_rows.append((case, cnum(refrow).real, cnum(refrow).imag, target,
                         dp, lp, tp, dp and lp and tp and target,
                         fmt_values(dv), fmt_values(lv), fmt_values(tv)))
    write_csv(root / 'known_oscillation_check.csv',
              ('case','reference_real','reference_imag','frequency_target_pass','dictionary_pass',
               'lag_pass','timestep_pass','recovered','dictionary_values','lag_values','timestep_values'),
              osc_rows)

    # Trivial mode table from the diagnostics files.
    trivial = []
    for case in CASES:
        for r in read_csv(root / f'{case}_gram_diagnostics.csv'):
            trivial.append((case,r['dictionary'],r['tau'],r['cutoff'],r['kept'],r['condition'],
                            r['trivial_mu_real'],r['trivial_mu_imag'],r['trivial_error'],
                            float(r['trivial_error']) <= 1e-8))
    write_csv(root / 'trivial_mode_check.csv',
              ('case','dictionary','tau','cutoff','kept','condition','mu_real','mu_imag','error','pass'),
              trivial)

    verdict = {}
    for bath in ('driven', 'equilibrium'):
        fine = bath + '_dt2p5e-4'
        coarse = bath + '_dt1e-3'
        fine_modes = [r for r in records if r['case'] == fine and r['reportable']]
        coarse_modes = [r for r in records if r['case'] == coarse and r['reportable']]
        verdict[bath] = {
            'fine_reportable_modes': [r['mode'] for r in fine_modes],
            'fine_rates': [[r['lambda'].real,r['lambda'].imag] for r in fine_modes],
            'coarse_reportable_modes': [r['mode'] for r in coarse_modes],
            'coarse_rates': [[r['lambda'].real,r['lambda'].imag] for r in coarse_modes],
        }
    verdict['oscillation_recovered_all'] = bool(osc_rows) and all(bool(r[7]) for r in osc_rows)
    verdict['trivial_all_pass'] = all(float(r[8]) <= 1e-8 for r in trivial)
    with (root / 'gate_verdict.json').open('w') as f:
        json.dump(verdict, f, indent=2)
    print(json.dumps(verdict, indent=2))


if __name__ == '__main__':
    main()
