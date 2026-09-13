#!/usr/bin/env python3
"""Apply the frozen gates to the targeted EDMD extension outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


CASES = ("driven_dt1e-3", "driven_dt2p5e-4",
         "equilibrium_dt1e-3", "equilibrium_dt2p5e-4")
DICTS = ("E1", "E2", "E3")
SHORT_LAGS = (0.02, 0.05, 0.10)
LONG_LAGS = (0.25, 0.50, 1.00)
REAL_LAGS = (0.10, 0.25, 0.50, 1.00)
CUTOFFS = (1e-8, 1e-10, 1e-12)
PRIMARY = 1e-10
OLD_CI = {
    "driven_dt1e-3": (-0.954628, -0.920824),
    "driven_dt2p5e-4": (-0.980007, -0.941942),
    "equilibrium_dt1e-3": (-0.972388, -0.932803),
    "equilibrium_dt2p5e-4": (-0.954830, -0.922965),
}
AUTOCORR_IMAG = {
    "driven_dt1e-3": 5.400912262786367,
    "driven_dt2p5e-4": 5.4996904921446035,
    "equilibrium_dt1e-3": 5.0809735907636755,
    "equilibrium_dt2p5e-4": 5.263668093703114,
}
OTHER_DT = {
    "driven_dt1e-3": "driven_dt2p5e-4",
    "driven_dt2p5e-4": "driven_dt1e-3",
    "equilibrium_dt1e-3": "equilibrium_dt2p5e-4",
    "equilibrium_dt2p5e-4": "equilibrium_dt1e-3",
}


def read_csv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def write_csv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def cnum(row, prefix="lambda"):
    return complex(float(row[prefix + "_real"]), float(row[prefix + "_imag"]))


def select(rows, dictionary, tau, cutoff=PRIMARY):
    return [r for r in rows
            if r["dictionary"] == dictionary
            and abs(float(r["tau"]) - tau) < 1e-12
            and abs(float(r["cutoff"]) - cutoff) <= 1e-18]


def nearest(reference, rows, real_only=False):
    candidates = []
    for r in rows:
        z = cnum(r)
        if not (np.isfinite(z.real) and np.isfinite(z.imag)) or z.imag < -1e-8:
            continue
        if real_only and abs(z.imag) > 1e-6:
            continue
        candidates.append((abs(z - reference), r))
    if not candidates:
        return None
    distance, row = min(candidates, key=lambda x: x[0])
    threshold = max(0.50, 0.35 * abs(reference))
    return row if distance <= threshold else None


def stability(values):
    if any(v is None for v in values):
        return False, math.nan, math.nan, math.nan, math.nan
    z = np.asarray(values, dtype=complex)
    real_range = float(np.ptp(z.real))
    imag_range = float(np.ptp(np.abs(z.imag)))
    real_tol = max(0.15, 0.20 * abs(float(np.median(z.real))))
    imag_tol = max(0.50, 0.20 * float(np.median(np.abs(z.imag))))
    return (real_range <= real_tol and imag_range <= imag_tol,
            real_range, real_tol, imag_range, imag_tol)


def bootstrap_stats(rows, mode=None):
    if mode is not None:
        rows = [r for r in rows if int(r["reference_mode"]) == mode]
    rows = [r for r in rows if int(r["matched"]) == 1]
    if not rows:
        return {"n": 0, "mean_real": math.nan, "lo_real": math.nan,
                "hi_real": math.nan, "mean_imag": math.nan,
                "lo_imag": math.nan, "hi_imag": math.nan}
    re = np.asarray([float(r["lambda_real"]) for r in rows])
    im = np.asarray([float(r["lambda_imag"]) for r in rows])
    return {"n": len(rows), "mean_real": float(re.mean()),
            "lo_real": float(np.quantile(re, .025)),
            "hi_real": float(np.quantile(re, .975)),
            "mean_imag": float(im.mean()),
            "lo_imag": float(np.quantile(im, .025)),
            "hi_imag": float(np.quantile(im, .975))}


def oscillatory_reference(rows):
    candidates = []
    for r in select(rows, "E3", 0.05):
        z = cnum(r)
        if 4.5 <= z.imag <= 6.0 and -20 <= z.real <= -0.05:
            candidates.append(r)
    return max(candidates, key=lambda r: cnum(r).real) if candidates else None


def leading_real_reference(rows):
    candidates = []
    for r in select(rows, "E3", 0.50):
        z, mu = cnum(r), cnum(r, "mu")
        if -3 <= z.real <= -0.05 and abs(z.imag) <= 1e-6 and 0 < abs(mu) <= 1.05:
            candidates.append(r)
    return max(candidates, key=lambda r: cnum(r).real) if candidates else None


def real_reference_rows(rows):
    out = []
    for r in select(rows, "E3", 0.50):
        z, mu = cnum(r), cnum(r, "mu")
        if -3 <= z.real <= -0.05 and -1e-8 <= z.imag <= 12 and 0 < abs(mu) <= 1.05:
            out.append(r)
        if len(out) == 12:
            break
    return out


def mode_number(rows, reference):
    refs = real_reference_rows(rows)
    if not refs:
        return None
    return int(np.argmin([abs(cnum(r) - reference) for r in refs]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--analysis", type=Path, required=True)
    args = p.parse_args()
    root = args.analysis
    spectra = {c: read_csv(root / f"{c}_spectra.csv") for c in CASES}

    osc_detail, osc_summary = [], []
    osc_records = {}
    for case in CASES:
        refrow = oscillatory_reference(spectra[case])
        if refrow is None:
            osc_records[case] = None
            osc_summary.append((case, "MISSING") + (math.nan,) * 23)
            continue
        ref = cnum(refrow)
        dict_rows = [nearest(ref, select(spectra[case], d, .05)) for d in DICTS]
        short_rows = [nearest(ref, select(spectra[case], "E3", t)) for t in SHORT_LAGS]
        long_rows = [nearest(ref, select(spectra[case], "E3", t)) for t in LONG_LAGS]
        cut_rows = [nearest(ref, select(spectra[case], "E3", .05, c)) for c in CUTOFFS]
        dict_vals = [None if r is None else cnum(r) for r in dict_rows]
        short_vals = [None if r is None else cnum(r) for r in short_rows]
        cut_vals = [None if r is None else cnum(r) for r in cut_rows]
        dstat = stability(dict_vals)
        sstat = stability(short_vals)
        cstat = stability(cut_vals)
        boot = bootstrap_stats(read_csv(root / f"{case}_bootstrap_oscillatory_raw.csv"))
        unaliased = all(r is not None and int(r["aliased_20pct"]) == 0 for r in short_rows)
        imag_excludes_zero = boot["lo_imag"] > 0 or boot["hi_imag"] < 0
        rec = {"reference": ref, "dict_rows": dict_rows, "short_rows": short_rows,
               "long_rows": long_rows, "cut_rows": cut_rows, "dstat": dstat,
               "sstat": sstat, "cstat": cstat, "boot": boot,
               "unaliased": unaliased, "imag_excludes_zero": imag_excludes_zero}
        osc_records[case] = rec
        for kind, labels, rows in (
                ("dictionary", DICTS, dict_rows),
                ("short_lag", SHORT_LAGS, short_rows),
                ("long_lag", LONG_LAGS, long_rows),
                ("cutoff", CUTOFFS, cut_rows)):
            for label, row in zip(labels, rows):
                if row is None:
                    osc_detail.append((case, kind, label, math.nan, math.nan,
                                       math.nan, math.nan, 1, 0))
                else:
                    z = cnum(row)
                    osc_detail.append((case, kind, label, z.real, z.imag,
                                       float(row["nyquist"]),
                                       abs(z - ref), int(row["aliased_20pct"]), 1))

    for case in CASES:
        rec = osc_records[case]
        if rec is None:
            continue
        other = osc_records[OTHER_DT[case]]
        dtstat = stability([rec["reference"], None if other is None else other["reference"]])
        rec["dtstat"] = dtstat
        overall = (rec["dstat"][0] and rec["sstat"][0] and rec["cstat"][0]
                   and dtstat[0] and rec["unaliased"] and rec["boot"]["n"] >= 450
                   and rec["imag_excludes_zero"])
        rec["pass"] = overall
        b = rec["boot"]
        z = rec["reference"]
        osc_summary.append((case, "PASS" if overall else "FAIL", z.real, z.imag,
                            b["mean_real"], b["lo_real"], b["hi_real"],
                            b["mean_imag"], b["lo_imag"], b["hi_imag"], b["n"],
                            AUTOCORR_IMAG[case], z.imag - AUTOCORR_IMAG[case],
                            rec["dstat"][0], rec["dstat"][1], rec["dstat"][2],
                            rec["dstat"][3], rec["dstat"][4],
                            rec["sstat"][0], rec["sstat"][1], rec["sstat"][2],
                            rec["sstat"][3], rec["sstat"][4],
                            dtstat[0], dtstat[1], dtstat[2], dtstat[3], dtstat[4],
                            rec["cstat"][0], rec["unaliased"],
                            rec["imag_excludes_zero"]))

    write_csv(root / "oscillatory_convergence.csv",
              ("case", "dimension", "setting", "lambda_real", "lambda_imag",
               "nyquist", "distance_from_reference", "aliased_20pct", "matched"),
              osc_detail)
    write_csv(root / "oscillatory_summary.csv",
              ("case", "overall_gate", "reference_real", "reference_imag",
               "bootstrap_mean_real", "bootstrap_ci_low_real", "bootstrap_ci_high_real",
               "bootstrap_mean_imag", "bootstrap_ci_low_imag", "bootstrap_ci_high_imag",
               "bootstrap_accept", "autocorr_imag", "imag_minus_autocorr",
               "dictionary_pass", "dictionary_real_range", "dictionary_real_tolerance",
               "dictionary_imag_range", "dictionary_imag_tolerance",
               "short_lag_pass", "short_lag_real_range", "short_lag_real_tolerance",
               "short_lag_imag_range", "short_lag_imag_tolerance",
               "timestep_pass", "timestep_real_range", "timestep_real_tolerance",
               "timestep_imag_range", "timestep_imag_tolerance", "cutoff_pass",
               "all_short_unaliased", "bootstrap_imag_ci_excludes_zero"), osc_summary)

    real_summary, real_detail, separation_rows = [], [], []
    real_records = {}
    for case in CASES:
        rr = leading_real_reference(spectra[case])
        if rr is None:
            continue
        ref = cnum(rr)
        dict_rows = [nearest(ref, select(spectra[case], d, .5), True) for d in DICTS]
        lag_rows = [nearest(ref, select(spectra[case], "E3", t), True) for t in REAL_LAGS]
        cut_rows = [nearest(ref, select(spectra[case], "E3", .5, c), True) for c in CUTOFFS]
        dt_row = nearest(ref, select(spectra[OTHER_DT[case]], "E3", .5), True)
        dstat = stability([None if r is None else cnum(r) for r in dict_rows])
        lstat = stability([None if r is None else cnum(r) for r in lag_rows])
        cstat = stability([None if r is None else cnum(r) for r in cut_rows])
        tstat = stability([ref, None if dt_row is None else cnum(dt_row)])
        mode = mode_number(spectra[case], ref)
        boot = bootstrap_stats(read_csv(root / f"{case}_bootstrap_real_raw.csv"), mode)
        # Nearest other real reference and its bootstrap CI for separation.
        refs = [cnum(r) for r in real_reference_rows(spectra[case]) if abs(cnum(r).imag) <= 1e-6]
        others = [z for z in refs if abs(z - ref) > 1e-7]
        neighbour = min(others, key=lambda z: abs(z.real - ref.real)) if others else None
        neighbour_boot = None
        if neighbour is not None:
            neighbour_mode = mode_number(spectra[case], neighbour)
            neighbour_boot = bootstrap_stats(read_csv(root / f"{case}_bootstrap_real_raw.csv"), neighbour_mode)
        separated = (neighbour_boot is not None and boot["n"] >= 450
                     and neighbour_boot["n"] >= 450
                     and max(boot["lo_real"], neighbour_boot["lo_real"])
                     > min(boot["hi_real"], neighbour_boot["hi_real"]))
        reportable = (dstat[0] and lstat[0] and cstat[0] and tstat[0]
                      and boot["n"] >= 450 and separated)
        oldlo, oldhi = OLD_CI[case]
        moved = boot["hi_real"] < oldlo or boot["lo_real"] > oldhi
        real_records[case] = {"reference": ref, "boot": boot, "reportable": reportable}
        real_summary.append((case, mode, ref.real, ref.imag, boot["mean_real"],
                             boot["lo_real"], boot["hi_real"], boot["n"],
                             dstat[0], lstat[0], tstat[0], cstat[0], separated,
                             reportable, oldlo, oldhi, moved,
                             math.nan if neighbour is None else neighbour.real,
                             math.nan if neighbour_boot is None else neighbour_boot["lo_real"],
                             math.nan if neighbour_boot is None else neighbour_boot["hi_real"]))
        for kind, labels, rows in (("dictionary", DICTS, dict_rows),
                                   ("lag", REAL_LAGS, lag_rows),
                                   ("cutoff", CUTOFFS, cut_rows),
                                   ("timestep", (OTHER_DT[case],), (dt_row,))):
            for label, row in zip(labels, rows):
                z = None if row is None else cnum(row)
                real_detail.append((case, kind, label,
                                    math.nan if z is None else z.real,
                                    math.nan if z is None else z.imag,
                                    0 if z is None else 1))

        # Explicitly audit the two faster driven real candidates closest to -1.8/-2.1.
        if case.startswith("driven"):
            refs_rows = real_reference_rows(spectra[case])
            real_refs = [(m, cnum(r)) for m, r in enumerate(refs_rows)
                         if abs(cnum(r).imag) <= 1e-6]
            selected = []
            used = set()
            for target in (-1.8, -2.1):
                choices = [(abs(z.real - target), m, z) for m, z in real_refs if m not in used]
                _, m, z = min(choices)
                used.add(m)
                b = bootstrap_stats(read_csv(root / f"{case}_bootstrap_real_raw.csv"), m)
                selected.append((target, m, z, b))
            b0, b1 = selected[0][3], selected[1][3]
            sep = (b0["n"] >= 450 and b1["n"] >= 450
                   and max(b0["lo_real"], b1["lo_real"])
                   > min(b0["hi_real"], b1["hi_real"]))
            for target, m, z, b in selected:
                separation_rows.append((case, target, m, z.real, b["mean_real"],
                                        b["lo_real"], b["hi_real"], b["n"], sep))

    write_csv(root / "leading_real_summary.csv",
              ("case", "reference_mode", "reference_real", "reference_imag",
               "bootstrap_mean_real", "bootstrap_ci_low", "bootstrap_ci_high",
               "bootstrap_accept", "dictionary_pass", "lag_pass", "timestep_pass",
               "cutoff_pass", "neighbour_separated", "reportable", "previous_ci_low",
               "previous_ci_high", "moved_outside_previous_ci", "nearest_real_neighbour",
               "neighbour_ci_low", "neighbour_ci_high"), real_summary)
    write_csv(root / "leading_real_convergence.csv",
              ("case", "dimension", "setting", "lambda_real", "lambda_imag", "matched"),
              real_detail)
    write_csv(root / "driven_faster_mode_separation.csv",
              ("case", "target", "reference_mode", "point_real", "bootstrap_mean_real",
               "ci_low", "ci_high", "bootstrap_accept", "pair_separated"), separation_rows)

    # Collect the leading-mode weights using the case-specific real reference number.
    weight_rows = []
    for case in CASES:
        rr = real_records.get(case)
        if rr is None:
            continue
        mode = mode_number(spectra[case], rr["reference"])
        for row in read_csv(root / f"{case}_observable_modal_weights.csv"):
            if int(row["reference_mode"]) == mode:
                weight_rows.append((case, row["observable"], float(row["weight_real"]),
                                    float(row["weight_imag"]), float(row["weight_abs"])))
    write_csv(root / "leading_modal_weights.csv",
              ("case", "observable", "weight_real", "weight_imag", "weight_abs"),
              weight_rows)

    all_osc_pass = all(osc_records[c] is not None and osc_records[c].get("pass", False)
                       for c in CASES)
    any_real_moved = any(bool(r[-4]) for r in real_summary)
    if all_osc_pass and not any_real_moved:
        verdict = "OSCILLATION_RECOVERED_LEADING_UNCHANGED"
    elif all_osc_pass and any_real_moved:
        verdict = "OSCILLATION_RECOVERED_LEADING_MOVED"
    else:
        verdict = "SPECTRUM_INCOMPLETE_OSCILLATION_GATE_FAILED"
    summary = {
        "verdict": verdict,
        "all_four_oscillatory_gates_pass": all_osc_pass,
        "any_leading_interval_moved": any_real_moved,
        "oscillatory": {c: None if osc_records[c] is None else {
            "reference_real": osc_records[c]["reference"].real,
            "reference_imag": osc_records[c]["reference"].imag,
            "pass": osc_records[c].get("pass", False),
            "dictionary_pass": osc_records[c]["dstat"][0],
            "short_lag_pass": osc_records[c]["sstat"][0],
            "timestep_pass": osc_records[c]["dtstat"][0],
            "cutoff_pass": osc_records[c]["cstat"][0],
            "bootstrap": osc_records[c]["boot"],
        } for c in CASES},
    }
    with (root / "primary_results.json").open("w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
