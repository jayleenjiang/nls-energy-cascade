#!/usr/bin/env python3
"""Fit Section-4.3 EDMD on train streams and select on validation streams."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from mode_common import (
    ROOT, CASES, CUTOFFS, DICTIONARIES, LAGS, aggregate, base, extended,
    fit, nearest, select_complex, select_real, split_paths, stability,
    train_log_stats, validation_residual,
)


OUT = ROOT/"modes"
MATRIX = OUT/"train_validation_matrices"
PRIMARY = {"complex": 0.05, "real": 0.50}
LAG_SET = {"complex": (0.02, 0.05, 0.10), "real": (0.25, 0.50, 1.00)}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def candidate_record(case, kind, dictionary, cutoff, selected, result,
                     validation, fits):
    if selected is None:
        return None
    index, mu, rate = selected
    coefficient = result.coeff[:, index]
    tau = PRIMARY[kind]
    residual = validation_residual(
        validation, dictionary, tau, mu, coefficient
    )
    dict_values = []
    for name in ("E2", "E3"):
        matched = nearest(rate, fits[(name, tau, cutoff)])
        if matched is None:
            dict_values = []
            break
        dict_values.append(matched[2])
    dictionary_pass, dictionary_detail = stability(dict_values)

    lag_values = []
    unaliased = True
    for lag in LAG_SET[kind]:
        matched = nearest(rate, fits[("E3", lag, cutoff)])
        if matched is None:
            lag_values = []
            break
        lag_values.append(matched[2])
        if kind == "complex" and abs(matched[2].imag) >= 0.8*np.pi/lag:
            unaliased = False
    lag_pass, lag_detail = stability(lag_values)
    lag_pass = lag_pass and unaliased

    cutoff_values = []
    for candidate_cutoff in CUTOFFS:
        matched = nearest(rate, fits[("E3", tau, candidate_cutoff)])
        if matched is None:
            cutoff_values = []
            break
        cutoff_values.append(matched[2])
    cutoff_pass, cutoff_detail = stability(cutoff_values)
    conditioning_pass = bool(
        result.trivial_error <= 1.0e-8 and result.condition <= 1.01/cutoff
    )
    return {
        "candidate_id": f"{case}:{kind}:{dictionary}:{cutoff:.0e}",
        "case": case, "kind": kind, "dictionary": dictionary,
        "nominal_K": len(extended.DICT_INDICES[dictionary]),
        "cutoff": cutoff, "tau": tau, "eigen_index": int(index),
        "mu": [float(mu.real), float(mu.imag)],
        "lambda": [float(rate.real), float(rate.imag)],
        "validation_residual": residual,
        "gram_kept": result.kept, "gram_condition": result.condition,
        "trivial_error": result.trivial_error,
        "conditioning_pass": conditioning_pass,
        "dictionary_pass": dictionary_pass,
        "dictionary_detail": dictionary_detail,
        "lag_pass": lag_pass, "lag_detail": lag_detail,
        "all_short_unaliased": unaliased if kind == "complex" else None,
        "cutoff_pass": cutoff_pass, "cutoff_detail": cutoff_detail,
        "local_admissible": bool(
            conditioning_pass and dictionary_pass and lag_pass and cutoff_pass
        ),
    }


def process_case(case: str):
    train_paths = split_paths(case, "train")
    validation_paths = split_paths(case, "validation")
    if len(train_paths) != 32 or len(validation_paths) != 16:
        raise RuntimeError(f"bad mode split counts for {case}")
    mean, sd = train_log_stats(train_paths)
    training = aggregate(train_paths, mean, sd, f"{case}/train")
    validation = aggregate(validation_paths, mean, sd, f"{case}/validation")
    MATRIX.mkdir(parents=True, exist_ok=True)
    np.savez(
        MATRIX/f"{case}.npz", log_mean=mean, log_sd=sd,
        train_A=training.A, train_B=training.B, train_C=training.C,
        validation_A=validation.A, validation_B=validation.B,
        validation_C=validation.C,
    )
    fits = {}
    spectra = []
    for dictionary in DICTIONARIES:
        for cutoff in CUTOFFS:
            for tau in LAGS:
                result = fit(training, dictionary, tau, cutoff, vectors=True)
                fits[(dictionary, tau, cutoff)] = result
                for rank, (index, mu, rate) in enumerate(base.canonical_modes(result)):
                    spectra.append({
                        "case": case, "dictionary": dictionary, "cutoff": cutoff,
                        "tau": tau, "rank": rank, "eigen_index": index,
                        "mu_real": mu.real, "mu_imag": mu.imag,
                        "lambda_real": rate.real, "lambda_imag": rate.imag,
                        "aliased_20pct": int(abs(rate.imag) >= 0.8*np.pi/tau),
                        "gram_kept": result.kept,
                        "gram_condition": result.condition,
                        "trivial_error": result.trivial_error,
                    })
    candidates = []
    for kind, selector in (("complex", select_complex), ("real", select_real)):
        tau = PRIMARY[kind]
        for dictionary in DICTIONARIES:
            for cutoff in CUTOFFS:
                result = fits[(dictionary, tau, cutoff)]
                record = candidate_record(
                    case, kind, dictionary, cutoff, selector(result), result,
                    validation, fits,
                )
                if record is not None:
                    candidates.append(record)
    with (OUT/f"{case}_spectra.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(spectra[0]))
        writer.writeheader(); writer.writerows(spectra)
    return candidates


def choose_pair(coarse: list[dict], fine: list[dict]):
    pairs = []
    for left in coarse:
        for right in fine:
            passed, detail = stability([
                complex(*left["lambda"]), complex(*right["lambda"])
            ])
            if passed:
                pairs.append((left["validation_residual"]+right["validation_residual"],
                              left, right, detail))
    pair_admissible = bool(pairs)
    if not pairs:
        left = min(coarse, key=lambda item: item["validation_residual"])
        right = min(fine, key=lambda item: item["validation_residual"])
        _, detail = stability([complex(*left["lambda"]), complex(*right["lambda"])])
        return left, right, False, detail
    best_score = min(item[0] for item in pairs)
    tie = [item for item in pairs if item[0] <= 1.05*best_score]
    order = {"E1": 0, "E2": 1, "E3": 2}
    tie.sort(key=lambda item: (
        order[item[1]["dictionary"]]+order[item[2]["dictionary"]],
        -item[1]["cutoff"]-item[2]["cutoff"], item[0],
    ))
    _, left, right, detail = tie[0]
    return left, right, pair_admissible, detail


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_candidates = []
    for case in CASES:
        all_candidates.extend(process_case(case))
    with (OUT/"validation_candidates.json").open("w") as handle:
        json.dump(all_candidates, handle, indent=2, sort_keys=True)

    selections = {}
    for bath in ("driven", "equilibrium"):
        coarse_case = f"{bath}_dt1e-3"
        fine_case = f"{bath}_dt2p5e-4"
        for kind in ("complex", "real"):
            coarse = [r for r in all_candidates if r["case"] == coarse_case
                      and r["kind"] == kind and r["local_admissible"]]
            fine = [r for r in all_candidates if r["case"] == fine_case
                    and r["kind"] == kind and r["local_admissible"]]
            if not coarse:
                coarse = [r for r in all_candidates if r["case"] == coarse_case
                          and r["kind"] == kind]
            if not fine:
                fine = [r for r in all_candidates if r["case"] == fine_case
                        and r["kind"] == kind]
            left, right, pair_pass, detail = choose_pair(coarse, fine)
            selections[left["case"]+":"+kind] = {
                **left, "timestep_pair_admissible": pair_pass,
                "timestep_detail": detail,
            }
            selections[right["case"]+":"+kind] = {
                **right, "timestep_pair_admissible": pair_pass,
                "timestep_detail": detail,
            }

    # Refit only selected training eigenvectors and freeze them before test.
    fit_dir = OUT/"frozen_fits"
    fit_dir.mkdir(exist_ok=True)
    for key, selected in selections.items():
        case, kind = key.split(":")
        data = np.load(MATRIX/f"{case}.npz")
        from mode_common import Aggregate
        training = Aggregate(data["train_A"], data["train_B"], data["train_C"],
                             32*2048, 32)
        result = fit(training, selected["dictionary"], selected["tau"],
                     selected["cutoff"], vectors=True)
        target = complex(*selected["lambda"])
        matched = nearest(target, result)
        if matched is None:
            raise RuntimeError(f"cannot refit frozen selection {key}")
        index, mu, rate = matched
        path = fit_dir/f"{case}_{kind}.npz"
        np.savez(
            path, coefficient=result.coeff[:, index],
            dictionary_indices=extended.DICT_INDICES[selected["dictionary"]],
            mu=np.asarray([mu.real, mu.imag]),
            rate=np.asarray([rate.real, rate.imag]),
            log_mean=data["log_mean"], log_sd=data["log_sd"],
        )
        selected["frozen_fit"] = str(path.relative_to(ROOT))
        selected["frozen_fit_sha256"] = sha256(path)

    selection_path = OUT/"MODE_SELECTION.json"
    with selection_path.open("w") as handle:
        json.dump({
            "protocol_version": "section4-v1",
            "test_streams_read": False,
            "selections": selections,
        }, handle, indent=2, sort_keys=True)
    (OUT/"MODE_SELECTION.sha256").write_text(
        f"{sha256(selection_path)}  {selection_path.name}\n"
    )
    print(json.dumps({
        key: {
            "lambda": value["lambda"],
            "validation_residual": value["validation_residual"],
            "local_admissible": value["local_admissible"],
            "timestep_pair_admissible": value["timestep_pair_admissible"],
        } for key, value in selections.items()
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

