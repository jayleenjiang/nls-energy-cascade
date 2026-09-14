#!/usr/bin/env python3
"""Apply frozen driven-validation gates and freeze one architecture."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import tensorflow as tf

from density_ratio_common import load_rows, target_reference_rows
from gibbs_ratio_model import load_bundle

HERE = Path(__file__).resolve().parent.parent
OLD = HERE.parent/"short_chain_section4_2026-09-14"
import sys
sys.path.insert(0, str(OLD/"stationary_density"))
from reduced_operator import energy  # noqa: E402

MOMENT_NAMES = (
    "I1", "I2", "I3", "M", "E", "sin_theta1", "cos_theta1",
    "sin_theta3", "cos_theta3", "I1_minus_I3",
    "I1_I2_sin_theta1", "I2_I3_sin_theta3",
)
BOOTSTRAPS = 1000


def observables(state):
    i1, i2, i3, th1, th3 = state.T
    return np.column_stack((
        i1, i2, i3, i1+i2+i3, energy(state),
        np.sin(th1), np.cos(th1), np.sin(th3), np.cos(th3), i1-i3,
        i1*i2*np.sin(th1), i2*i3*np.sin(th3),
    ))


def log_ratio(model, state, batch=8192):
    return np.concatenate([
        model.log_ratio(tf.constant(state[s:s+batch], tf.float64)).numpy()
        for s in range(0, len(state), batch)
    ])


def bootstrap_moments(target_streams, reference_streams, reference_weights, seed):
    p_sum = np.asarray([observables(x).sum(axis=0) for x in target_streams])
    p_count = np.asarray([len(x) for x in target_streams], dtype=np.float64)
    q_weight = np.asarray([w.sum() for w in reference_weights])
    q_sum = np.asarray([
        np.sum(w[:, None]*observables(x), axis=0)
        for x, w in zip(reference_streams, reference_weights)
    ])
    p_mean = p_sum.sum(axis=0)/p_count.sum()
    q_mean = q_sum.sum(axis=0)/q_weight.sum()
    difference = q_mean-p_mean
    rng = np.random.default_rng(seed)
    draws = np.empty((BOOTSTRAPS, len(MOMENT_NAMES)))
    for b in range(BOOTSTRAPS):
        ip = rng.integers(0, len(p_sum), len(p_sum))
        iq = rng.integers(0, len(q_sum), len(q_sum))
        bp = p_sum[ip].sum(axis=0)/p_count[ip].sum()
        bq = q_sum[iq].sum(axis=0)/q_weight[iq].sum()
        draws[b] = bq-bp
    se = draws.std(axis=0, ddof=1)
    return p_mean, q_mean, difference, se


def tv(a, b):
    pa = a/a.sum(); pb = b/b.sum()
    return float(0.5*np.sum(np.abs(pa-pb)))


def marginal_metrics(target, reference, weights, action_edges, angle_edges):
    out = []
    for j, name in enumerate(("I1", "I2", "I3")):
        hp, _ = np.histogram(np.log(target[:, j]), action_edges[j])
        hq, _ = np.histogram(np.log(reference[:, j]), action_edges[j], weights=weights)
        out.append((name, tv(hp, hq)))
    for j, name in ((3, "theta1"), (4, "theta3")):
        hp, _ = np.histogram(target[:, j], angle_edges)
        hq, _ = np.histogram(reference[:, j], angle_edges, weights=weights)
        out.append((name, tv(hp, hq)))
    hp, _, _ = np.histogram2d(target[:, 3], target[:, 4], (angle_edges, angle_edges))
    hq, _, _ = np.histogram2d(
        reference[:, 3], reference[:, 4], (angle_edges, angle_edges), weights=weights
    )
    out.append(("theta1_theta3", tv(hp, hq)))
    return out


def write_csv(path, rows):
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def main():
    equilibrium = json.loads((HERE/"EQUILIBRIUM_BLIND_GATE.json").read_text())
    if not equilibrium["pass"]:
        raise SystemExit("equilibrium blind gate did not pass")
    admissible = json.loads((HERE/"EQUILIBRIUM_ADMISSIBLE_SET.json").read_text())[
        "admissible_architectures"
    ]
    target_rows, reference_rows = target_reference_rows("driven", "fine", "validation")
    target, target_streams = load_rows(target_rows)
    reference, reference_streams = load_rows(reference_rows)
    action_edges = []
    for j in range(3):
        low, high = np.quantile(np.log(target[:, j]), (0.001, 0.999))
        action_edges.append(np.r_[-np.inf, np.linspace(low, high, 73)[1:-1], np.inf])
    angle_edges = np.linspace(-np.pi, np.pi, 73)
    weak = json.loads((OLD/"analysis/weak_identity_calibration.json").read_text())
    z_threshold = float(weak["familywise_z_threshold"])

    analysis = HERE/"analysis"
    analysis.mkdir(exist_ok=True)
    moment_rows = []; marginal_rows = []; summaries = []
    normalized_predictions = {}
    for architecture in admissible:
        arch_rows = []
        for seed in (5201, 5202, 5203):
            directory = HERE/"driven_fine_candidates"/f"{architecture}_seed{seed}"
            metrics = json.loads((directory/"validation_metrics.json").read_text())
            config = json.loads((directory/"model_config.json").read_text())
            model = load_bundle(directory)
            fp = log_ratio(model, target)
            fq = log_ratio(model, reference)
            logz = float(config["validation_logz"])
            normalized_predictions[(architecture, seed)] = fp[:50_000]-logz
            weights = np.exp(fq-logz)
            weight_streams = []
            start = 0
            for stream in reference_streams:
                weight_streams.append(weights[start:start+len(stream)])
                start += len(stream)
            salt = int(hashlib.sha256(
                f"section4-density-ratio-v1|moments|{architecture}|{seed}".encode()
            ).hexdigest()[:16], 16)
            p_mean, q_mean, diff, se = bootstrap_moments(
                target_streams, reference_streams, weight_streams, salt
            )
            moment_pass = True
            for j, name in enumerate(MOMENT_NAMES):
                z = diff[j]/se[j] if se[j] > 0 else np.inf
                passed = bool(abs(z) <= z_threshold)
                moment_pass &= passed
                moment_rows.append({
                    "architecture": architecture, "seed": seed, "observable": name,
                    "target_mean": p_mean[j], "ratio_mean": q_mean[j],
                    "difference": diff[j], "stream_bootstrap_se": se[j],
                    "z": z, "familywise_z_threshold": z_threshold, "pass": int(passed),
                })
            marginals = marginal_metrics(
                target, reference, weights, action_edges, angle_edges
            )
            marginal_pass = all(value <= 0.05 for _, value in marginals)
            for name, value in marginals:
                marginal_rows.append({
                    "architecture": architecture, "seed": seed, "marginal": name,
                    "tv": value, "gate": 0.05, "pass": int(value <= 0.05),
                })
            item = {
                "architecture": architecture, "seed": seed,
                "relative_validation_nll": float(-np.mean(fp-logz)),
                "ess_fraction": metrics["reference_weight_ess_fraction"],
                "max_weight_share": metrics["max_weight_share"],
                "fp_abs_median": metrics["fp_abs_median"],
                "fp_abs_p90": metrics["fp_abs_p90"],
                "moment_pass": moment_pass, "marginal_pass": marginal_pass,
                "finite": metrics["finite"],
            }
            item["individual_pass"] = bool(
                item["ess_fraction"] >= 0.10 and item["fp_abs_median"] <= 0.10
                and item["fp_abs_p90"] <= 0.50 and moment_pass and marginal_pass
                and item["finite"]
            )
            arch_rows.append(item)
        pairwise = []
        for left, right in ((5201, 5202), (5201, 5203), (5202, 5203)):
            delta = normalized_predictions[(architecture, left)]-normalized_predictions[(architecture, right)]
            pairwise.append({"left": left, "right": right,
                             "centered_log_ratio_rms": float(np.sqrt(np.mean(delta*delta)))})
        pair_pass = all(x["centered_log_ratio_rms"] <= 0.10 for x in pairwise)
        summaries.append({
            "architecture": architecture, "seeds": arch_rows, "pairwise": pairwise,
            "pairwise_pass": pair_pass,
            "all_seed_admissible": bool(pair_pass and all(x["individual_pass"] for x in arch_rows)),
            "median_relative_validation_nll": float(np.median([
                x["relative_validation_nll"] for x in arch_rows
            ])),
            "relative_validation_nll_se": float(np.std([
                x["relative_validation_nll"] for x in arch_rows
            ], ddof=1)/math.sqrt(3)),
        })
    eligible = [x for x in summaries if x["all_seed_admissible"]]
    selected = None
    if eligible:
        eligible.sort(key=lambda x: (x["median_relative_validation_nll"],
                                    ("linear", "mlp32", "mlp64").index(x["architecture"])))
        selected = eligible[0]["architecture"]
    result = {
        "protocol_version": "section4-density-ratio-v1",
        "status": "FROZEN" if selected else "NO_ADMISSIBLE_DRIVEN_CANDIDATE",
        "selected_architecture": selected, "candidate_summaries": summaries,
        "familywise_z_threshold": z_threshold,
        "target_validation_streams": [int(x["stream_id"]) for x in target_rows],
        "reference_validation_streams": [int(x["stream_id"]) for x in reference_rows],
        "action_log_edges": [x.tolist() for x in action_edges],
        "angle_edges": angle_edges.tolist(),
    }
    write_csv(analysis/"driven_validation_moments.csv", moment_rows)
    write_csv(analysis/"driven_validation_marginals.csv", marginal_rows)
    (HERE/"FROZEN_DRIVEN_CHOICE.json").write_text(
        json.dumps(result, indent=2, sort_keys=True)+"\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if selected else 2)


if __name__ == "__main__":
    tf.keras.backend.set_floatx("float64")
    main()
