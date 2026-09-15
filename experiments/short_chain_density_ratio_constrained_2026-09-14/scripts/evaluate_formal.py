#!/usr/bin/env python3
"""One-shot formal evaluation of hybrid five-dimensional density ratios."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp

from common_v2 import fixed_indices, formal_evaluation_roles, load_rows
from hybrid_ratio_model import load_bundle

HERE = Path(__file__).resolve().parent.parent
OLD = HERE.parent / "short_chain_section4_2026-09-14"
sys.path.insert(0, str(OLD / "stationary_density"))
from fp_tensorflow import relative_adjoint_residual  # noqa: E402
from reduced_operator import energy  # noqa: E402

MOMENT_NAMES = (
    "I1", "I2", "I3", "M", "E", "sin_theta1", "cos_theta1",
    "sin_theta3", "cos_theta3", "I1_minus_I3",
    "I1_I2_sin_theta1", "I2_I3_sin_theta3",
)
BOOTSTRAPS = 1000
Z_GATE = 3.748286557303567


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, action="append", required=True)
    parser.add_argument("--phase", choices=("driven", "equilibrium"), required=True)
    parser.add_argument("--timestep", choices=("fine", "coarse"), required=True)
    parser.add_argument("--split", choices=("validation", "test"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--edges-from", type=Path)
    parser.add_argument("--validation-verdict", type=Path)
    return parser.parse_args()


def observables(state):
    i1, i2, i3, th1, th3 = state.T
    return np.column_stack((
        i1, i2, i3, i1+i2+i3, energy(state),
        np.sin(th1), np.cos(th1), np.sin(th3), np.cos(th3), i1-i3,
        i1*i2*np.sin(th1), i2*i3*np.sin(th3),
    ))


def log_ratio(model, state, batch=8192):
    return np.concatenate([
        model.log_ratio(tf.constant(state[start:start+batch], tf.float64)).numpy()
        for start in range(0, len(state), batch)
    ])


def residuals(model, state, t_left, t_right, batch=256):
    return np.concatenate([
        relative_adjoint_residual(
            model, tf.constant(state[start:start+batch], tf.float64),
            t_left, t_right, 0.1,
        ).numpy()
        for start in range(0, len(state), batch)
    ])


def bootstrap_moments(target_streams, reference_streams, reference_weights, salt):
    target_sum = np.asarray([observables(stream).sum(axis=0) for stream in target_streams])
    target_count = np.asarray([len(stream) for stream in target_streams], np.float64)
    weight_sum = np.asarray([weight.sum() for weight in reference_weights])
    ratio_sum = np.asarray([
        np.sum(weight[:, None]*observables(stream), axis=0)
        for stream, weight in zip(reference_streams, reference_weights)
    ])
    target_mean = target_sum.sum(axis=0)/target_count.sum()
    ratio_mean = ratio_sum.sum(axis=0)/weight_sum.sum()
    difference = ratio_mean-target_mean
    seed = int(hashlib.sha256(salt.encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    draws = np.empty((BOOTSTRAPS, len(MOMENT_NAMES)))
    for b in range(BOOTSTRAPS):
        it = rng.integers(0, len(target_streams), len(target_streams))
        ir = rng.integers(0, len(reference_streams), len(reference_streams))
        tm = target_sum[it].sum(axis=0)/target_count[it].sum()
        rm = ratio_sum[ir].sum(axis=0)/weight_sum[ir].sum()
        draws[b] = rm-tm
    return target_mean, ratio_mean, difference, draws.std(axis=0, ddof=1)


def tv(left, right):
    left = np.asarray(left, np.float64)
    right = np.asarray(right, np.float64)
    return float(0.5*np.sum(np.abs(left/left.sum()-right/right.sum())))


def create_edges(target):
    result = {}
    for index in range(3):
        low, high = np.quantile(np.log(target[:, index]), (0.001, 0.999))
        result[f"action_{index}"] = np.r_[-np.inf, np.linspace(low, high, 73)[1:-1], np.inf]
    result["angle"] = np.linspace(-np.pi, np.pi, 73)
    return result


def marginal_tvs(target, reference, weights, edges):
    values = {}
    for index, name in enumerate(("I1", "I2", "I3")):
        bins = edges[f"action_{index}"]
        ht, _ = np.histogram(np.log(target[:, index]), bins)
        hr, _ = np.histogram(np.log(reference[:, index]), bins, weights=weights)
        values[name] = tv(ht, hr)
    angle_edges = edges["angle"]
    for index, name in ((3, "theta1"), (4, "theta3")):
        ht, _ = np.histogram(target[:, index], angle_edges)
        hr, _ = np.histogram(reference[:, index], angle_edges, weights=weights)
        values[name] = tv(ht, hr)
    ht, _, _ = np.histogram2d(target[:, 3], target[:, 4], (angle_edges, angle_edges))
    hr, _, _ = np.histogram2d(
        reference[:, 3], reference[:, 4], (angle_edges, angle_edges), weights=weights
    )
    values["theta1_theta3"] = tv(ht, hr)
    return values


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    if args.split == "validation" and args.edges_from:
        raise ValueError("validation creates its own frozen edges")
    if args.split == "test":
        if not args.edges_from or not args.validation_verdict:
            raise ValueError("test requires --edges-from and --validation-verdict")
        validation = json.loads(args.validation_verdict.read_text())
        if not validation.get("complete_pass", False):
            raise SystemExit("blind test remains sealed because validation did not pass")

    tf.keras.backend.set_floatx("float64")
    target_rows, reference_rows = formal_evaluation_roles(
        args.phase, args.timestep, args.split
    )
    target, target_streams = load_rows(target_rows)
    reference, reference_streams = load_rows(reference_rows)
    if args.edges_from:
        with np.load(args.edges_from) as archive:
            edges = {key: archive[key] for key in archive.files}
    else:
        edges = create_edges(target)
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez(args.output/"marginal_edges.npz", **edges)

    t_left, t_right = (2.0, 8.0) if args.phase == "driven" else (5.0, 5.0)
    fp_index = fixed_indices(
        len(target), 10_000,
        f"formal|{args.phase}|{args.timestep}|{args.split}|fp",
    )
    agreement_index = fixed_indices(
        len(target), 50_000,
        f"formal|{args.phase}|{args.timestep}|{args.split}|agreement",
    )
    fp_state = target[fp_index]
    agreement_state = target[agreement_index]

    summaries = []
    moment_rows = []
    marginal_rows = []
    agreement_values = []
    normalizers = []
    for model_index, directory in enumerate(args.model):
        model = load_bundle(directory)
        config = json.loads((directory/"model_config.json").read_text())
        fq = log_ratio(model, reference)
        logz = float(logsumexp(fq)-math.log(len(fq)))
        shifted = fq-logz
        if not np.isfinite(shifted).all():
            raise RuntimeError(f"non-finite log weight for {directory}")
        weights = np.exp(shifted)
        if not np.isfinite(weights).all():
            raise RuntimeError(f"non-finite weight for {directory}")
        split_weights = []
        start = 0
        for stream in reference_streams:
            split_weights.append(weights[start:start+len(stream)])
            start += len(stream)

        target_mean, ratio_mean, difference, se = bootstrap_moments(
            target_streams, reference_streams, split_weights,
            f"formal-moments|{args.phase}|{args.timestep}|{args.split}|{model_index}",
        )
        z = np.divide(difference, se, out=np.full_like(difference, np.inf), where=se>0)
        label = f"seed_{int(config['seed'])}"
        for index, name in enumerate(MOMENT_NAMES):
            moment_rows.append({
                "model": label, "observable": name,
                "target_mean": target_mean[index], "ratio_mean": ratio_mean[index],
                "difference": difference[index], "stream_bootstrap_se": se[index],
                "z": z[index], "gate": Z_GATE,
                "pass": int(abs(z[index]) <= Z_GATE),
            })
        tvs = marginal_tvs(target, reference, weights, edges)
        for name, value in tvs.items():
            marginal_rows.append({
                "model": label, "marginal": name, "tv": value,
                "gate": 0.05, "pass": int(value <= 0.05),
            })
        rr = residuals(model, fp_state, t_left, t_right)
        fp_median = float(np.median(np.abs(rr)))
        fp_p90 = float(np.quantile(np.abs(rr), 0.90))
        ess_fraction = float(weights.sum()**2/np.sum(weights*weights)/len(weights))
        max_z = float(np.max(np.abs(z)))
        max_tv = float(max(tvs.values()))
        per_model_pass = bool(
            ess_fraction >= 0.10 and fp_median <= 0.10 and fp_p90 <= 0.50
            and max_z <= Z_GATE and max_tv <= 0.05
        )
        summaries.append({
            "model": label, "model_path": str(directory),
            "training_base_seed": config.get("training_base_seed"),
            "model_seed": int(config["seed"]), "log_normalizer": logz,
            "ess_fraction": ess_fraction, "fp_abs_median": fp_median,
            "fp_abs_p90": fp_p90, "max_abs_moment_z": max_z,
            "max_marginal_tv": max_tv, "per_model_pass": per_model_pass,
        })
        values = log_ratio(model, agreement_state)
        agreement_values.append(values-values.mean())
        normalizers.append({"model": label, "log_normalizer": logz})

    agreement_rows = []
    maximum_rms = 0.0
    for left in range(len(agreement_values)):
        for right in range(left+1, len(agreement_values)):
            rms = float(np.sqrt(np.mean(np.square(
                agreement_values[left]-agreement_values[right]
            ))))
            maximum_rms = max(maximum_rms, rms)
            agreement_rows.append({
                "left_model": summaries[left]["model"],
                "right_model": summaries[right]["model"],
                "centered_log_ratio_rms": rms, "gate": 0.10,
                "pass": int(rms <= 0.10),
            })
    if len(agreement_values) < 2:
        raise RuntimeError("formal evaluation requires at least two independent models")

    complete_pass = bool(
        all(row["per_model_pass"] for row in summaries) and maximum_rms <= 0.10
    )
    write_csv(args.output/"model_summary.csv", summaries)
    write_csv(args.output/"moment_details.csv", moment_rows)
    write_csv(args.output/"marginal_details.csv", marginal_rows)
    write_csv(args.output/"seed_agreement.csv", agreement_rows)
    write_csv(args.output/"normalizers.csv", normalizers)
    result = {
        "protocol_version": "section4-density-ratio-constrained-formal-v1",
        "phase": args.phase, "timestep": args.timestep, "split": args.split,
        "target_streams": [int(row["stream_id"]) for row in target_rows],
        "reference_streams": [int(row["stream_id"]) for row in reference_rows],
        "bootstrap_replicates": BOOTSTRAPS,
        "fp_sample_count": len(fp_index),
        "seed_agreement_sample_count": len(agreement_index),
        "maximum_pairwise_centered_log_ratio_rms": maximum_rms,
        "model_summaries": summaries,
        "complete_pass": complete_pass,
    }
    (args.output/"verdict.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
