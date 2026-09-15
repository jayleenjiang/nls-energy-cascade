#!/usr/bin/env python3
"""Blind equilibrium evaluation of the validation-selected normalized density."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE/"stationary_density"))

from fp_tensorflow import relative_adjoint_residual  # noqa: E402
from normalized_density import load_bundle  # noqa: E402
from reduced_operator import energy  # noqa: E402

CASE = "equilibrium_dt2p5e-4"
N_SNAP = 100_001
N_COL = 12
MODEL_SAMPLES = 1_000_000
MODEL_SAMPLE_SEED = 2026091410
FP_POINTS = 20_000
FP_SELECTION_SEED = 2026091411
PAIRWISE_POINTS = 50_000
PAIRWISE_SELECTION_SEED = 2026091412
SEEDS = (4101, 4102, 4103)
MOMENT_NAMES = (
    "I1", "I2", "I3", "M", "E", "sin_theta1", "cos_theta1",
    "sin_theta3", "cos_theta3", "I1_minus_I3",
    "I1_I2_sin_theta1", "I2_I3_sin_theta3",
)


def candidate_dir(choice: dict, seed: int) -> Path:
    widths = "x".join(str(x) for x in choice["widths"])
    lam = str(choice["lambda_fp"]).replace(".", "p")
    return (HERE/"stationary_density"/"equilibrium_candidates"/
            f"K{choice['components']}_W{widths}_L{lam}_S{seed}")


def split_rows(split: str) -> list[dict[str, str]]:
    with (HERE/"STREAM_SPLITS.csv").open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle)
                if row["branch"] == "density" and row["case"] == CASE
                and row["split"] == split]
    return sorted(rows, key=lambda row: int(row["stream_id"]))


def load_state(row: dict[str, str]) -> np.ndarray:
    raw = np.fromfile(REPO/row["relative_path"], dtype="<f4").reshape(N_SNAP, N_COL)
    out = np.empty((N_SNAP, 5), dtype=np.float64)
    out[:, 0] = 0.5*(raw[:, 6].astype(np.float64)+raw[:, 7].astype(np.float64))
    out[:, 1] = raw[:, 5]
    out[:, 2] = 0.5*(raw[:, 6].astype(np.float64)-raw[:, 7].astype(np.float64))
    out[:, 3] = np.arctan2(raw[:, 1], raw[:, 0])
    out[:, 4] = np.arctan2(raw[:, 3], raw[:, 2])
    return out


def log_prob(model, state: np.ndarray, batch: int = 16_384) -> np.ndarray:
    values = []
    for start in range(0, len(state), batch):
        x = tf.constant(state[start:start+batch], tf.float64)
        values.append(model.log_prob(x).numpy())
    return np.concatenate(values)


def fp_residual(model, state: np.ndarray, batch: int = 512) -> np.ndarray:
    values = []
    for start in range(0, len(state), batch):
        x = tf.constant(state[start:start+batch], tf.float64)
        values.append(relative_adjoint_residual(model, x, 5.0, 5.0, 0.1).numpy())
    return np.concatenate(values)


def observables(state: np.ndarray) -> np.ndarray:
    i1, i2, i3, th1, th3 = state.T
    return np.column_stack((
        i1, i2, i3, i1+i2+i3, energy(state),
        np.sin(th1), np.cos(th1), np.sin(th3), np.cos(th3), i1-i3,
        i1*i2*np.sin(th1), i2*i3*np.sin(th3),
    ))


def mean_se_by_stream(stream_values: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    means = np.asarray([x.mean(axis=0) for x in stream_values])
    return means.mean(axis=0), means.std(axis=0, ddof=1)/math.sqrt(len(means))


def tv_from_counts(a: np.ndarray, b: np.ndarray) -> float:
    pa = a/a.sum(); pb = b/b.sum()
    return float(0.5*np.sum(np.abs(pa-pb)))


def fixed_marginals(validation: np.ndarray, test: np.ndarray,
                    generated: np.ndarray) -> list[dict]:
    rows = []
    for j, name in enumerate(("I1", "I2", "I3")):
        lo, hi = np.quantile(np.log(validation[:, j]), (0.001, 0.999))
        core = np.linspace(lo, hi, 73)
        edges = np.concatenate(([-np.inf], core[1:-1], [np.inf]))
        ct, _ = np.histogram(np.log(test[:, j]), edges)
        cm, _ = np.histogram(np.log(generated[:, j]), edges)
        rows.append({
            "marginal": name, "kind": "log_action_72_bins_with_flow",
            "validation_log_low": lo, "validation_log_high": hi,
            "test_underflow": int(np.sum(np.log(test[:, j]) < lo)),
            "test_overflow": int(np.sum(np.log(test[:, j]) > hi)),
            "model_underflow": int(np.sum(np.log(generated[:, j]) < lo)),
            "model_overflow": int(np.sum(np.log(generated[:, j]) > hi)),
            "tv": tv_from_counts(ct, cm),
        })
    angle_edges = np.linspace(-np.pi, np.pi, 73)
    for j, name in ((3, "theta1"), (4, "theta3")):
        ct, _ = np.histogram(test[:, j], angle_edges)
        cm, _ = np.histogram(generated[:, j], angle_edges)
        rows.append({
            "marginal": name, "kind": "periodic_angle_72_bins",
            "validation_log_low": "", "validation_log_high": "",
            "test_underflow": 0, "test_overflow": 0,
            "model_underflow": 0, "model_overflow": 0,
            "tv": tv_from_counts(ct, cm),
        })
    ct, _, _ = np.histogram2d(test[:, 3], test[:, 4], (angle_edges, angle_edges))
    cm, _, _ = np.histogram2d(
        generated[:, 3], generated[:, 4], (angle_edges, angle_edges)
    )
    rows.append({
        "marginal": "theta1_theta3", "kind": "periodic_angle_72x72_bins",
        "validation_log_low": "", "validation_log_high": "",
        "test_underflow": 0, "test_overflow": 0,
        "model_underflow": 0, "model_overflow": 0,
        "tv": tv_from_counts(ct, cm),
    })
    return rows


def write_rows(path: Path, rows: list[dict]) -> None:
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    tf.keras.backend.set_floatx("float64")
    choice_document = json.loads(
        (HERE/"stationary_density"/"FROZEN_DENSITY_CHOICE.json").read_text()
    )
    if choice_document["status"] != "FROZEN":
        raise RuntimeError("density choice is not frozen")
    choice = choice_document["selected"]
    weak = json.loads((HERE/"analysis"/"weak_identity_calibration.json").read_text())
    z_threshold = float(weak["familywise_z_threshold"])
    weak_rows = list(csv.DictReader(
        (HERE/"analysis"/"weak_identity_blind_tests.csv").open()
    ))
    weak_pass = all(int(row["pass"]) for row in weak_rows if row["case"] == CASE)

    validation_streams = [load_state(row) for row in split_rows("validation")]
    test_streams = [load_state(row) for row in split_rows("test")]
    validation = np.concatenate(validation_streams)
    test = np.concatenate(test_streams)
    test_observables = [observables(x) for x in test_streams]
    trajectory_mean, trajectory_se = mean_se_by_stream(test_observables)
    trajectory_second, trajectory_second_se = mean_se_by_stream(
        [x*x for x in test_observables]
    )

    fp_rng = np.random.default_rng(FP_SELECTION_SEED)
    fp_index = fp_rng.choice(len(test), size=FP_POINTS, replace=False)
    pair_rng = np.random.default_rng(PAIRWISE_SELECTION_SEED)
    pair_index = pair_rng.choice(len(test), size=PAIRWISE_POINTS, replace=False)

    seed_summaries = []
    moment_rows = []
    marginal_rows = []
    nll_rows = []
    density_error_rows = []
    pair_logp = []
    for seed in SEEDS:
        directory = candidate_dir(choice, seed)
        model = load_bundle(directory)
        validation_metrics = json.loads((directory/"validation_metrics.json").read_text())
        generated = model.sample_numpy(MODEL_SAMPLES, MODEL_SAMPLE_SEED+seed)
        generated_obs = observables(generated)
        generated_mean = generated_obs.mean(axis=0)
        generated_se = generated_obs.std(axis=0, ddof=1)/math.sqrt(MODEL_SAMPLES)
        generated_second = np.mean(generated_obs*generated_obs, axis=0)
        generated_second_se = np.std(
            generated_obs*generated_obs, axis=0, ddof=1
        )/math.sqrt(MODEL_SAMPLES)
        moment_pass = True
        for stat, ref, ref_se, estimate, estimate_se in (
            ("mean", trajectory_mean, trajectory_se, generated_mean, generated_se),
            ("second_moment", trajectory_second, trajectory_second_se,
             generated_second, generated_second_se),
        ):
            combined = np.sqrt(ref_se*ref_se+estimate_se*estimate_se)
            z = (estimate-ref)/combined
            for q, name in enumerate(MOMENT_NAMES):
                passed = bool(abs(z[q]) <= z_threshold)
                moment_pass &= passed
                moment_rows.append({
                    "seed": seed, "observable": name, "statistic": stat,
                    "trajectory": ref[q], "trajectory_se": ref_se[q],
                    "model": estimate[q], "model_se": estimate_se[q],
                    "difference": estimate[q]-ref[q], "combined_z": z[q],
                    "familywise_z_threshold": z_threshold, "pass": int(passed),
                })

        marginals = fixed_marginals(validation, test, generated)
        circular_pass = all(
            row["tv"] <= 0.05 for row in marginals
            if row["marginal"] in ("theta1", "theta3", "theta1_theta3")
        )
        for row in marginals:
            marginal_rows.append({"seed": seed, **row,
                                  "tv_gate": 0.05,
                                  "gate_applies": int(row["marginal"].startswith("theta")),
                                  "pass": int(row["tv"] <= 0.05)})

        stream_nll = []
        for row, state in zip(split_rows("test"), test_streams):
            value = float(-np.mean(log_prob(model, state)))
            stream_nll.append(value)
            nll_rows.append({"seed": seed, "stream_id": row["stream_id"],
                             "test_nll": value, "finite": int(np.isfinite(value))})
        lp = log_prob(model, test)
        exact_shape = -energy(test)/5.0
        bounds = validation_metrics["gibbs_bulk_bounds"]
        bulk = (exact_shape >= bounds[0]) & (exact_shape <= bounds[1])
        slope, intercept = np.polyfit(exact_shape[bulk], lp[bulk], 1)
        additive = float(validation_metrics["gibbs_validation_additive_constant"])
        error = lp-exact_shape-additive
        bulk_rmse = float(np.sqrt(np.mean(error[bulk]*error[bulk])))
        percentiles = {
            int(key): float(value)
            for key, value in validation_metrics["gibbs_validation_percentiles"].items()
        }
        for low, high in ((1, 5), (5, 25), (25, 75), (75, 95), (95, 99)):
            selected = ((exact_shape >= percentiles[low]) &
                        (exact_shape < percentiles[high] if high < 99
                         else exact_shape <= percentiles[high]))
            band_error = error[selected]
            density_error_rows.append({
                "seed": seed, "exact_shape_percentile_low": low,
                "exact_shape_percentile_high": high,
                "points": int(selected.sum()),
                "centered_log_error_mean": float(np.mean(band_error)),
                "centered_log_error_rms": float(np.sqrt(np.mean(band_error*band_error))),
                "centered_log_error_p90_abs": float(np.quantile(np.abs(band_error), 0.90)),
            })
        residual = fp_residual(model, test[fp_index])
        fp_median = float(np.median(np.abs(residual)))
        fp_p90 = float(np.quantile(np.abs(residual), 0.90))
        pair_logp.append(log_prob(model, test[pair_index]))
        seed_pass = bool(
            np.isfinite(lp).mean() >= 0.99
            and 0.98 <= slope <= 1.02
            and bulk_rmse <= 0.15
            and fp_median <= 0.10 and fp_p90 <= 0.50
            and moment_pass and circular_pass and weak_pass
            and np.isfinite(stream_nll).all()
        )
        seed_summaries.append({
            "seed": seed, "test_rows": len(test),
            "finite_log_density_fraction": float(np.isfinite(lp).mean()),
            "gibbs_test_slope": float(slope),
            "gibbs_test_intercept": float(intercept),
            "validation_frozen_additive_constant": additive,
            "gibbs_bulk_centered_log_rmse": bulk_rmse,
            "fp_abs_median": fp_median, "fp_abs_p90": fp_p90,
            "moment_family_pass": moment_pass,
            "circular_marginals_pass": circular_pass,
            "weak_identities_pass": weak_pass,
            "all_seed_gates_pass": seed_pass,
        })
        del model, generated, generated_obs, lp

    pairwise_rows = []
    pairwise_pass = True
    for i in range(len(SEEDS)):
        for j in range(i+1, len(SEEDS)):
            difference = pair_logp[i]-pair_logp[j]
            difference -= difference.mean()
            rms = float(np.sqrt(np.mean(difference*difference)))
            pairwise_pass &= rms <= 0.10
            pairwise_rows.append({
                "seed_a": SEEDS[i], "seed_b": SEEDS[j],
                "centered_log_density_rms": rms, "gate": 0.10,
                "pass": int(rms <= 0.10),
            })

    all_pass = bool(all(x["all_seed_gates_pass"] for x in seed_summaries)
                    and pairwise_pass)
    output = HERE/"analysis"
    output.mkdir(exist_ok=True)
    write_rows(output/"equilibrium_density_seed_summary.csv", seed_summaries)
    write_rows(output/"equilibrium_density_moments.csv", moment_rows)
    write_rows(output/"equilibrium_density_marginals.csv", marginal_rows)
    write_rows(output/"equilibrium_density_test_nll_by_stream.csv", nll_rows)
    write_rows(output/"equilibrium_density_error_bands.csv", density_error_rows)
    write_rows(output/"equilibrium_density_seed_pairwise.csv", pairwise_rows)
    verdict = {
        "protocol_version": "section4-v1",
        "case": CASE,
        "selected": choice,
        "model_samples_per_seed": MODEL_SAMPLES,
        "model_sample_seed_base": MODEL_SAMPLE_SEED,
        "fp_test_points": FP_POINTS,
        "familywise_z_threshold": z_threshold,
        "seed_summaries": seed_summaries,
        "pairwise": pairwise_rows,
        "blind_equilibrium_density_pass": all_pass,
    }
    (output/"equilibrium_density_verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True)+"\n"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
