#!/usr/bin/env python3
"""One-shot blind known-answer gate for equilibrium-admissible models."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp
from sklearn.metrics import roc_auc_score

from density_ratio_common import fixed_indices, load_rows, target_reference_rows
from gibbs_ratio_model import load_bundle

HERE = Path(__file__).resolve().parent.parent
OLD = HERE.parent/"short_chain_section4_2026-09-14"
import sys
sys.path.insert(0, str(OLD/"stationary_density"))
from fp_tensorflow import relative_adjoint_residual  # noqa: E402


def values(model, state, batch=8192):
    return np.concatenate([
        model.log_ratio(tf.constant(state[s:s+batch], tf.float64)).numpy()
        for s in range(0, len(state), batch)
    ])


def residual(model, state, batch=256):
    return np.concatenate([
        relative_adjoint_residual(
            model, tf.constant(state[s:s+batch], tf.float64), 5.0, 5.0, 0.1
        ).numpy() for s in range(0, len(state), batch)
    ])


def main():
    choice = json.loads((HERE/"EQUILIBRIUM_ADMISSIBLE_SET.json").read_text())
    if choice["status"] != "FROZEN":
        raise SystemExit("equilibrium validation gate did not pass")
    target_rows, reference_rows = target_reference_rows("equilibrium", "fine", "test")
    target, _ = load_rows(target_rows); reference, _ = load_rows(reference_rows)
    target = target[fixed_indices(len(target), 100_000, "blind-eq-target")]
    reference = reference[fixed_indices(len(reference), 100_000, "blind-eq-reference")]
    fp_state = reference[fixed_indices(len(reference), 10_000, "blind-eq-fp")]
    rows = []
    for architecture in choice["admissible_architectures"]:
        for seed in (5201, 5202, 5203):
            directory = HERE/"equilibrium_candidates"/f"{architecture}_seed{seed}"
            config = json.loads((directory/"model_config.json").read_text())
            model = load_bundle(directory)
            ft = values(model, target); fq = values(model, reference)
            frozen_logz = float(config["validation_logz"])
            test_logz = float(logsumexp(fq)-math.log(len(fq)))
            w = np.exp(fq-frozen_logz)
            rr = residual(model, fp_state)
            normalized = ft-frozen_logz
            auc = float(roc_auc_score(
                np.r_[np.ones(len(ft)), np.zeros(len(fq))], np.r_[ft, fq]
            ))
            item = {
                "architecture": architecture, "seed": seed,
                "centered_log_ratio_rms": float(np.sqrt(np.mean(normalized*normalized))),
                "auc": auc, "validation_logz": frozen_logz,
                "blind_logz": test_logz,
                "ess_fraction": float(w.sum()**2/np.sum(w*w)/len(w)),
                "max_weight_share": float(w.max()/w.sum()),
                "fp_abs_median": float(np.median(np.abs(rr))),
                "fp_abs_p90": float(np.quantile(np.abs(rr), 0.90)),
                "finite": bool(np.isfinite(np.r_[ft, fq, rr, w]).all()),
            }
            item["pass"] = bool(
                item["centered_log_ratio_rms"] <= 0.05
                and 0.48 <= auc <= 0.52 and item["ess_fraction"] >= 0.50
                and item["fp_abs_median"] <= 0.10 and item["fp_abs_p90"] <= 0.50
                and item["finite"]
            )
            rows.append(item)
    result = {
        "protocol_version": "section4-density-ratio-v1",
        "blind_opened": True,
        "target_streams": [int(x["stream_id"]) for x in target_rows],
        "reference_streams": [int(x["stream_id"]) for x in reference_rows],
        "rows": rows,
        "pass": all(x["pass"] for x in rows),
    }
    (HERE/"EQUILIBRIUM_BLIND_GATE.json").write_text(
        json.dumps(result, indent=2, sort_keys=True)+"\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if result["pass"] else 2)


if __name__ == "__main__":
    tf.keras.backend.set_floatx("float64")
    main()

