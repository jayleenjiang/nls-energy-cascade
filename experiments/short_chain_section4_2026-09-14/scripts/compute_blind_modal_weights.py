#!/usr/bin/env python3
"""Complete the frozen blind modal-weight diagnostic without reselection."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg

from mode_common import (
    ROOT, N_ORIGINS, Aggregate, extended, feature_matrix, fit, load_state,
    nearest, origin_indices, split_paths,
)


OUT = ROOT / "modes"
SELECTION = OUT / "MODE_SELECTION.json"
BLIND_RESULTS = OUT / "BLIND_MODE_RESULTS.json"
MARKER = OUT / "MODAL_WEIGHT_READ_MARKER.json"
BOOTSTRAPS = 500
BOOTSTRAP_SEED = 2026091405
OBSERVABLES = (
    "I2", "I1_plus_I3", "cos_theta3", "cos_theta1_plus_cos_theta3",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def observables(state: np.ndarray) -> np.ndarray:
    return np.column_stack((
        state[:, 1],
        state[:, 0] + state[:, 2],
        np.cos(state[:, 4]),
        np.cos(state[:, 3]) + np.cos(state[:, 4]),
    ))


def contribution_weights(A: np.ndarray, b: np.ndarray, gmean: np.ndarray,
                         result, eigen_index: int) -> np.ndarray:
    # This is the frozen formula in the prior targeted EDMD implementation.
    centered = b - A[:, [0]] * gmean[None, :]
    coordinates = result.transform.T @ centered
    vectors = result.reduced_vectors
    amplitudes = linalg.inv(vectors, check_finite=False) @ coordinates
    inner = centered.T @ (result.transform @ vectors)
    contributions = inner.T * amplitudes
    reconstructed = np.sum(contributions, axis=0)
    rate = result.lam[eigen_index]
    if abs(rate.imag) > 1.0e-8:
        partner = int(np.argmin(np.abs(result.lam - np.conj(rate))))
        selected = contributions[eigen_index] + contributions[partner]
    else:
        partner = -1
        selected = contributions[eigen_index]
    weights = np.divide(
        selected, reconstructed,
        out=np.full_like(selected, np.nan, dtype=np.complex128),
        where=np.abs(reconstructed) > 0,
    )
    return weights, partner, reconstructed, selected


def main() -> None:
    if MARKER.exists():
        raise RuntimeError("modal-weight blind pass has already been recorded")
    selection = json.loads(SELECTION.read_text())
    selection_hash = sha256(SELECTION)
    blind_hash = sha256(BLIND_RESULTS)
    # Preflight every non-test dependency before the additional test read.
    for key, chosen in sorted(selection["selections"].items()):
        case, _ = key.split(":")
        if len(split_paths(case, "test")) != 16:
            raise RuntimeError(f"bad test split for {case}")
        if not (ROOT / "modes" / "train_validation_matrices" / f"{case}.npz").is_file():
            raise RuntimeError(f"missing frozen train matrix for {case}")
        if not (ROOT / chosen["frozen_fit"]).is_file():
            raise RuntimeError(f"missing frozen selected fit for {key}")
    MARKER.write_text(json.dumps({
        "status": "started",
        "selection_sha256": selection_hash,
        "blind_results_sha256": blind_hash,
        "additional_test_read": "modal weights only",
    }, indent=2, sort_keys=True) + "\n")

    summaries: list[dict] = []
    bootstrap_rows: list[dict] = []
    stream_rows: list[dict] = []
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    cases = sorted({key.split(":")[0] for key in selection["selections"]})
    for case in cases:
        chosen_for_case = {
            key.split(":")[1]: value
            for key, value in selection["selections"].items()
            if key.split(":")[0] == case
        }
        matrix = np.load(ROOT / "modes" / "train_validation_matrices" / f"{case}.npz")
        log_mean, log_sd = matrix["log_mean"], matrix["log_sd"]
        origins = origin_indices()
        test_paths = split_paths(case, "test")
        K = len(extended.DICT_INDICES["E3"])
        A_stream = np.empty((16, K, K), dtype=np.float64)
        b_stream = np.empty((16, K, len(OBSERVABLES)), dtype=np.float64)
        g_stream = np.empty((16, len(OBSERVABLES)), dtype=np.float64)
        for stream_order, path in enumerate(test_paths):
            state = load_state(path)[origins]
            X = feature_matrix(state, log_mean, log_sd)
            g = observables(state)
            A_stream[stream_order] = X.T @ X / N_ORIGINS
            b_stream[stream_order] = X.T @ g / N_ORIGINS
            g_stream[stream_order] = g.mean(axis=0)
            stream_id = int(path.stem.rsplit("_", 1)[-1])
            for q, observable in enumerate(OBSERVABLES):
                stream_rows.append({
                    "case": case, "stream_order": stream_order + 1,
                    "stream_id": stream_id, "observable": observable,
                    "mean": float(g_stream[stream_order, q]),
                })
            print(f"modal weights {case}: stream {stream_order+1}/16", flush=True)

        for kind, chosen in sorted(chosen_for_case.items()):
            indices = extended.DICT_INDICES[chosen["dictionary"]]
            training = Aggregate(
                matrix["train_A"], matrix["train_B"], matrix["train_C"],
                32 * N_ORIGINS, 32,
            )
            result = fit(training, chosen["dictionary"], chosen["tau"],
                         chosen["cutoff"], vectors=True)
            matched = nearest(complex(*chosen["lambda"]), result)
            if matched is None:
                raise RuntimeError(f"cannot recover frozen mode {case}:{kind}")
            eigen_index, _, rate = matched
            A = A_stream[:, indices][:, :, indices].mean(axis=0)
            b = b_stream[:, indices].mean(axis=0)
            gmean = g_stream.mean(axis=0)
            weights, partner, reconstructed, contribution = contribution_weights(
                A, b, gmean, result, eigen_index
            )
            bootstrap = np.empty((BOOTSTRAPS, len(OBSERVABLES)), dtype=np.float64)
            for replicate in range(BOOTSTRAPS):
                pick = rng.integers(0, 16, size=16)
                w, _, _, _ = contribution_weights(
                    A_stream[pick][:, indices][:, :, indices].mean(axis=0),
                    b_stream[pick][:, indices].mean(axis=0),
                    g_stream[pick].mean(axis=0), result, eigen_index,
                )
                bootstrap[replicate] = np.abs(w)
                for q, observable in enumerate(OBSERVABLES):
                    bootstrap_rows.append({
                        "case": case, "kind": kind,
                        "replicate": replicate, "observable": observable,
                        "weight_abs": float(bootstrap[replicate, q]),
                    })
            for q, observable in enumerate(OBSERVABLES):
                lo, hi = np.nanpercentile(bootstrap[:, q], (2.5, 97.5))
                summaries.append({
                    "case": case, "kind": kind,
                    "dictionary": chosen["dictionary"],
                    "cutoff": chosen["cutoff"], "tau": chosen["tau"],
                    "lambda_real": float(rate.real),
                    "lambda_imag": float(rate.imag),
                    "eigen_index": int(eigen_index),
                    "conjugate_index": int(partner),
                    "observable": observable,
                    "contribution_real": float(contribution[q].real),
                    "contribution_imag": float(contribution[q].imag),
                    "reconstructed_variance_real": float(reconstructed[q].real),
                    "reconstructed_variance_imag": float(reconstructed[q].imag),
                    "weight_real": float(weights[q].real),
                    "weight_imag": float(weights[q].imag),
                    "weight_abs": float(abs(weights[q])),
                    "weight_abs_ci_low": float(lo),
                    "weight_abs_ci_high": float(hi),
                    "bootstrap_replicates": BOOTSTRAPS,
                })

    def write_csv(path: Path, rows: list[dict]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)

    write_csv(OUT / "blind_modal_weights.csv", summaries)
    write_csv(OUT / "blind_modal_weights_bootstrap.csv", bootstrap_rows)
    write_csv(OUT / "blind_modal_weight_stream_means.csv", stream_rows)
    MARKER.write_text(json.dumps({
        "status": "complete",
        "selection_sha256": selection_hash,
        "blind_results_sha256": blind_hash,
        "additional_test_read": "modal weights only",
        "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "summary_sha256": sha256(OUT / "blind_modal_weights.csv"),
        "bootstrap_sha256": sha256(OUT / "blind_modal_weights_bootstrap.csv"),
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "rows": len(summaries)}, indent=2))


if __name__ == "__main__":
    main()
