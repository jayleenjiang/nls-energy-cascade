#!/usr/bin/env python3
"""Direct-sampling SCGF and importance-weight diagnostics.

The input is the audited ``NLS_entropy_ft.cpp`` block output.  Consecutive
base blocks are summed within each independent stream to form non-overlapping
windows.  For the medium entropy ``Sigma`` this script estimates

    psi_t(k) = t^{-1} log E[exp(-k Sigma_t)]

and the finite-time Gallavotti--Cohen residual

    D_t(k) = psi_t(k) - psi_t(1-k).

The exponential average can look precise while being dominated by one rare
sample.  We therefore report both sample-weight and independent-stream ESS,
maximum sample/stream weight fractions, and paired stream-bootstrap intervals.
The predeclared reliability gate is diagnostic rather than a theorem:

* sample-weight ESS >= 1000;
* independent-stream ESS >= 32;
* maximum sample weight fraction <= 0.01;
* maximum stream weight fraction <= 0.10.

No plus-four or other pseudocount enters the SCGF estimate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ESS_MIN = 1000.0
STREAM_ESS_MIN = 32.0
MAX_SAMPLE_FRACTION = 0.01
MAX_STREAM_FRACTION = 0.10


@dataclass
class TauEstimate:
    rows: list[dict[str, float | int | str]]
    psi_bootstrap: np.ndarray
    reliable: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--base-time", type=float, default=20.0)
    parser.add_argument("--max-multiple", type=int, default=10)
    parser.add_argument("--k-min", type=float, default=0.0)
    parser.add_argument("--k-max", type=float, default=1.0)
    parser.add_argument("--k-step", type=float, default=0.025)
    parser.add_argument("--bootstraps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def logsumexp(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    maximum = np.max(values, axis=axis, keepdims=True)
    shifted = np.exp(values - maximum)
    result = maximum + np.log(np.sum(shifted, axis=axis, keepdims=True))
    if axis is None:
        return np.asarray(result.squeeze())
    return np.squeeze(result, axis=axis)


def infer_n(path: Path) -> int:
    match = re.search(r"n(\d+)_blocks", path.name)
    if not match:
        raise ValueError(f"cannot infer chain length from {path}")
    return int(match.group(1))


def load_entropy_matrix(path: Path) -> np.ndarray:
    raw = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 1, 5))
    if raw.ndim != 2 or raw.shape[1] != 3 or raw.shape[0] == 0:
        raise ValueError(f"unexpected raw shape in {path}: {raw.shape}")
    if not np.all(np.isfinite(raw)):
        raise ValueError(f"nonfinite values in {path}")

    stream_ids = raw[:, 0].astype(np.int64)
    block_ids = raw[:, 1].astype(np.int64)
    streams, counts = np.unique(stream_ids, return_counts=True)
    if not np.array_equal(streams, np.arange(streams.size)):
        raise ValueError(f"stream IDs are not contiguous in {path}")
    if np.unique(counts).size != 1:
        raise ValueError(f"streams have unequal block counts in {path}")

    blocks_per_stream = int(counts[0])
    expected_streams = np.repeat(streams, blocks_per_stream)
    expected_blocks = np.tile(np.arange(blocks_per_stream), streams.size)
    if not np.array_equal(stream_ids, expected_streams):
        raise ValueError(f"rows are not grouped by stream in {path}")
    if not np.array_equal(block_ids, expected_blocks):
        raise ValueError(f"block IDs are not ordered within streams in {path}")
    return raw[:, 2].reshape(streams.size, blocks_per_stream)


def aggregate_nonoverlapping(matrix: np.ndarray, multiple: int) -> np.ndarray:
    usable = (matrix.shape[1] // multiple) * multiple
    if usable == 0:
        raise ValueError("aggregation multiple exceeds available blocks")
    return matrix[:, :usable].reshape(matrix.shape[0], -1, multiple).sum(axis=2)


def estimate_tau(
    entropy: np.ndarray,
    tau: float,
    k_values: np.ndarray,
    bootstraps: int,
    seed: int,
) -> TauEstimate:
    n_streams, groups_per_stream = entropy.shape
    sample_count = int(entropy.size)
    rng = np.random.default_rng(seed)
    bootstrap_counts = rng.multinomial(
        n_streams,
        np.full(n_streams, 1.0 / n_streams),
        size=bootstraps,
    )

    rows: list[dict[str, float | int | str]] = []
    psi_bootstrap = np.empty((bootstraps, k_values.size), dtype=float)
    reliable = np.zeros(k_values.size, dtype=bool)

    for index, k in enumerate(k_values):
        log_weights = -k * entropy
        stream_log_sums = logsumexp(log_weights, axis=1)
        total_log_sum = float(logsumexp(stream_log_sums))
        total_log_sum_sq = float(logsumexp(2.0 * log_weights))
        maximum_log_weight = float(np.max(log_weights))
        maximum_stream_log_weight = float(np.max(stream_log_sums))

        psi = (total_log_sum - math.log(sample_count)) / tau
        log_ess = 2.0 * total_log_sum - total_log_sum_sq
        sample_ess = min(float(sample_count), math.exp(min(log_ess, 700.0)))
        max_sample_fraction = math.exp(maximum_log_weight - total_log_sum)

        stream_log_sum_sq = float(logsumexp(2.0 * stream_log_sums))
        stream_log_ess = 2.0 * total_log_sum - stream_log_sum_sq
        stream_ess = min(float(n_streams), math.exp(min(stream_log_ess, 700.0)))
        max_stream_fraction = math.exp(
            maximum_stream_log_weight - total_log_sum
        )

        normalized_stream_sums = np.exp(
            stream_log_sums - maximum_stream_log_weight
        )
        bootstrap_sums = bootstrap_counts @ normalized_stream_sums
        if np.any(bootstrap_sums <= 0.0):
            raise RuntimeError("nonpositive bootstrap weight sum")
        psi_b = (
            maximum_stream_log_weight
            + np.log(bootstrap_sums)
            - math.log(sample_count)
        ) / tau
        psi_bootstrap[:, index] = psi_b
        ci_low, ci_high = np.quantile(psi_b, [0.025, 0.975])

        gate = (
            sample_ess >= ESS_MIN
            and stream_ess >= STREAM_ESS_MIN
            and max_sample_fraction <= MAX_SAMPLE_FRACTION
            and max_stream_fraction <= MAX_STREAM_FRACTION
        )
        reliable[index] = gate
        rows.append(
            {
                "tau": tau,
                "k": float(k),
                "psi": psi,
                "psi_ci_low": float(ci_low),
                "psi_ci_high": float(ci_high),
                "sample_count": sample_count,
                "n_streams": n_streams,
                "groups_per_stream": groups_per_stream,
                "sample_weight_ess": sample_ess,
                "sample_weight_relative_ess": sample_ess / sample_count,
                "stream_weight_ess": stream_ess,
                "max_sample_weight_fraction": max_sample_fraction,
                "max_stream_weight_fraction": max_stream_fraction,
                "direct_reliability_gate": int(gate),
            }
        )
    return TauEstimate(rows, psi_bootstrap, reliable)


def symmetry_rows(
    n: int,
    tau: float,
    k_values: np.ndarray,
    estimate: TauEstimate,
) -> list[dict[str, float | int]]:
    result: list[dict[str, float | int]] = []
    psi = np.asarray([float(row["psi"]) for row in estimate.rows])
    for index, k in enumerate(k_values):
        partner = int(np.argmin(np.abs(k_values - (1.0 - k))))
        if not np.isclose(k_values[partner], 1.0 - k, atol=1.0e-12):
            raise ValueError("k grid is not symmetric around 1/2")
        bootstrap_difference = (
            estimate.psi_bootstrap[:, index]
            - estimate.psi_bootstrap[:, partner]
        )
        ci_low, ci_high = np.quantile(bootstrap_difference, [0.025, 0.975])
        pair_gate = bool(estimate.reliable[index] and estimate.reliable[partner])
        result.append(
            {
                "n": n,
                "tau": tau,
                "k": float(k),
                "one_minus_k": float(k_values[partner]),
                "symmetry_residual": float(psi[index] - psi[partner]),
                "residual_ci_low": float(ci_low),
                "residual_ci_high": float(ci_high),
                "pair_reliability_gate": int(pair_gate),
                "ci_contains_zero": int(ci_low <= 0.0 <= ci_high),
            }
        )
    return result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize(
    n: int,
    tau: float,
    k_values: np.ndarray,
    estimate: TauEstimate,
    symmetry: list[dict[str, float | int]],
) -> dict[str, float | int | str]:
    reliable_indices = np.flatnonzero(estimate.reliable)
    reliable_k_min = (
        float(k_values[reliable_indices[0]]) if reliable_indices.size else math.nan
    )
    reliable_k_max = (
        float(k_values[reliable_indices[-1]]) if reliable_indices.size else math.nan
    )
    pair_rows = [
        row
        for row in symmetry
        if int(row["pair_reliability_gate"]) == 1
        and not math.isclose(float(row["k"]), 0.5)
    ]
    residuals = np.asarray(
        [float(row["symmetry_residual"]) for row in pair_rows], dtype=float
    )
    endpoint_index = int(np.argmin(np.abs(k_values - 1.0)))
    endpoint = estimate.rows[endpoint_index]
    return {
        "n": n,
        "tau": tau,
        "reliable_k_min": reliable_k_min,
        "reliable_k_max": reliable_k_max,
        "reliable_k_count": int(reliable_indices.size),
        "reliable_symmetry_pair_count": len(pair_rows),
        "reliable_residual_max_abs": (
            float(np.max(np.abs(residuals))) if residuals.size else math.nan
        ),
        "reliable_residual_rms": (
            float(np.sqrt(np.mean(residuals**2))) if residuals.size else math.nan
        ),
        "psi_k1": float(endpoint["psi"]),
        "psi_k1_ci_low": float(endpoint["psi_ci_low"]),
        "psi_k1_ci_high": float(endpoint["psi_ci_high"]),
        "k1_sample_weight_ess": float(endpoint["sample_weight_ess"]),
        "k1_stream_weight_ess": float(endpoint["stream_weight_ess"]),
        "k1_max_sample_weight_fraction": float(
            endpoint["max_sample_weight_fraction"]
        ),
        "k1_max_stream_weight_fraction": float(
            endpoint["max_stream_weight_fraction"]
        ),
        "k1_direct_reliability_gate": int(endpoint["direct_reliability_gate"]),
    }


def plot_n(
    output: Path,
    n: int,
    k_values: np.ndarray,
    tau_estimates: dict[float, TauEstimate],
    symmetry_by_tau: dict[float, list[dict[str, float | int]]],
) -> None:
    selected = [
        tau
        for tau in sorted(tau_estimates)
        if tau in {20.0, 40.0, 80.0, 120.0, 160.0, 200.0}
    ]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(selected)))
    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.7))

    for color, tau in zip(colors, selected):
        estimate = tau_estimates[tau]
        psi = np.asarray([float(row["psi"]) for row in estimate.rows])
        axes[0].plot(k_values, psi, color=color, label=rf"$t={tau:g}$")
        symmetry = symmetry_by_tau[tau]
        residual = np.asarray(
            [float(row["symmetry_residual"]) for row in symmetry]
        )
        gate = np.asarray(
            [bool(int(row["pair_reliability_gate"])) for row in symmetry]
        )
        axes[1].plot(k_values, residual, color=color, alpha=0.28)
        axes[1].plot(
            k_values[gate], residual[gate], "o-", ms=2.8, lw=1.1,
            color=color, label=rf"$t={tau:g}$",
        )

    axes[0].set_xlabel(r"tilt $k$")
    axes[0].set_ylabel(r"$\psi_t(k)$")
    axes[0].set_title("Direct SCGF")
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].grid(alpha=0.22)

    axes[1].axhline(0.0, color="black", lw=0.9, ls="--")
    axes[1].set_xlabel(r"tilt $k$")
    axes[1].set_ylabel(r"$\psi_t(k)-\psi_t(1-k)$")
    axes[1].set_title("GC residual (markers pass ESS gate)")
    axes[1].grid(alpha=0.22)

    taus = np.asarray(sorted(tau_estimates), dtype=float)
    log_ess = np.empty((taus.size, k_values.size), dtype=float)
    for row_index, tau in enumerate(taus):
        log_ess[row_index] = np.log10(
            [
                max(float(row["sample_weight_ess"]), 1.0)
                for row in tau_estimates[tau].rows
            ]
        )
    image = axes[2].imshow(
        log_ess,
        aspect="auto",
        origin="lower",
        extent=[k_values[0], k_values[-1], taus[0], taus[-1]],
        cmap="magma",
    )
    axes[2].set_xlabel(r"tilt $k$")
    axes[2].set_ylabel(r"averaging time $t$")
    axes[2].set_title(r"$\log_{10}$ sample-weight ESS")
    figure.colorbar(image, ax=axes[2], pad=0.02)

    figure.suptitle(rf"Medium-entropy direct-sampling pilot, $n={n}$")
    figure.tight_layout()
    figure.savefig(output / f"scgf_ess_n{n}.png", dpi=220)
    figure.savefig(output / f"scgf_ess_n{n}.pdf")
    plt.close(figure)


def self_test() -> None:
    values = np.asarray([[-1.0, 0.0], [1.0, 2.0]])
    estimate = estimate_tau(
        values, tau=1.0, k_values=np.asarray([0.0, 0.5, 1.0]),
        bootstraps=50, seed=17,
    )
    flattened = values.ravel()
    for row in estimate.rows:
        k = float(row["k"])
        weights = np.exp(-k * flattened)
        expected_psi = math.log(float(np.mean(weights)))
        expected_ess = float(weights.sum() ** 2 / np.square(weights).sum())
        if not math.isclose(float(row["psi"]), expected_psi, rel_tol=1e-13):
            raise AssertionError("SCGF algebra self-test failed")
        if not math.isclose(
            float(row["sample_weight_ess"]), expected_ess, rel_tol=1e-13
        ):
            raise AssertionError("ESS algebra self-test failed")
    print("SCGF/ESS self-test: PASS")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if not args.inputs or args.output_dir is None:
        raise SystemExit("inputs and --output-dir are required unless --self-test")
    if args.bootstraps < 20 or args.max_multiple < 1 or args.k_step <= 0.0:
        raise ValueError("invalid analysis configuration")

    count = int(round((args.k_max - args.k_min) / args.k_step))
    k_values = args.k_min + args.k_step * np.arange(count + 1)
    if not np.isclose(k_values[-1], args.k_max):
        raise ValueError("k range must be divisible by k step")
    if not np.allclose(k_values, 1.0 - k_values[::-1]):
        raise ValueError("k grid must be symmetric under k -> 1-k")

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    scgf_all: list[dict[str, object]] = []
    symmetry_all: list[dict[str, object]] = []
    summary_all: list[dict[str, object]] = []

    for path in args.inputs:
        n = infer_n(path)
        print(f"Loading n={n}: {path}", flush=True)
        base_entropy = load_entropy_matrix(path)
        tau_estimates: dict[float, TauEstimate] = {}
        symmetry_by_tau: dict[float, list[dict[str, float | int]]] = {}
        for multiple in range(1, args.max_multiple + 1):
            tau = args.base_time * multiple
            aggregated = aggregate_nonoverlapping(base_entropy, multiple)
            estimate = estimate_tau(
                aggregated,
                tau=tau,
                k_values=k_values,
                bootstraps=args.bootstraps,
                seed=args.seed + n * 1000 + multiple,
            )
            tau_estimates[tau] = estimate
            for row in estimate.rows:
                scgf_all.append({"source": str(path), "n": n, **row})
            symmetry = symmetry_rows(n, tau, k_values, estimate)
            symmetry_by_tau[tau] = symmetry
            symmetry_all.extend(symmetry)
            summary_all.append(
                summarize(n, tau, k_values, estimate, symmetry)
            )
            k1 = estimate.rows[-1]
            print(
                f"  t={tau:g}: k=1 ESS={float(k1['sample_weight_ess']):.2f}, "
                f"stream ESS={float(k1['stream_weight_ess']):.2f}, "
                f"gate={int(k1['direct_reliability_gate'])}",
                flush=True,
            )
        plot_n(output, n, k_values, tau_estimates, symmetry_by_tau)

    write_csv(output / "scgf_points.csv", scgf_all)
    write_csv(output / "gc_symmetry_residuals.csv", symmetry_all)
    write_csv(output / "scgf_reliability_summary.csv", summary_all)
    metadata = {
        "observable": "medium entropy production",
        "definition": "psi_t(k)=log(E[exp(-k Sigma_medium)])/t",
        "symmetry_test": "psi_t(k)=psi_t(1-k)",
        "inputs": [str(path) for path in args.inputs],
        "base_time": args.base_time,
        "max_multiple": args.max_multiple,
        "k_min": args.k_min,
        "k_max": args.k_max,
        "k_step": args.k_step,
        "bootstraps": args.bootstraps,
        "bootstrap_unit": "independent stream",
        "seed": args.seed,
        "reliability_gate": {
            "sample_weight_ess_min": ESS_MIN,
            "stream_weight_ess_min": STREAM_ESS_MIN,
            "max_sample_weight_fraction": MAX_SAMPLE_FRACTION,
            "max_stream_weight_fraction": MAX_STREAM_FRACTION,
        },
    }
    (output / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(f"Wrote SCGF pilot to {output}")


if __name__ == "__main__":
    main()
