#!/usr/bin/env python3
# analyze_lte_v2.py -- the frozen Sec. 3.2 analysis of fixed_reference_2026-09-18 (analyze_lte.py),
# applied unchanged to the nls_cartesian v2 histograms in raw/cartesian_v2.  Differences from the
# original: integrator list, samples per batch (200000: v2 samples at 0.01, 0.02, ..., 2000, without the
# extra sample at measurement time zero), output names, and an extra file with the batch-summed
# histograms (analysis/aggregated_hist_v2.npz) for the residual-mesh figure.

from __future__ import annotations

import csv
import hashlib
import json
import math
import struct
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
NB = 80
N_BOOT = 1000
BOOT_SEED = 2026091800
MIN_COUNT = 50
CORE_MAX = 2.5
EXPECTED_TOTAL = 16 * 200_000
EXPECTED_BATCHES = 32


def read_matrix() -> list[dict[str, str]]:
    with (ROOT / "RUN_MATRIX.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_histogram(path: Path) -> tuple[np.ndarray, dict[str, int]]:
    with path.open("rb") as handle:
        magic = handle.read(8)
        nb, site = struct.unpack("<Ii", handle.read(8))
        total, overflow = struct.unpack("<QQ", handle.read(16))
        bins = np.fromfile(handle, dtype="<u4", count=nb**3)
    if magic != b"LTEH80V1" or nb != NB or bins.size != NB**3:
        raise RuntimeError(f"invalid histogram {path}")
    if total != EXPECTED_TOTAL or int(bins.sum()) != total - overflow:
        raise RuntimeError(f"histogram count mismatch {path}")
    return bins.reshape((NB, NB, NB)), {
        "site": site,
        "total": total,
        "overflow": overflow,
    }


def prefixes(temp_label: str, n: int, integrator: str) -> list[Path]:
    result: list[Path] = []
    for replicate in range(2):
        for batch in range(16):
            name = f"{temp_label}_n{n}_rep{replicate}_batch{batch:02d}"
            result.append(ROOT / "raw" / integrator / name)
    if len(result) != EXPECTED_BATCHES:
        raise AssertionError("internal batch count mismatch")
    return result


def load_batches(temp_label: str, n: int, site: int, integrator: str,
                 mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    arrays = []
    overflows = []
    for prefix in prefixes(temp_label, n, integrator):
        histogram, metadata = read_histogram(Path(f"{prefix}_j{site}.histbin"))
        if metadata["site"] != site:
            raise RuntimeError(f"site mismatch for {prefix}")
        arrays.append(histogram[mask])
        overflows.append(metadata["overflow"] / metadata["total"])
    return np.asarray(arrays, dtype=np.float64), np.asarray(overflows)


def normalize(values: np.ndarray) -> np.ndarray:
    total = values.sum()
    if not total > 0:
        raise RuntimeError("cannot normalize empty histogram")
    return values / total


def js_distance(left: np.ndarray, right: np.ndarray) -> float:
    p = normalize(left)
    q = normalize(right)
    middle = 0.5 * (p + q)
    left_mask = p > 0
    right_mask = q > 0
    divergence = 0.5 * np.sum(p[left_mask] * np.log(p[left_mask] / middle[left_mask]))
    divergence += 0.5 * np.sum(q[right_mask] * np.log(q[right_mask] / middle[right_mask]))
    return float(np.sqrt(max(divergence, 0.0)))


def weighted_fit(q_counts: np.ndarray, p_counts: np.ndarray,
                 reliable: np.ndarray) -> tuple[float, float, float]:
    q = normalize(q_counts)
    p = normalize(p_counts)
    x = np.log(p[reliable])
    y = np.log(q[reliable])
    weights = q_counts[reliable] + p_counts[reliable]
    design = np.column_stack((np.ones_like(x), x))
    root_w = np.sqrt(weights)
    coefficients, *_ = np.linalg.lstsq(
        design * root_w[:, None], y * root_w, rcond=None
    )
    intercept, slope = coefficients
    fitted = intercept + slope * x
    weighted_mean = np.average(y, weights=weights)
    residual = np.sum(weights * (y - fitted) ** 2)
    total = np.sum(weights * (y - weighted_mean) ** 2)
    r_squared = 1.0 - residual / total if total > 0 else math.nan
    return float(slope), float(intercept), float(r_squared)


def percentile_interval(values: np.ndarray) -> tuple[float, float]:
    return tuple(np.quantile(values, [0.025, 0.975]).tolist())


def bootstrap_metrics(q_batches: np.ndarray, p_batches: np.ndarray,
                      reliable: np.ndarray,
                      q_weights: np.ndarray,
                      p_weights: np.ndarray) -> dict[str, np.ndarray]:
    slopes = np.empty(N_BOOT)
    r_squared = np.empty(N_BOOT)
    distances = np.empty(N_BOOT)
    chunk = 20
    for start in range(0, N_BOOT, chunk):
        stop = min(N_BOOT, start + chunk)
        q_chunk = q_weights[start:stop] @ q_batches
        p_chunk = p_weights[start:stop] @ p_batches
        for offset in range(stop - start):
            index = start + offset
            slopes[index], _, r_squared[index] = weighted_fit(
                q_chunk[offset], p_chunk[offset], reliable
            )
            distances[index] = js_distance(q_chunk[offset], p_chunk[offset])
    return {"slope": slopes, "r2": r_squared, "js": distances}


def theta_mirror_distance(batch_histograms: np.ndarray) -> float:
    theta = batch_histograms.sum(axis=0).reshape(-1, NB).sum(axis=0)
    return js_distance(theta, theta[::-1])


def bootstrap_theta_floor(batch_histograms: np.ndarray,
                          weights: np.ndarray) -> np.ndarray:
    theta_batches = batch_histograms.reshape(EXPECTED_BATCHES, -1, NB).sum(axis=1)
    result = np.empty(N_BOOT)
    for start in range(0, N_BOOT, 100):
        stop = min(N_BOOT, start + 100)
        theta = weights[start:stop] @ theta_batches
        for offset in range(stop - start):
            result[start + offset] = js_distance(theta[offset], theta[offset, ::-1])
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    analysis = ROOT / "analysis"
    analysis.mkdir(exist_ok=True)
    edges = np.linspace(0.0, 4.0, NB + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    core2 = (centers[:, None] < CORE_MAX) & (centers[None, :] < CORE_MAX)
    core3 = np.broadcast_to(core2[:, :, None], (NB, NB, NB))

    rng = np.random.default_rng(BOOT_SEED)
    q_weights = rng.multinomial(
        EXPECTED_BATCHES,
        np.full(EXPECTED_BATCHES, 1.0 / EXPECTED_BATCHES),
        size=N_BOOT,
    ).astype(np.float64)
    p_weights = rng.multinomial(
        EXPECTED_BATCHES,
        np.full(EXPECTED_BATCHES, 1.0 / EXPECTED_BATCHES),
        size=N_BOOT,
    ).astype(np.float64)

    rows: list[dict[str, object]] = []
    aggregated: dict[str, np.ndarray] = {}
    for integrator in ("cartesian_v2",):
        for n in (25, 50, 100):
            sites = [int(0.24 * n), int(0.48 * n), int(0.72 * n)]
            for site in sites:
                q_batches, q_overflow = load_batches(
                    "NESS_T10_T5", n, site, integrator, core3
                )
                p_batches, p_overflow = load_batches(
                    "EQ_T7p5", n, site, integrator, core3
                )
                q_full = q_batches.sum(axis=0)
                p_full = p_batches.sum(axis=0)
                reliable = (q_full >= MIN_COUNT) & (p_full >= MIN_COUNT)
                if int(reliable.sum()) < 2:
                    raise RuntimeError(
                        f"too few reliable bins for {integrator}, n={n}, site={site}"
                    )
                slope, intercept, r_squared = weighted_fit(
                    q_full, p_full, reliable
                )
                distance = js_distance(q_full, p_full)
                bootstrap = bootstrap_metrics(
                    q_batches, p_batches, reliable, q_weights, p_weights
                )

                eq_batches_full, _ = load_batches(
                    "EQ_T7p5", n, site, integrator,
                    np.ones((NB, NB, NB), dtype=bool),
                )
                ness_batches_full, _ = load_batches(
                    "NESS_T10_T5", n, site, integrator,
                    np.ones((NB, NB, NB), dtype=bool),
                )
                aggregated[f"NESS_n{n}_j{site}"] = ness_batches_full.sum(axis=0).reshape(NB, NB, NB).astype(np.uint32)
                aggregated[f"EQ_n{n}_j{site}"] = eq_batches_full.sum(axis=0).reshape(NB, NB, NB).astype(np.uint32)
                mirror = theta_mirror_distance(eq_batches_full)
                mirror_boot = bootstrap_theta_floor(eq_batches_full, p_weights)

                slope_lo, slope_hi = percentile_interval(bootstrap["slope"])
                r2_lo, r2_hi = percentile_interval(bootstrap["r2"])
                js_lo, js_hi = percentile_interval(bootstrap["js"])
                mirror_lo, mirror_hi = percentile_interval(mirror_boot)
                rows.append(
                    {
                        "integrator": integrator,
                        "n": n,
                        "site_left_0based": site,
                        "position_over_n": site / n,
                        "slope_x": slope,
                        "slope_ci_low": slope_lo,
                        "slope_ci_high": slope_hi,
                        "intercept": intercept,
                        "weighted_r2": r_squared,
                        "r2_ci_low": r2_lo,
                        "r2_ci_high": r2_hi,
                        "js_distance": distance,
                        "js_ci_low": js_lo,
                        "js_ci_high": js_hi,
                        "reliable_bins": int(reliable.sum()),
                        "ness_core_mass_fraction":
                            q_full.sum() / (EXPECTED_BATCHES * EXPECTED_TOTAL),
                        "eq_core_mass_fraction":
                            p_full.sum() / (EXPECTED_BATCHES * EXPECTED_TOTAL),
                        "ness_max_batch_overflow_fraction": float(q_overflow.max()),
                        "eq_max_batch_overflow_fraction": float(p_overflow.max()),
                        "eq_theta_mirror_js": mirror,
                        "eq_theta_mirror_js_ci_low": mirror_lo,
                        "eq_theta_mirror_js_ci_high": mirror_hi,
                    }
                )
                print(
                    f"{integrator} n={n} site={site}: x={slope:.6f}, "
                    f"R2={r_squared:.6f}, JS={distance:.6f}",
                    flush=True,
                )

    np.savez_compressed(analysis / "aggregated_hist_v2.npz", **aggregated)
    output = analysis / "lte_summary_v2.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "analysis_source_sha256": sha256(Path(__file__)),
        "bootstrap_replicates": N_BOOT,
        "bootstrap_seed": BOOT_SEED,
        "bootstrap_unit": "independent 16-trajectory batch",
        "min_count_per_aggregate_side": MIN_COUNT,
        "core_action_max": CORE_MAX,
        "histogram_bins_per_axis": NB,
        "summary_rows": len(rows),
    }
    (analysis / "analysis_manifest_v2.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
