#!/usr/bin/env python3
"""Recompute Section 4.2 diagnostics from frozen densities and V6 trajectories."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[2]
DENSITY_ROOT = ROOT / "experiments/short_chain_density_ratio_constrained_2026-09-14"
HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"
FIGURES = HERE / "figures"
SCRIPTS = DENSITY_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from final_ness_density import FinalNESSDensity, energy  # noqa: E402
from v6_common import evaluation_roles, load_rows  # noqa: E402


BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260915042
ACTION_QUANTILES = (0.001, 0.999)
ANGLE_BINS = 60
LOCK_BINS = 72
GAMMA = 0.1
THETA0 = math.pi - math.asin(
    (math.sqrt(4.0 * GAMMA * GAMMA + 3.0) - GAMMA)
    / (2.0 * (1.0 + GAMMA * GAMMA))
)
THETA_LOW = -2.0 * math.pi / 3.0
LOCK_HALF_WIDTH = math.pi / 6.0
NEW_SEEDS = (8201, 8202, 8203)

OLD_DRIVEN_MODEL = ROOT / "KDE/4:15_NN/h5_files/final.keras"
OLD_EQUILIBRIUM_MODEL = ROOT / "KDE/4:15_NN/h5_files_eq/final.keras"
OLD_DRIVEN_DENSITY = ROOT / "KDE/4:15_NN/NLS_FP_density.txt"
OLD_EQUILIBRIUM_DENSITY = ROOT / "KDE/4:15_NN/eq/NLS_FP_density.txt"


@tf.keras.utils.register_keras_serializable()
def periodic_encode(x: tf.Tensor) -> tf.Tensor:
    return tf.concat(
        [x[:, :3], tf.cos(x[:, 3:4]), tf.sin(x[:, 3:4]),
         tf.cos(x[:, 4:5]), tf.sin(x[:, 4:5])], axis=-1
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def derived_seed(label: str) -> int:
    value = hashlib.sha256(f"{BOOTSTRAP_SEED}|{label}".encode()).hexdigest()
    return int(value[:16], 16)


def wrap(theta):
    return (np.asarray(theta) + np.pi) % (2.0 * np.pi) - np.pi


def angular_distance(theta, reference):
    return np.abs(wrap(np.asarray(theta) - reference))


def ci(values):
    values = np.asarray(values, np.float64)
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def stream_draws(n_streams: int, label: str) -> np.ndarray:
    rng = np.random.default_rng(derived_seed(label))
    return rng.multinomial(
        n_streams, np.full(n_streams, 1.0 / n_streams),
        size=BOOTSTRAP_REPLICATES,
    ).astype(np.float64)


def split_like(array: np.ndarray, streams: list[np.ndarray]) -> list[np.ndarray]:
    lengths = [len(stream) for stream in streams]
    cuts = np.cumsum(lengths)[:-1]
    return list(np.split(array, cuts))


@dataclass(frozen=True)
class Support:
    endpoint_low: float
    endpoint_high: float
    middle_low: float
    middle_high: float
    energy_low_mid: float
    energy_mid_high: float

    def mask(self, state: np.ndarray) -> np.ndarray:
        return (
            (state[:, 0] >= self.endpoint_low)
            & (state[:, 0] <= self.endpoint_high)
            & (state[:, 2] >= self.endpoint_low)
            & (state[:, 2] <= self.endpoint_high)
            & (state[:, 1] >= self.middle_low)
            & (state[:, 1] <= self.middle_high)
        )

    def sector(self, state: np.ndarray) -> np.ndarray:
        e = energy(state)
        return np.where(e < self.energy_low_mid, 0,
                        np.where(e < self.energy_mid_high, 1, 2))


def freeze_support(validation_target: np.ndarray) -> Support:
    pooled = np.concatenate((validation_target[:, 0], validation_target[:, 2]))
    endpoint_low, endpoint_high = np.quantile(pooled, ACTION_QUANTILES)
    middle_low, middle_high = np.quantile(validation_target[:, 1], ACTION_QUANTILES)
    preliminary = Support(
        float(endpoint_low), float(endpoint_high),
        float(middle_low), float(middle_high), 0.0, 0.0,
    )
    mask = preliminary.mask(validation_target)
    low_mid, mid_high = np.quantile(energy(validation_target[mask]), (1.0 / 3.0, 2.0 / 3.0))
    return Support(
        float(endpoint_low), float(endpoint_high),
        float(middle_low), float(middle_high),
        float(low_mid), float(mid_high),
    )


def old_log_normalization(path: Path) -> tuple[float, float]:
    density = np.loadtxt(path, dtype=np.float32)
    logs = np.log(density[density > 0])
    return float(logs.mean()), float(logs.std())


def predict_old_log_density(model, state: np.ndarray, mu: float, sigma: float,
                            batch: int = 32768) -> np.ndarray:
    result = []
    for start in range(0, len(state), batch):
        raw = model.predict(
            state[start:start + batch].astype(np.float32),
            batch_size=batch, verbose=0,
        ).reshape(-1)
        result.append(sigma * raw.astype(np.float64) + mu)
    return np.concatenate(result)


def scaled_weights(log_weight: np.ndarray) -> np.ndarray:
    log_weight = np.asarray(log_weight, np.float64)
    if not np.isfinite(log_weight).all():
        raise RuntimeError("non-finite log importance weight")
    return np.exp(log_weight - np.max(log_weight))


def effective_sample_fraction(weights: np.ndarray) -> float:
    weights = np.asarray(weights, np.float64)
    return float(weights.sum() ** 2 / (len(weights) * np.square(weights).sum()))


def probability_with_ci(stream_num: np.ndarray, stream_den: np.ndarray, label: str):
    estimate = float(stream_num.sum() / stream_den.sum())
    draws = stream_draws(len(stream_num), label)
    values = (draws @ stream_num) / (draws @ stream_den)
    low, high = ci(values)
    return estimate, low, high, values


def support_rows(case: str, support: Support, target_streams: list[np.ndarray],
                 reference_streams: list[np.ndarray], weight_sources: dict[str, list[np.ndarray]]):
    rows = []
    target_num = np.asarray([support.mask(s).sum() for s in target_streams], np.float64)
    target_den = np.asarray([len(s) for s in target_streams], np.float64)
    est, lo, hi, _ = probability_with_ci(target_num, target_den, f"support|{case}|direct")
    rows.append({
        "case": case, "source": "direct_trajectory", "mass": est,
        "ci_low": lo, "ci_high": hi, "ess_fraction_full": 1.0,
        "endpoint_low": support.endpoint_low, "endpoint_high": support.endpoint_high,
        "middle_low": support.middle_low, "middle_high": support.middle_high,
        "energy_low_mid": support.energy_low_mid,
        "energy_mid_high": support.energy_mid_high,
    })
    for source, streams in weight_sources.items():
        num = np.asarray([(w * support.mask(s)).sum() for s, w in zip(reference_streams, streams)])
        den = np.asarray([w.sum() for w in streams])
        est, lo, hi, _ = probability_with_ci(num, den, f"support|{case}|{source}")
        flat = np.concatenate(streams)
        rows.append({
            "case": case, "source": source, "mass": est,
            "ci_low": lo, "ci_high": hi,
            "ess_fraction_full": effective_sample_fraction(flat),
            "endpoint_low": support.endpoint_low, "endpoint_high": support.endpoint_high,
            "middle_low": support.middle_low, "middle_high": support.middle_high,
            "energy_low_mid": support.energy_low_mid,
            "energy_mid_high": support.energy_mid_high,
        })
    return rows


def exchange(state: np.ndarray) -> np.ndarray:
    return state[:, [2, 1, 0, 4, 3]]


def relative_asymmetry(log_a: np.ndarray, log_b: np.ndarray) -> np.ndarray:
    return np.abs(np.tanh(0.5 * (np.asarray(log_a) - np.asarray(log_b))))


def asymmetry_summary(case: str, source: str, values: np.ndarray) -> dict:
    return {
        "case": case, "source": source, "sample_count": len(values),
        "mean": float(np.mean(values)), "median": float(np.median(values)),
        "p90": float(np.quantile(values, 0.90)), "maximum": float(np.max(values)),
        "seed_min_mean": "", "seed_max_mean": "",
        "seed_min_median": "", "seed_max_median": "",
        "seed_min_p90": "", "seed_max_p90": "",
        "seed_min_maximum": "", "seed_max_maximum": "",
    }


def add_seed_ranges(rows: list[dict], case: str, ensemble_source: str,
                    seed_sources: list[str]) -> None:
    ensemble = next(row for row in rows if row["case"] == case and row["source"] == ensemble_source)
    selected = [row for row in rows if row["case"] == case and row["source"] in seed_sources]
    for metric in ("mean", "median", "p90", "maximum"):
        values = [float(row[metric]) for row in selected]
        ensemble[f"seed_min_{metric}"] = min(values)
        ensemble[f"seed_max_{metric}"] = max(values)


def angular_hist_streams(streams: list[np.ndarray], support: Support,
                         weights: list[np.ndarray] | None) -> np.ndarray:
    edges = np.linspace(-np.pi, np.pi, ANGLE_BINS + 1)
    output = []
    for index, state in enumerate(streams):
        mask = support.mask(state)
        weight = None if weights is None else weights[index][mask]
        hist, _, _ = np.histogram2d(
            state[mask, 3], state[mask, 4], bins=(edges, edges), weights=weight,
        )
        output.append(hist)
    return np.asarray(output, np.float64)


def normalize_hist(hist: np.ndarray) -> np.ndarray:
    total = hist.sum(axis=(-2, -1), keepdims=True)
    return hist / total


def tv_from_hist(hist: np.ndarray) -> np.ndarray:
    prob = normalize_hist(hist)
    return 0.5 * np.abs(prob - np.swapaxes(prob, -2, -1)).sum(axis=(-2, -1))


def angular_tv_row(case: str, source: str, per_stream_hist: np.ndarray):
    aggregate = per_stream_hist.sum(axis=0)
    estimate = float(tv_from_hist(aggregate[None, ...])[0])
    draws = stream_draws(len(per_stream_hist), f"angular-tv|{case}|{source}")
    boot_hist = (draws @ per_stream_hist.reshape(len(per_stream_hist), -1)).reshape(
        BOOTSTRAP_REPLICATES, ANGLE_BINS, ANGLE_BINS
    )
    values = tv_from_hist(boot_hist)
    low, high = ci(values)
    return {
        "case": case, "source": source, "tv_exchange": estimate,
        "ci_low": low, "ci_high": high,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed_min": "", "seed_max": "",
    }, aggregate / aggregate.sum(), values


def phase_stream_summaries(streams: list[np.ndarray], support: Support,
                           weights: list[np.ndarray] | None):
    denom = np.zeros((len(streams), 3))
    cos_sum = np.zeros((len(streams), 3, 2))
    sin_sum = np.zeros((len(streams), 3, 2))
    hist = np.zeros((len(streams), 3, 2, LOCK_BINS))
    high = np.zeros((len(streams), 3, 2))
    low3 = np.zeros((len(streams), 3))
    joint = np.zeros((len(streams), 3))
    edges = np.linspace(-np.pi, np.pi, LOCK_BINS + 1)
    for s_index, state in enumerate(streams):
        mask = support.mask(state)
        state = state[mask]
        weight = np.ones(len(state)) if weights is None else weights[s_index][mask]
        sector = support.sector(state)
        for q in range(3):
            selected = sector == q
            w = weight[selected]
            x = state[selected]
            denom[s_index, q] = w.sum()
            for angle in range(2):
                theta = x[:, 3 + angle]
                cos_sum[s_index, q, angle] = np.sum(w * np.cos(theta))
                sin_sum[s_index, q, angle] = np.sum(w * np.sin(theta))
                hist[s_index, q, angle], _ = np.histogram(theta, bins=edges, weights=w)
                high[s_index, q, angle] = np.sum(w * (angular_distance(theta, THETA0) <= LOCK_HALF_WIDTH))
            low3[s_index, q] = np.sum(
                w * (angular_distance(x[:, 4], THETA_LOW) <= LOCK_HALF_WIDTH)
            )
            joint[s_index, q] = np.sum(
                w * (angular_distance(x[:, 3], THETA0) <= LOCK_HALF_WIDTH)
                * (angular_distance(x[:, 4], THETA0) <= LOCK_HALF_WIDTH)
            )
    return {
        "denom": denom, "cos": cos_sum, "sin": sin_sum, "hist": hist,
        "high": high, "low3": low3, "joint": joint,
    }


def circular_boot_interval(point: float, draws: np.ndarray):
    difference = wrap(draws - point)
    lo, hi = ci(difference)
    return float(wrap(point + lo)), float(wrap(point + hi))


def phase_rows(case: str, source: str, summary: dict):
    rows = []
    sectors = ("low", "middle", "high")
    edges = np.linspace(-np.pi, np.pi, LOCK_BINS + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    draws = stream_draws(len(summary["denom"]), f"phase|{case}|{source}")
    denom_b = draws @ summary["denom"]
    for q, sector in enumerate(sectors):
        for a, angle in enumerate(("theta1", "theta3")):
            c = summary["cos"][:, q, a].sum()
            s = summary["sin"][:, q, a].sum()
            point_mean = float(math.atan2(s, c))
            point_r = float(math.hypot(c, s) / summary["denom"][:, q].sum())
            cb = draws @ summary["cos"][:, q, a]
            sb = draws @ summary["sin"][:, q, a]
            mean_b = np.arctan2(sb, cb)
            r_b = np.hypot(cb, sb) / denom_b[:, q]
            mean_lo, mean_hi = circular_boot_interval(point_mean, mean_b)
            r_lo, r_hi = ci(r_b)
            h = summary["hist"][:, q, a].sum(axis=0)
            point_mode = float(centers[int(np.argmax(h))])
            hb = draws @ summary["hist"][:, q, a]
            mode_b = centers[np.argmax(hb, axis=1)]
            mode_lo, mode_hi = circular_boot_interval(point_mode, mode_b)
            p = float(summary["high"][:, q, a].sum() / summary["denom"][:, q].sum())
            pb = (draws @ summary["high"][:, q, a]) / denom_b[:, q]
            p_lo, p_hi = ci(pb)
            for metric, estimate, low, high in (
                ("circular_mean", point_mean, mean_lo, mean_hi),
                ("circular_resultant", point_r, r_lo, r_hi),
                ("histogram_mode", point_mode, mode_lo, mode_hi),
                ("prob_within_pi_over_6_of_theta0", p, p_lo, p_hi),
            ):
                rows.append({
                    "case": case, "source": source, "sector": sector,
                    "angle": angle, "metric": metric, "estimate": estimate,
                    "ci_low": low, "ci_high": high,
                    "seed_min": "", "seed_max": "",
                })
            if angle == "theta3":
                p_low = float(summary["low3"][:, q].sum() / summary["denom"][:, q].sum())
                p_low_b = (draws @ summary["low3"][:, q]) / denom_b[:, q]
                low, high = ci(p_low_b)
                rows.append({
                    "case": case, "source": source, "sector": sector,
                    "angle": angle, "metric": "prob_within_pi_over_6_of_minus_2pi_over_3",
                    "estimate": p_low, "ci_low": low, "ci_high": high,
                    "seed_min": "", "seed_max": "",
                })
        p_joint = float(summary["joint"][:, q].sum() / summary["denom"][:, q].sum())
        p_joint_b = (draws @ summary["joint"][:, q]) / denom_b[:, q]
        low, high = ci(p_joint_b)
        rows.append({
            "case": case, "source": source, "sector": sector,
            "angle": "joint", "metric": "prob_both_within_pi_over_6_of_theta0",
            "estimate": p_joint, "ci_low": low, "ci_high": high,
            "seed_min": "", "seed_max": "",
        })
    return rows


def phase_hist_rows(case: str, source: str, summary: dict):
    rows = []
    edges = np.linspace(-np.pi, np.pi, LOCK_BINS + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    for q, sector in enumerate(("low", "middle", "high")):
        for a, angle in enumerate(("theta1", "theta3")):
            values = summary["hist"][:, q, a].sum(axis=0)
            values = values / (values.sum() * width)
            for b, center in enumerate(centers):
                rows.append({
                    "case": case, "source": source, "sector": sector,
                    "angle": angle, "bin": b, "theta": center,
                    "density": values[b],
                })
    return rows


def add_phase_seed_ranges(rows: list[dict], case: str):
    seed_names = [f"new_seed_{seed}" for seed in NEW_SEEDS]
    for row in rows:
        if row["case"] != case or row["source"] != "new_ensemble":
            continue
        selected = [candidate for candidate in rows
                    if candidate["case"] == case
                    and candidate["source"] in seed_names
                    and candidate["sector"] == row["sector"]
                    and candidate["angle"] == row["angle"]
                    and candidate["metric"] == row["metric"]]
        values = [float(candidate["estimate"]) for candidate in selected]
        row["seed_min"] = min(values)
        row["seed_max"] = max(values)


def current_stream_summaries(streams: list[np.ndarray], support: Support,
                             weights: list[np.ndarray] | None):
    sums = np.zeros((len(streams), 3))
    denom = np.zeros(len(streams))
    for index, state in enumerate(streams):
        mask = support.mask(state)
        x = state[mask]
        w = np.ones(len(x)) if weights is None else weights[index][mask]
        j12 = x[:, 0] * x[:, 1] * np.sin(x[:, 3])
        j23 = x[:, 1] * x[:, 2] * np.sin(x[:, 4])
        sums[index, 0] = np.sum(w * j12)
        sums[index, 1] = np.sum(w * j23)
        sums[index, 2] = np.sum(w * (j12 + j23))
        denom[index] = w.sum()
    return sums, denom


def current_rows(case: str, source: str, sums: np.ndarray, denom: np.ndarray):
    labels = ("J12", "J23", "J12_plus_J23")
    point = sums.sum(axis=0) / denom.sum()
    draws = stream_draws(len(denom), f"current|{case}|{source}")
    boot = (draws @ sums) / (draws @ denom)[:, None]
    rows = []
    for index, observable in enumerate(labels):
        low, high = ci(boot[:, index])
        rows.append({
            "case": case, "source": source, "observable": observable,
            "estimate": point[index], "ci_low": low, "ci_high": high,
            "bootstrap_se": float(np.std(boot[:, index], ddof=1)),
            "seed_min": "", "seed_max": "",
        })
    return rows, boot


def add_current_seed_ranges(rows: list[dict], case: str):
    seed_names = [f"new_seed_{seed}" for seed in NEW_SEEDS]
    for row in rows:
        if row["case"] != case or row["source"] != "new_ensemble":
            continue
        values = [float(candidate["estimate"]) for candidate in rows
                  if candidate["case"] == case and candidate["source"] in seed_names
                  and candidate["observable"] == row["observable"]]
        row["seed_min"] = min(values)
        row["seed_max"] = max(values)


def interval_overlap(a_low, a_high, b_low, b_high):
    return max(float(a_low), float(b_low)) <= min(float(a_high), float(b_high))


def plot_phase(hist_rows: list[dict]):
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.1), sharex=True)
    for q, sector in enumerate(("low", "middle", "high")):
        for a, angle in enumerate(("theta1", "theta3")):
            ax = axes[a, q]
            for source, color, label, style in (
                ("direct_trajectory", "black", "held-out trajectory", "-"),
                ("new_ensemble", "#0072B2", "validated density", "-"),
                ("legacy_density", "#D55E00", "legacy density", "--"),
            ):
                selected = [r for r in hist_rows if r["case"] == "driven"
                            and r["source"] == source and r["sector"] == sector
                            and r["angle"] == angle]
                ax.plot([r["theta"] for r in selected], [r["density"] for r in selected],
                        color=color, ls=style, lw=1.35, label=label)
            ax.axvline(THETA0, color="#009E73", lw=1.0, ls=":", label=r"$\theta_0$")
            if angle == "theta3":
                ax.axvline(THETA_LOW, color="#CC79A7", lw=1.0, ls=":",
                           label=r"$-2\pi/3$")
            angle_label = r"$\theta_1$" if angle == "theta1" else r"$\theta_3$"
            ax.set_title(f"{sector} energy, {angle_label}")
            ax.set_xlim(-np.pi, np.pi)
            ax.set_xticks((-np.pi, 0, np.pi), (r"$-\pi$", "0", r"$\pi$"))
            ax.grid(alpha=0.18)
    axes[0, 0].set_ylabel("conditional density")
    axes[1, 0].set_ylabel("conditional density")
    axes[1, 1].set_xlabel("angle")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=4, loc="upper center")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIGURES / "phase_locking_rebuilt.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "phase_locking_rebuilt.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_asymmetry(histograms: dict[tuple[str, str], np.ndarray]):
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.5), constrained_layout=True)
    extent = (-np.pi, np.pi, -np.pi, np.pi)
    for row, case in enumerate(("driven", "equilibrium")):
        direct = histograms[(case, "direct_trajectory")]
        new = histograms[(case, "new_ensemble")]
        denom = new + new.T
        asym = np.divide(new - new.T, denom, out=np.zeros_like(new), where=denom > 0)
        vmax = max(direct.max(), new.max())
        for col, data, title in (
            (0, direct.T, "held-out angular marginal"),
            (1, new.T, "density angular marginal"),
        ):
            image = axes[row, col].imshow(data, origin="lower", extent=extent,
                                          cmap="viridis", vmin=0, vmax=vmax,
                                          aspect="equal")
            axes[row, col].set_title(f"{case}: {title}")
            fig.colorbar(image, ax=axes[row, col], shrink=0.75)
        limit = max(float(np.max(np.abs(asym))), 1e-6)
        image = axes[row, 2].imshow(asym.T, origin="lower", extent=extent,
                                    cmap="coolwarm", vmin=-limit, vmax=limit,
                                    aspect="equal")
        axes[row, 2].set_title(f"{case}: exchange asymmetry")
        fig.colorbar(image, ax=axes[row, 2], shrink=0.75)
        for col in range(3):
            axes[row, col].set_xlabel(r"$\theta_1$")
            axes[row, col].set_ylabel(r"$\theta_3$")
    fig.savefig(FIGURES / "exchange_asymmetry_rebuilt.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "exchange_asymmetry_rebuilt.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_current(rows: list[dict]):
    selected = [r for r in rows if r["case"] == "driven"
                and r["source"] in ("direct_trajectory", "new_ensemble", "legacy_density")]
    sources = ("direct_trajectory", "new_ensemble", "legacy_density")
    labels = ("held-out", "validated density", "legacy density")
    colors = ("black", "#0072B2", "#D55E00")
    observables = ("J12", "J23", "J12_plus_J23")
    x = np.arange(len(observables))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    width = 0.23
    for index, (source, label, color) in enumerate(zip(sources, labels, colors)):
        source_rows = [next(r for r in selected if r["source"] == source
                            and r["observable"] == observable)
                       for observable in observables]
        value = np.asarray([float(r["estimate"]) for r in source_rows])
        lower = value - np.asarray([float(r["ci_low"]) for r in source_rows])
        upper = np.asarray([float(r["ci_high"]) for r in source_rows]) - value
        ax.bar(x + (index - 1) * width, value, width, color=color, alpha=0.82, label=label)
        ax.errorbar(x + (index - 1) * width, value, yerr=np.vstack((lower, upper)),
                    fmt="none", ecolor=color, capsize=2.5, lw=1)
    ax.axhline(0, color="0.35", lw=0.8)
    ax.set_xticks(x, (r"$J_{12}$", r"$J_{23}$", r"$J_{12}+J_{23}$"))
    ax.set_ylabel("stationary expectation")
    ax.legend(frameon=False, ncol=3)
    ax.grid(axis="y", alpha=0.18)
    fig.tight_layout()
    fig.savefig(FIGURES / "current_balance_rebuilt.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "current_balance_rebuilt.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def process_case(case: str, new_density: FinalNESSDensity, old_model,
                 old_mu: float, old_sigma: float):
    validation_rows, _ = evaluation_roles(case, "fine", "validation")
    validation_target, _ = load_rows(validation_rows)
    support = freeze_support(validation_target)
    del validation_target

    target_rows, reference_rows = evaluation_roles(case, "fine", "test")
    target, target_streams = load_rows(target_rows)
    reference, reference_streams = load_rows(reference_rows)
    target_mask = support.mask(target)

    proposal_logp = -energy(reference) / 5.0 - new_density.log_gibbs_partition
    if case == "driven":
        component_logp = new_density.component_log_density(reference)
    else:
        component_logp = np.repeat(proposal_logp[:, None], 3, axis=1)
    new_weights = [scaled_weights(component_logp[:, index] - proposal_logp)
                   for index in range(3)]
    ensemble_weights = np.mean(
        [weights / np.mean(weights) for weights in new_weights], axis=0
    )

    old_logp_reference = predict_old_log_density(old_model, reference, old_mu, old_sigma)
    old_weights = scaled_weights(old_logp_reference - proposal_logp)
    weight_arrays = {
        **{f"new_seed_{seed}": split_like(weights, reference_streams)
           for seed, weights in zip(NEW_SEEDS, new_weights)},
        "new_ensemble": split_like(ensemble_weights, reference_streams),
        "legacy_density": split_like(old_weights, reference_streams),
    }

    support_output = support_rows(case, support, target_streams, reference_streams, weight_arrays)

    # Pointwise endpoint exchange asymmetry on every masked held-out target state.
    x = target[target_mask]
    sx = exchange(x)
    point_rows = []
    if case == "driven":
        new_x = new_density.component_log_density(x)
        new_sx = new_density.component_log_density(sx)
    else:
        exact_x = -energy(x) / 5.0 - new_density.log_gibbs_partition
        exact_sx = -energy(sx) / 5.0 - new_density.log_gibbs_partition
        new_x = np.repeat(exact_x[:, None], 3, axis=1)
        new_sx = np.repeat(exact_sx[:, None], 3, axis=1)
    for index, seed in enumerate(NEW_SEEDS):
        point_rows.append(asymmetry_summary(
            case, f"new_seed_{seed}", relative_asymmetry(new_x[:, index], new_sx[:, index])
        ))
    ensemble_x = logsumexp(new_x, axis=1) - math.log(3.0)
    ensemble_sx = logsumexp(new_sx, axis=1) - math.log(3.0)
    point_rows.append(asymmetry_summary(
        case, "new_ensemble", relative_asymmetry(ensemble_x, ensemble_sx)
    ))
    old_x = predict_old_log_density(old_model, x, old_mu, old_sigma)
    old_sx = predict_old_log_density(old_model, sx, old_mu, old_sigma)
    point_rows.append(asymmetry_summary(
        case, "legacy_density", relative_asymmetry(old_x, old_sx)
    ))
    add_seed_ranges(point_rows, case, "new_ensemble", [f"new_seed_{s}" for s in NEW_SEEDS])

    # Integrated angular asymmetry.
    tv_rows = []
    histograms = {}
    direct_hist = angular_hist_streams(target_streams, support, None)
    row, hist, direct_tv_boot = angular_tv_row(case, "direct_trajectory", direct_hist)
    tv_rows.append(row); histograms[(case, "direct_trajectory")] = hist
    tv_boot = {"direct_trajectory": direct_tv_boot}
    for source, weights in weight_arrays.items():
        per_stream = angular_hist_streams(reference_streams, support, weights)
        row, hist, values = angular_tv_row(case, source, per_stream)
        tv_rows.append(row); histograms[(case, source)] = hist; tv_boot[source] = values
    ensemble_row = next(r for r in tv_rows if r["source"] == "new_ensemble")
    seed_values = [next(r for r in tv_rows if r["source"] == f"new_seed_{s}")["tv_exchange"]
                   for s in NEW_SEEDS]
    ensemble_row["seed_min"] = min(seed_values); ensemble_row["seed_max"] = max(seed_values)

    # Conditional phase locking and current balance.
    phase_output = []
    phase_hist_output = []
    phase_summaries = {"direct_trajectory": phase_stream_summaries(target_streams, support, None)}
    for source, weights in weight_arrays.items():
        phase_summaries[source] = phase_stream_summaries(reference_streams, support, weights)
    for source, summary in phase_summaries.items():
        phase_output += phase_rows(case, source, summary)
        if source in ("direct_trajectory", "new_ensemble", "legacy_density"):
            phase_hist_output += phase_hist_rows(case, source, summary)
    add_phase_seed_ranges(phase_output, case)

    current_output = []
    current_boot = {}
    sums, denom = current_stream_summaries(target_streams, support, None)
    rows, values = current_rows(case, "direct_trajectory", sums, denom)
    current_output += rows; current_boot["direct_trajectory"] = values
    for source, weights in weight_arrays.items():
        sums, denom = current_stream_summaries(reference_streams, support, weights)
        rows, values = current_rows(case, source, sums, denom)
        current_output += rows; current_boot[source] = values
    add_current_seed_ranges(current_output, case)

    return {
        "support": support_output,
        "pointwise": point_rows,
        "tv": tv_rows,
        "tv_boot": tv_boot,
        "histograms": histograms,
        "phase": phase_output,
        "phase_hist": phase_hist_output,
        "phase_summaries": phase_summaries,
        "current": current_output,
        "current_boot": current_boot,
    }


def main():
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    tf.keras.backend.set_floatx("float64")
    plt.rcParams.update({
        "font.family": "serif", "font.size": 8.5, "axes.labelsize": 9,
        "axes.titlesize": 9, "legend.fontsize": 7.5, "pdf.fonttype": 42,
        "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
    })

    new_density = FinalNESSDensity("fine")
    old_driven = tf.keras.models.load_model(
        OLD_DRIVEN_MODEL, custom_objects={"periodic_encode": periodic_encode}
    )
    old_equilibrium = tf.keras.models.load_model(
        OLD_EQUILIBRIUM_MODEL, custom_objects={"periodic_encode": periodic_encode}
    )
    old_driven_mu, old_driven_sigma = old_log_normalization(OLD_DRIVEN_DENSITY)
    old_eq_mu, old_eq_sigma = old_log_normalization(OLD_EQUILIBRIUM_DENSITY)

    driven = process_case("driven", new_density, old_driven, old_driven_mu, old_driven_sigma)
    equilibrium = process_case(
        "equilibrium", new_density, old_equilibrium, old_eq_mu, old_eq_sigma
    )
    results = (driven, equilibrium)

    support_rows_all = sum((result["support"] for result in results), [])
    point_rows = sum((result["pointwise"] for result in results), [])
    tv_rows = sum((result["tv"] for result in results), [])
    phase_rows_all = sum((result["phase"] for result in results), [])
    phase_hist_rows_all = sum((result["phase_hist"] for result in results), [])
    current_rows_all = sum((result["current"] for result in results), [])

    # Direct comparisons and fixed interpretation gates.
    comparison_rows = []
    for case, result in zip(("driven", "equilibrium"), results):
        for row in [r for r in phase_rows_all if r["case"] == case and r["source"] == "new_ensemble"]:
            direct = next(r for r in phase_rows_all if r["case"] == case
                          and r["source"] == "direct_trajectory"
                          and r["sector"] == row["sector"] and r["angle"] == row["angle"]
                          and r["metric"] == row["metric"])
            if row["metric"] in ("circular_mean", "histogram_mode"):
                difference = float(wrap(float(row["estimate"]) - float(direct["estimate"])))
                # Circular interval overlap is represented by whether zero lies in
                # the bootstrap-scale difference interval.  Conservative SE uses
                # the sum of the two marginal bootstrap variances.
                se = math.sqrt(
                    ((float(row["ci_high"]) - float(row["ci_low"])) / 3.92) ** 2
                    + ((float(direct["ci_high"]) - float(direct["ci_low"])) / 3.92) ** 2
                )
                z = difference / se if se > 0 else math.inf
                agrees = abs(z) <= 1.96
            else:
                difference = float(row["estimate"]) - float(direct["estimate"])
                se = math.sqrt(float(row.get("bootstrap_se", 0) or 0) ** 2)
                agrees = interval_overlap(row["ci_low"], row["ci_high"],
                                          direct["ci_low"], direct["ci_high"])
                width = math.sqrt(
                    ((float(row["ci_high"]) - float(row["ci_low"])) / 3.92) ** 2
                    + ((float(direct["ci_high"]) - float(direct["ci_low"])) / 3.92) ** 2
                )
                z = difference / width if width > 0 else math.inf
            comparison_rows.append({
                "case": case, "diagnostic": "phase_locking",
                "sector": row["sector"], "angle": row["angle"],
                "metric": row["metric"], "density_estimate": row["estimate"],
                "trajectory_estimate": direct["estimate"], "difference": difference,
                "difference_z_approx": z, "intervals_overlap": bool(agrees),
            })

        density_balance = next(r for r in current_rows_all if r["case"] == case
                               and r["source"] == "new_ensemble"
                               and r["observable"] == "J12_plus_J23")
        direct_balance = next(r for r in current_rows_all if r["case"] == case
                              and r["source"] == "direct_trajectory"
                              and r["observable"] == "J12_plus_J23")
        # Independent stream-bootstrap draws already use source-specific seeds.
        difference_boot = (result["current_boot"]["new_ensemble"][:, 2]
                           - result["current_boot"]["direct_trajectory"][:, 2])
        difference = float(density_balance["estimate"] - direct_balance["estimate"])
        z = difference / float(np.std(difference_boot, ddof=1))
        new_zero = float(density_balance["ci_low"]) <= 0 <= float(density_balance["ci_high"])
        direct_zero = float(direct_balance["ci_low"]) <= 0 <= float(direct_balance["ci_high"])
        comparison_rows.append({
            "case": case, "diagnostic": "current_balance", "sector": "all",
            "angle": "none", "metric": "E_J12_plus_J23",
            "density_estimate": density_balance["estimate"],
            "trajectory_estimate": direct_balance["estimate"],
            "difference": difference, "difference_z_approx": z,
            "intervals_overlap": bool(new_zero and direct_zero and abs(z) <= 1.96),
        })

    driven_tv = next(r for r in tv_rows if r["case"] == "driven" and r["source"] == "new_ensemble")
    eq_tv = next(r for r in tv_rows if r["case"] == "equilibrium" and r["source"] == "new_ensemble")
    asymmetry_gate = (
        min(float(next(r for r in tv_rows if r["case"] == "driven"
                       and r["source"] == f"new_seed_{s}")["tv_exchange"])
            for s in NEW_SEEDS) > float(eq_tv["ci_high"])
        and float(driven_tv["ci_low"]) > 3.0 * float(eq_tv["ci_high"])
    )
    current_gate = next(r for r in comparison_rows if r["case"] == "driven"
                        and r["diagnostic"] == "current_balance")["intervals_overlap"]
    phase_failures = [r for r in comparison_rows if r["case"] == "driven"
                      and r["diagnostic"] == "phase_locking"
                      and not r["intervals_overlap"]]
    verdict = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_only": True,
        "new_simulation": False,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "theta0_radians": THETA0,
        "theta0_degrees": math.degrees(THETA0),
        "asymmetry_gate_pass": bool(asymmetry_gate),
        "current_balance_gate_pass": bool(current_gate),
        "phase_comparison_count": len([r for r in comparison_rows if r["case"] == "driven"
                                       and r["diagnostic"] == "phase_locking"]),
        "phase_comparison_failures": len(phase_failures),
        "phase_density_trajectory_gate_pass": len(phase_failures) == 0,
        "overall_pass": bool(asymmetry_gate and current_gate and len(phase_failures) == 0),
        "claim_boundary": "validated on the explicit sampled-support mask only",
    }

    write_csv(ANALYSIS / "support_summary.csv", support_rows_all)
    write_csv(ANALYSIS / "asymmetry_pointwise.csv", point_rows)
    write_csv(ANALYSIS / "asymmetry_angular_tv.csv", tv_rows)
    write_csv(ANALYSIS / "phase_locking.csv", phase_rows_all)
    write_csv(ANALYSIS / "phase_histograms.csv", phase_hist_rows_all)
    write_csv(ANALYSIS / "current_balance.csv", current_rows_all)
    write_csv(ANALYSIS / "density_trajectory_comparisons.csv", comparison_rows)
    (ANALYSIS / "verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")

    histograms = {**driven["histograms"], **equilibrium["histograms"]}
    plot_asymmetry(histograms)
    plot_phase(phase_hist_rows_all)
    plot_current(current_rows_all)

    input_paths = [
        Path(__file__).resolve(), HERE / "PROTOCOL.md",
        DENSITY_ROOT / "V6_HOLDOUT_SPLITS.csv",
        DENSITY_ROOT / "final_analysis/density_manifest.json",
        OLD_DRIVEN_MODEL, OLD_EQUILIBRIUM_MODEL,
        OLD_DRIVEN_DENSITY, OLD_EQUILIBRIUM_DENSITY,
    ]
    for seed in NEW_SEEDS:
        input_paths += [
            DENSITY_ROOT / f"recovery_r11_fine/models/seed_{seed}/model.weights.h5",
            DENSITY_ROOT / f"recovery_r11_fine/models/seed_{seed}/model_config.json",
        ]
    hash_rows = [{
        "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path),
    } for path in input_paths]
    write_csv(ANALYSIS / "input_hashes.csv", hash_rows)
    provenance = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "tensorflow": tf.__version__,
        "numpy": np.__version__, "command": (
            "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 "
            "analyze_stabilization.py"
        ),
        "v6_base_seeds": {"driven": 2026091701, "equilibrium": 2026091702},
        "density_training_seeds": list(NEW_SEEDS),
        "bootstrap_seed": BOOTSTRAP_SEED,
    }
    (ANALYSIS / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
