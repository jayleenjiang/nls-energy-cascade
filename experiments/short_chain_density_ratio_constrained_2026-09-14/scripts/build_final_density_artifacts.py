#!/usr/bin/env python3
"""Build final normalized-density tables and publication figures."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from scipy.special import logsumexp

from triple_calibrated_transport_ratio_model import load_bundle
from v4_common import evaluation_roles as v4_roles, load_rows as v4_load
from v6_common import evaluation_roles as v6_roles, load_rows as v6_load

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "final_analysis"
FIGURES = ROOT / "final_figures"
COLORS = ("#0072B2", "#D55E00", "#009E73")
FINE_MODELS = [ROOT / f"recovery_r11_fine/models/seed_{s}" for s in (8201, 8202, 8203)]
COARSE_MODELS = [ROOT / f"recovery_r10_coarse/models/seed_{s}" for s in (8201, 8202, 8203)]


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def log_ratio(model, state, batch=8192):
    return np.concatenate([
        model.log_ratio(tf.constant(state[start:start + batch], tf.float64)).numpy()
        for start in range(0, len(state), batch)
    ])


def energy(state):
    i1, i2, i3, th1, th3 = state.T
    mass = i1 + i2 + i3
    return (0.5 * mass * mass - 0.25 * (i1 * i1 + i2 * i2 + i3 * i3)
            + i1 * i2 * np.cos(th1) + i2 * i3 * np.cos(th3))


def load_test(timestep):
    if timestep == "fine":
        tr, rr = v6_roles("driven", "fine", "test")
        target, target_streams = v6_load(tr); reference, reference_streams = v6_load(rr)
        paths = FINE_MODELS
    else:
        tr, rr = v4_roles("driven", "coarse", "test")
        target, target_streams = v4_load(tr); reference, reference_streams = v4_load(rr)
        paths = COARSE_MODELS
    return target, target_streams, reference, reference_streams, [load_bundle(p) for p in paths]


def normalizers(models, reference_streams, timestep):
    rows = []; summaries = []
    for index, model in enumerate(models):
        stream_means = []
        for stream in reference_streams:
            stream_means.append(float(np.exp(logsumexp(log_ratio(model, stream)) - math.log(len(stream)))))
        stream_means = np.asarray(stream_means)
        point = float(stream_means.mean()); log_point = math.log(point)
        seed = int(hashlib.sha256(f"normalizer|{timestep}|{index}".encode()).hexdigest()[:16], 16)
        rng = np.random.default_rng(seed)
        draws = np.empty(1000)
        for b in range(1000):
            draws[b] = math.log(stream_means[rng.integers(0, len(stream_means), len(stream_means))].mean())
        lo, hi = np.quantile(draws, (0.025, 0.975))
        summaries.append({
            "timestep": timestep, "model": index + 1,
            "ratio_normalizer": point, "log_ratio_normalizer": log_point,
            "log_ratio_ci_low": lo, "log_ratio_ci_high": hi,
            "bootstrap_replicates": 1000,
        })
        for stream, value in enumerate(stream_means):
            rows.append({"timestep": timestep, "model": index + 1,
                         "stream_index": stream, "mean_exp_log_ratio": value})
    return summaries, rows


def marginal_data(target, reference, model_weights, timestep):
    rows = []
    specs = [
        ("log_I1", np.log(target[:, 0]), np.log(reference[:, 0]), None),
        ("log_I2", np.log(target[:, 1]), np.log(reference[:, 1]), None),
        ("log_I3", np.log(target[:, 2]), np.log(reference[:, 2]), None),
        ("theta1", target[:, 3], reference[:, 3], np.linspace(-np.pi, np.pi, 73)),
        ("theta3", target[:, 4], reference[:, 4], np.linspace(-np.pi, np.pi, 73)),
    ]
    for name, tv, rv, fixed_edges in specs:
        edges = fixed_edges
        if edges is None:
            low, high = np.quantile(tv, (0.001, 0.999))
            edges = np.linspace(low, high, 73)
        ht, _ = np.histogram(tv, edges, density=True)
        seed_hist = []
        for weights in model_weights:
            h, _ = np.histogram(rv, edges, weights=weights, density=True)
            seed_hist.append(h)
        seed_hist = np.asarray(seed_hist)
        centers = 0.5 * (edges[:-1] + edges[1:])
        for j, center in enumerate(centers):
            rows.append({
                "timestep": timestep, "marginal": name, "bin": j,
                "center": center, "left": edges[j], "right": edges[j + 1],
                "target_density": ht[j], "reconstructed_density": seed_hist[:, j].mean(),
                "seed_min_density": seed_hist[:, j].min(),
                "seed_max_density": seed_hist[:, j].max(),
            })
    return rows


def plot_marginals(rows, timestep):
    labels = {"log_I1": r"$\log I_1$", "log_I2": r"$\log I_2$",
              "log_I3": r"$\log I_3$", "theta1": r"$\theta_1$",
              "theta3": r"$\theta_3$"}
    fig, axes = plt.subplots(2, 3, figsize=(9.0, 5.2))
    for axis, name in zip(axes.flat, labels):
        selected = [r for r in rows if r["timestep"] == timestep and r["marginal"] == name]
        x = np.asarray([r["center"] for r in selected])
        target = np.asarray([r["target_density"] for r in selected])
        fit = np.asarray([r["reconstructed_density"] for r in selected])
        lo = np.asarray([r["seed_min_density"] for r in selected])
        hi = np.asarray([r["seed_max_density"] for r in selected])
        axis.plot(x, target, color="black", lw=1.4, label="direct NESS")
        axis.plot(x, fit, color=COLORS[0], lw=1.4, label="density reconstruction")
        axis.fill_between(x, lo, hi, color=COLORS[0], alpha=0.18, lw=0)
        axis.set_xlabel(labels[name]); axis.set_ylabel("density"); axis.grid(alpha=0.18)
    axes.flat[-1].axis("off")
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Independent {timestep}-timestep test", y=0.995, fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / f"ness_marginals_{timestep}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"ness_marginals_{timestep}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def angle_joint(target, reference, weights, timestep):
    edges = np.linspace(-np.pi, np.pi, 61)
    ht, _, _ = np.histogram2d(target[:, 3], target[:, 4], bins=(edges, edges), density=True)
    recon = []
    for w in weights:
        h, _, _ = np.histogram2d(reference[:, 3], reference[:, 4], bins=(edges, edges), weights=w, density=True)
        recon.append(h)
    hr = np.mean(recon, axis=0); diff = hr - ht
    rows = []
    for i in range(len(edges) - 1):
        for j in range(len(edges) - 1):
            rows.append({
                "timestep": timestep, "i": i, "j": j,
                "theta1": 0.5 * (edges[i] + edges[i + 1]),
                "theta3": 0.5 * (edges[j] + edges[j + 1]),
                "target_density": ht[i, j], "reconstructed_density": hr[i, j],
                "difference": diff[i, j],
            })
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 2.9), constrained_layout=True)
    vmax = max(ht.max(), hr.max()); extent=(-np.pi, np.pi, -np.pi, np.pi)
    for axis, data, title in zip(axes[:2], (ht.T, hr.T), ("direct NESS", "reconstruction")):
        image = axis.imshow(data, origin="lower", extent=extent, cmap="viridis", vmin=0, vmax=vmax, aspect="equal")
        axis.set_title(title); axis.set_xlabel(r"$\theta_1$"); axis.set_ylabel(r"$\theta_3$")
    limit = np.max(np.abs(diff))
    image_diff = axes[2].imshow(diff.T, origin="lower", extent=extent, cmap="coolwarm",
                                vmin=-limit, vmax=limit, aspect="equal")
    axes[2].set_title("reconstruction $-$ direct"); axes[2].set_xlabel(r"$\theta_1$")
    fig.colorbar(image, ax=axes[:2], shrink=0.78, label="density")
    fig.colorbar(image_diff, ax=axes[2], shrink=0.78, label="difference")
    fig.savefig(FIGURES / f"ness_angle_joint_{timestep}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"ness_angle_joint_{timestep}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return rows


def density_slices(models, logz, logz_eq):
    grid = np.linspace(-np.pi, np.pi, 121)
    th1, th3 = np.meshgrid(grid, grid, indexing="ij")
    rows = []; panels = []
    for action in (1.0, 2.0, 4.0):
        state = np.column_stack((
            np.full(th1.size, action), np.full(th1.size, action),
            np.full(th1.size, action), th1.ravel(), th3.ravel(),
        ))
        components = []
        for model, normalization in zip(models, logz):
            components.append(-energy(state) / 5.0 - logz_eq
                              + log_ratio(model, state) - normalization)
        log_density = logsumexp(np.asarray(components), axis=0) - math.log(len(models))
        panel = log_density.reshape(th1.shape); panels.append(panel)
        for index, value in enumerate(log_density):
            rows.append({
                "I1_I2_I3": action, "theta1": state[index, 3],
                "theta3": state[index, 4], "log_density": value,
                "density": math.exp(value),
            })
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.9), constrained_layout=True)
    vmin = min(p.min() for p in panels); vmax = max(p.max() for p in panels)
    for axis, action, panel in zip(axes, (1, 2, 4), panels):
        image = axis.imshow(panel.T, origin="lower", extent=(-np.pi, np.pi, -np.pi, np.pi),
                            cmap="magma", vmin=vmin, vmax=vmax, aspect="equal")
        axis.set_title(rf"$I_1=I_2=I_3={action}$")
        axis.set_xlabel(r"$\theta_1$"); axis.set_ylabel(r"$\theta_3$")
    fig.colorbar(image, ax=axes, shrink=0.78, label=r"absolute $\log\rho_{\rm ss}$")
    fig.savefig(FIGURES / "ness_density_angle_slices.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "ness_density_angle_slices.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return rows


def main():
    tf.keras.backend.set_floatx("float64")
    ANALYSIS.mkdir(exist_ok=True); FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "serif", "font.size": 9, "axes.labelsize": 9,
        "axes.titlesize": 9, "legend.fontsize": 8, "pdf.fonttype": 42,
        "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
    })
    all_normalizers = []; all_stream_normalizers = []; all_marginals = []; all_joint = []
    fine_models = None; fine_logz = None
    for timestep in ("fine", "coarse"):
        target, _, reference, reference_streams, models = load_test(timestep)
        summaries, stream_rows = normalizers(models, reference_streams, timestep)
        all_normalizers += summaries; all_stream_normalizers += stream_rows
        logz = [row["log_ratio_normalizer"] for row in summaries]
        model_weights = [np.exp(log_ratio(model, reference) - z) for model, z in zip(models, logz)]
        rows = marginal_data(target, reference, model_weights, timestep)
        all_marginals += rows; plot_marginals(rows, timestep)
        all_joint += angle_joint(target, reference, model_weights, timestep)
        if timestep == "fine":
            fine_models = models; fine_logz = logz
        del target, reference, model_weights
    logz_eq = json.loads((ROOT / "absolute_normalization/result.json").read_text())["log_partition"]
    for row in all_normalizers:
        row["gibbs_log_partition"] = logz_eq
        row["absolute_log_partition"] = logz_eq + row["log_ratio_normalizer"]
        row["absolute_partition"] = math.exp(row["absolute_log_partition"])
    write_csv(ANALYSIS / "normalizer_summary.csv", all_normalizers)
    write_csv(ANALYSIS / "normalizer_stream_values.csv", all_stream_normalizers)
    write_csv(ANALYSIS / "marginal_histograms.csv", all_marginals)
    write_csv(ANALYSIS / "angle_joint_histograms.csv", all_joint)
    slices = density_slices(fine_models, fine_logz, logz_eq)
    write_csv(ANALYSIS / "density_angle_slices.csv", slices)
    result = {
        "density_measure": "dI1 dI2 dI3 dtheta1 dtheta3",
        "gibbs_temperature": 5.0, "gibbs_log_partition": logz_eq,
        "fine_models": [str(path) for path in FINE_MODELS],
        "coarse_models": [str(path) for path in COARSE_MODELS],
        "ensemble_definition": "arithmetic mean of three separately normalized densities",
    }
    (ANALYSIS / "density_manifest.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
