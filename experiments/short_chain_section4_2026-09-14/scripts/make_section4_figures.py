#!/usr/bin/env python3
"""Create reproducible publication figures for the frozen Section-4 verdict."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
FIGURES = ROOT / "figures"
ANALYSIS = ROOT / "analysis"
COLORS = ("#0072B2", "#D55E00", "#009E73")


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def density_figure() -> None:
    choice = json.loads((ROOT / "stationary_density/FROZEN_DENSITY_CHOICE.json").read_text())
    rows = []
    labels = []
    for family_index, family in enumerate(choice["candidate_summaries"]):
        widths = "x".join(str(x) for x in family["widths"])
        label = rf"$K={family['components']}$" + "\n" + widths + "\n" + rf"$\lambda={family['lambda_fp']:g}$"
        labels.append(label)
        for seed_record in family["seeds"]:
            rows.append({
                "family_index": family_index,
                "components": family["components"],
                "widths": widths,
                "lambda_fp": family["lambda_fp"],
                "seed": seed_record["seed"],
                "gibbs_slope": seed_record["gibbs_bulk_log_slope"],
                "gibbs_rmse": seed_record["gibbs_bulk_centered_log_rmse"],
                "fp_median": seed_record["fp_abs_median"],
                "fp_p90": seed_record["fp_abs_p90"],
                "validation_nll": seed_record["best_validation_nll"],
                "admissible": int(seed_record["selection_admissible"]),
            })
    write_csv(ANALYSIS / "density_candidate_metrics.csv", rows)
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.35), sharex=True)
    offsets = (-0.17, 0.0, 0.17)
    for seed_index, seed in enumerate((4101, 4102, 4103)):
        selected = [r for r in rows if r["seed"] == seed]
        x = np.asarray([r["family_index"] for r in selected]) + offsets[seed_index]
        axes[0].scatter(x, [r["gibbs_slope"] for r in selected], s=28,
                        color=COLORS[seed_index], label=f"seed {seed}", zorder=3)
        axes[1].scatter(x, [r["fp_median"] for r in selected], s=28,
                        color=COLORS[seed_index], zorder=3)
        axes[2].scatter(x, [r["fp_p90"] for r in selected], s=28,
                        color=COLORS[seed_index], zorder=3)
    axes[0].axhspan(0.95, 1.05, color="#999999", alpha=0.16)
    axes[0].axhline(0.95, color="#555555", linestyle="--", linewidth=0.8)
    axes[0].axhline(1.05, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].axhline(0.20, color="black", linestyle="--", linewidth=1.0)
    axes[2].axhline(1.00, color="black", linestyle="--", linewidth=1.0)
    axes[0].set_ylabel(r"Gibbs log-density slope")
    axes[1].set_ylabel(r"median $|L^\dagger\rho/\rho|$")
    axes[2].set_ylabel(r"90th percentile $|L^\dagger\rho/\rho|$")
    axes[1].set_yscale("log"); axes[2].set_yscale("log")
    for axis in axes:
        axis.set_xticks(range(len(labels)), labels, fontsize=7)
        axis.grid(axis="y", alpha=0.22)
        axis.tick_params(direction="out")
    axes[0].legend(frameon=False, fontsize=7, loc="best")
    axes[1].annotate("gate = 0.20", (0.03, 0.07), xycoords="axes fraction", fontsize=7)
    axes[2].annotate("gate = 1.00", (0.03, 0.07), xycoords="axes fraction", fontsize=7)
    fig.tight_layout(w_pad=1.4)
    fig.savefig(FIGURES / "equilibrium_density_candidate_gates.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "equilibrium_density_candidate_gates.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def mode_figure() -> None:
    results = json.loads((ROOT / "modes/BLIND_MODE_RESULTS.json").read_text())["summaries"]
    selection = json.loads((ROOT / "modes/MODE_SELECTION.json").read_text())["selections"]
    rows = []
    order = (
        "driven_dt1e-3", "driven_dt2p5e-4",
        "equilibrium_dt1e-3", "equilibrium_dt2p5e-4",
    )
    labels = ("driven\n$10^{-3}$", "driven\n$2.5\\times10^{-4}$",
              "equilibrium\n$10^{-3}$", "equilibrium\n$2.5\\times10^{-4}$")
    for case in order:
        for kind in ("real", "complex"):
            record = results[f"{case}:{kind}"]
            selected = selection[f"{case}:{kind}"]
            rows.append({
                "case": case, "kind": kind,
                "blind_lambda_R": record["lambda_R"],
                "blind_lambda_R_ci_low": record["lambda_R_ci"][0],
                "blind_lambda_R_ci_high": record["lambda_R_ci"][1],
                "blind_lambda_I": record.get("lambda_I", 0.0),
                "blind_lambda_I_ci_low": record.get("lambda_I_ci", [0.0, 0.0])[0],
                "blind_lambda_I_ci_high": record.get("lambda_I_ci", [0.0, 0.0])[1],
                "selected_edmd_lambda_R": selected["lambda"][0],
                "selected_edmd_lambda_I": selected["lambda"][1],
                "independent_lambda_I": record.get("independent_lambda_I", ""),
                "test_residual": record["test_residual"],
                "validation_residual": record["validation_residual"],
                "reportable": int(record.get("reportable_edmd_visible_mode",
                                               record.get("reportable_oscillatory_mode", False))),
            })
    write_csv(ANALYSIS / "heldout_mode_summary.csv", rows)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    x = np.arange(4)
    real = [next(r for r in rows if r["case"] == case and r["kind"] == "real") for case in order]
    y = np.asarray([r["blind_lambda_R"] for r in real])
    lo = np.asarray([r["blind_lambda_R_ci_low"] for r in real])
    hi = np.asarray([r["blind_lambda_R_ci_high"] for r in real])
    axes[0].errorbar(x, y, yerr=[y-lo, hi-y], fmt="o", color=COLORS[0],
                     capsize=3, label="blind autocorrelation")
    axes[0].scatter(x, [r["selected_edmd_lambda_R"] for r in real], marker="x",
                    s=42, color=COLORS[1], label="training EDMD eigenvalue")
    axes[0].set_ylabel(r"slow real rate $\operatorname{Re}\lambda$")
    axes[0].legend(frameon=False, fontsize=7)

    complex_rows = [next(r for r in rows if r["case"] == case and r["kind"] == "complex") for case in order]
    yi = np.asarray([r["blind_lambda_I"] for r in complex_rows])
    loi = np.asarray([r["blind_lambda_I_ci_low"] for r in complex_rows])
    hii = np.asarray([r["blind_lambda_I_ci_high"] for r in complex_rows])
    axes[1].errorbar(x, yi, yerr=[yi-loi, hii-yi], fmt="o", color=COLORS[0],
                     capsize=3, label="blind damped fit")
    axes[1].scatter(x, [r["selected_edmd_lambda_I"] for r in complex_rows],
                    marker="x", s=42, color=COLORS[1], label="training EDMD")
    axes[1].scatter(x, [r["independent_lambda_I"] for r in complex_rows],
                    marker="s", facecolors="none", edgecolors=COLORS[2], s=38,
                    label="independent observable fit")
    axes[1].set_ylabel(r"oscillation frequency $|\operatorname{Im}\lambda|$")
    axes[1].set_ylim(0, 7.5)
    axes[1].legend(frameon=False, fontsize=7)
    for axis in axes:
        axis.set_xticks(x, labels)
        axis.grid(axis="y", alpha=0.22)
        axis.axvline(1.5, color="#AAAAAA", linewidth=0.7)
    fig.tight_layout(w_pad=1.8)
    fig.savefig(FIGURES / "heldout_relaxation_modes.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "heldout_relaxation_modes.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def modal_weight_figure() -> None:
    with (ROOT / "modes/blind_modal_weights.csv").open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["kind"] == "real"]
    observables = ("I2", "I1_plus_I3", "cos_theta3", "cos_theta1_plus_cos_theta3")
    labels = (r"$I_2$", r"$I_1+I_3$", r"$\cos\theta_3$",
              r"$\cos\theta_1+\cos\theta_3$")
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.35), sharey=True)
    for axis, bath in zip(axes, ("driven", "equilibrium")):
        x = np.arange(len(observables))
        width = 0.34
        for timestep_index, suffix in enumerate(("dt1e-3", "dt2p5e-4")):
            selected = [next(r for r in rows if r["case"] == f"{bath}_{suffix}"
                             and r["observable"] == observable)
                        for observable in observables]
            y = np.asarray([float(r["weight_abs"]) for r in selected])
            lo = np.asarray([float(r["weight_abs_ci_low"]) for r in selected])
            hi = np.asarray([float(r["weight_abs_ci_high"]) for r in selected])
            position = x + (timestep_index-0.5)*width
            axis.bar(position, y, width=width, color=COLORS[timestep_index], alpha=0.86,
                     label=(r"$dt=10^{-3}$" if timestep_index == 0
                            else r"$dt=2.5\times10^{-4}$"))
            axis.errorbar(position, y, yerr=[y-lo, hi-y], fmt="none", ecolor="black",
                          elinewidth=0.8, capsize=2)
        axis.set_xticks(x, labels, rotation=17, ha="right")
        axis.set_title(bath.capitalize())
        axis.grid(axis="y", alpha=0.22)
        axis.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("absolute modal weight")
    fig.tight_layout(w_pad=1.2)
    fig.savefig(FIGURES / "blind_slow_mode_weights.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "blind_slow_mode_weights.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    ANALYSIS.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "serif", "font.size": 9, "axes.labelsize": 9,
        "axes.titlesize": 9, "legend.fontsize": 8,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    density_figure(); mode_figure(); modal_weight_figure()


if __name__ == "__main__":
    main()
