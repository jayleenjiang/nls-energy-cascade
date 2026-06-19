#!/usr/bin/env python3
"""Reproduce and stress-test the backward-Monte-Carlo relaxation-rate fit.

The historical eigenfunction workflow fixes lambda_R near -0.934.  That value
comes from pure-exponential fits

    y(t) = A exp(lambda_R t) + c

to the MC conditional expectations stored in KDE/NLS_backward_Y_train.txt, with
only fits satisfying -5 < lambda_R < 0 retained.  This script reproduces that
number and measures its sensitivity to the time window used in the fit.

Outputs are written inside Paper/revision_2026-06-19 so that the manuscript can
cite a source-traced diagnostic without modifying any original data files.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parents[3]
REVISION_DIR = ROOT / "Paper" / "revision_2026-06-19"
SOURCE = ROOT / "KDE" / "NLS_backward_Y_train.txt"
DT_EFF = 0.01

# Labels are manuscript-facing.  Indices use Python half-open convention.
WINDOWS = [
    ("0--5", 0, 501),
    ("0--2", 0, 201),
    ("0.1--2", 10, 201),
    ("0.25--3", 25, 301),
    ("0.5--5", 50, 501),
]


def pure_exp(t: np.ndarray, amplitude: float, lam_r: float, offset: float) -> np.ndarray:
    return amplitude * np.exp(lam_r * t) + offset


def fit_window(y_curves: np.ndarray, start: int, stop: int) -> dict[str, object]:
    """Fit one window and return retained rates plus summary statistics."""

    t = np.arange(start, stop, dtype=float) * DT_EFF
    rates: list[float] = []
    ssr: list[float] = []
    failures = 0

    for y_full in y_curves:
        y = y_full[start:stop]
        tail = y[-min(50, len(y)) :].mean()
        try:
            params, _ = curve_fit(
                pure_exp,
                t,
                y,
                p0=[0.5, -1.0, tail],
                maxfev=5000,
            )
        except Exception:
            failures += 1
            continue

        lam_r = float(params[1])
        if -5.0 < lam_r < 0.0 and np.isfinite(lam_r):
            fit = pure_exp(t, *params)
            rates.append(lam_r)
            ssr.append(float(np.mean((y - fit) ** 2)))

    rates_np = np.asarray(rates, dtype=float)
    ssr_np = np.asarray(ssr, dtype=float)
    if len(rates_np) == 0:
        raise RuntimeError(f"No valid fits for window {start}:{stop}")

    q25, q75 = np.percentile(rates_np, [25, 75])
    return {
        "start_index": start,
        "stop_index": stop,
        "time_start": start * DT_EFF,
        "time_stop": (stop - 1) * DT_EFF,
        "valid_fits": int(len(rates_np)),
        "total_curves": int(len(y_curves)),
        "failed_optimizations": int(failures),
        "median_lambda_R": float(np.median(rates_np)),
        "mean_lambda_R": float(np.mean(rates_np)),
        "q25_lambda_R": float(q25),
        "q75_lambda_R": float(q75),
        "median_ssr": float(np.median(ssr_np)),
        "rates": rates_np,
    }


def main() -> None:
    y_curves = np.loadtxt(SOURCE)
    summaries = []
    fit_outputs = []

    for label, start, stop in WINDOWS:
        result = fit_window(y_curves, start, stop)
        result["label"] = label
        fit_outputs.append(result)
        serializable = {k: v for k, v in result.items() if k != "rates"}
        summaries.append(serializable)

    output_json = REVISION_DIR / "eigen_fit_sensitivity.json"
    output_json.write_text(
        json.dumps(
            {
                "source": str(SOURCE.relative_to(ROOT)),
                "dt_eff": DT_EFF,
                "model": "A * exp(lambda_R * t) + c",
                "retention_rule": "-5 < lambda_R < 0",
                "summaries": summaries,
            },
            indent=2,
        )
        + "\n"
    )

    full = fit_outputs[0]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))

    axes[0].hist(full["rates"], bins=70, color="#4C78A8", edgecolor="white")
    axes[0].axvline(
        full["median_lambda_R"],
        color="#D62728",
        linestyle="--",
        linewidth=2,
        label=rf"median $={full['median_lambda_R']:.3f}$",
    )
    axes[0].set_xlabel(r"fitted $\lambda_R$")
    axes[0].set_ylabel("count")
    axes[0].set_title(r"Pure-exponential fits, $t\in[0,5]$")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    labels = [str(r["label"]) for r in fit_outputs]
    medians = np.array([float(r["median_lambda_R"]) for r in fit_outputs])
    q25 = np.array([float(r["q25_lambda_R"]) for r in fit_outputs])
    q75 = np.array([float(r["q75_lambda_R"]) for r in fit_outputs])
    x = np.arange(len(labels))
    axes[1].errorbar(
        x,
        medians,
        yerr=np.vstack([medians - q25, q75 - medians]),
        fmt="o",
        color="#4C78A8",
        ecolor="#9ECae9",
        elinewidth=4,
        capsize=4,
    )
    axes[1].axhline(full["median_lambda_R"], color="#D62728", linestyle="--", linewidth=1.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=25, ha="right")
    axes[1].set_xlabel("fit window")
    axes[1].set_ylabel(r"median $\lambda_R$ with IQR")
    axes[1].set_title("Window sensitivity")
    axes[1].grid(alpha=0.2)

    fig.suptitle(r"Backward Monte Carlo relaxation-rate diagnostic ($\gamma=0.1$)")
    fig.tight_layout()
    fig.savefig(REVISION_DIR / "eigenvalue_scatter.png", dpi=220)
    fig.savefig(REVISION_DIR / "eigenvalue_scatter.pdf")

    print(f"Wrote {output_json}")
    for row in summaries:
        print(
            f"{row['label']:>8s}: valid={row['valid_fits']:4d}/"
            f"{row['total_curves']} median={row['median_lambda_R']:.6f} "
            f"IQR=[{row['q25_lambda_R']:.6f}, {row['q75_lambda_R']:.6f}]"
        )


if __name__ == "__main__":
    main()
