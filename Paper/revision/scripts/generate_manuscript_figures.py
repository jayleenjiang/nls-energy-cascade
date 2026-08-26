#!/usr/bin/env python3
"""Generate manuscript figures from local NLS numerical-study artifacts.

The script is intentionally read-mostly: it consumes profile text files,
CKSTT embedding CSVs, and dense LTE histograms already present in the
repository, then writes figure files into ``Paper/revision_2026-06-19``.
It does not modify the source data.
"""

from __future__ import annotations

import csv
import json
from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]


def read_profile(path: Path) -> tuple[int | None, np.ndarray, np.ndarray]:
    samples: int | None = None
    rows: list[tuple[int, float]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "total_samples=" in line:
                token = line.split("total_samples=", 1)[1].split()[0]
                try:
                    samples = int(float(token))
                except ValueError:
                    samples = None
            continue
        parts = line.split()
        if len(parts) >= 2:
            rows.append((int(parts[0]), float(parts[1])))
    if not rows:
        raise ValueError(f"No profile rows found in {path}")
    idx = np.array([r[0] for r in rows], dtype=float)
    vals = np.array([r[1] for r in rows], dtype=float)
    return samples, idx, vals


def make_action_profile_figure() -> dict[str, object]:
    profiles = [
        ("n=25", ROOT / "experiments/lte/test_profile.txt"),
        ("n=50", ROOT / "experiments/lte/simd_n50_profile.txt"),
        ("n=100", ROOT / "experiments/lte/n100_dt25_profile.txt"),
    ]

    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    metrics: dict[str, object] = {}
    for label, path in profiles:
        samples, idx, vals = read_profile(path)
        x = idx / idx[-1]
        ax.plot(x, vals, marker="o", ms=3.4, lw=1.7, label=label)
        metrics[label] = {
            "source": str(path.relative_to(ROOT)),
            "samples": samples,
            "first": float(vals[0]),
            "middle": float(vals[len(vals) // 2]),
            "last": float(vals[-1]),
            "min": float(vals.min()),
            "max": float(vals.max()),
        }

    ax.set_xlabel(r"normalized mode index $j/(n-1)$")
    ax.set_ylabel(r"mean action $\langle I_j\rangle$")
    ax.set_title(r"Closed-chain nonequilibrium action profiles ($T_1=10$, $T_n=2$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()

    out_pdf = REVISION / "action_profiles.pdf"
    out_png = REVISION / "action_profiles.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {"files": [str(out_pdf.relative_to(ROOT)), str(out_png.relative_to(ROOT))], "metrics": metrics}


def frac_to_float(value: str) -> float:
    return float(Fraction(value))


def read_embedding_transition(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "family_id": int(row["family_id"]),
                    "point_type": row["point_type"],
                    "real": frac_to_float(row["real"]),
                    "imag": frac_to_float(row["imag"]),
                }
            )
    return rows


def make_embedding_figure() -> dict[str, object]:
    base = ROOT / "Energy Cascade/Rect Plots"
    files = [
        ("$\\Lambda_1\\to\\Lambda_2$", base / "rectangles_N5_R100_f1_to_f2.csv"),
        ("$\\Lambda_2\\to\\Lambda_3$", base / "rectangles_N5_R100_f2_to_f3.csv"),
        ("$\\Lambda_3\\to\\Lambda_4$", base / "rectangles_N5_R100_f3_to_f4.csv"),
        ("$\\Lambda_4\\to\\Lambda_5$", base / "rectangles_N5_R100_f4_to_f5.csv"),
    ]

    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    colors = plt.cm.viridis(np.linspace(0.08, 0.88, len(files)))
    all_points: list[tuple[float, float]] = []

    for (label, path), color in zip(files, colors):
        rows = read_embedding_transition(path)
        by_family: dict[int, list[dict[str, object]]] = {}
        for row in rows:
            by_family.setdefault(int(row["family_id"]), []).append(row)
            all_points.append((float(row["real"]), float(row["imag"])))
        for family_rows in by_family.values():
            parents = [r for r in family_rows if r["point_type"] == "parent"]
            children = [r for r in family_rows if r["point_type"] == "child"]
            pts = parents + children
            xs = [float(p["real"]) for p in pts]
            ys = [float(p["imag"]) for p in pts]
            if len(xs) == 4:
                order = [0, 2, 1, 3, 0]
                ax.plot([xs[i] for i in order], [ys[i] for i in order], color=color, lw=1.0, alpha=0.55)
            ax.scatter(xs, ys, s=16, color=color, edgecolor="black", linewidth=0.25, label=label)

    # Deduplicate legend labels.
    handles, labels = ax.get_legend_handles_labels()
    seen: set[str] = set()
    uniq = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
    ax.legend([h for h, _ in uniq], [l for _, l in uniq], frameon=False, fontsize=9, loc="upper right")
    ax.axhline(0, color="0.65", lw=0.8)
    ax.axvline(0, color="0.65", lw=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$\operatorname{Re} K$")
    ax.set_ylabel(r"$\operatorname{Im} K$")
    ax.set_title("CKSTT resonant-family embedding, first five generations")
    ax.grid(alpha=0.22)
    fig.tight_layout()

    out_pdf = REVISION / "cascade_embedding.pdf"
    out_png = REVISION / "cascade_embedding.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    xs = [p[0] for p in all_points]
    ys = [p[1] for p in all_points]
    return {
        "files": [str(out_pdf.relative_to(ROOT)), str(out_png.relative_to(ROOT))],
        "metrics": {
            "source_files": [str(path.relative_to(ROOT)) for _, path in files],
            "point_count_with_repetitions": len(all_points),
            "real_range": [float(min(xs)), float(max(xs))],
            "imag_range": [float(min(ys)), float(max(ys))],
        },
    }


def load_hist(path: Path) -> dict[str, object]:
    meta: dict[str, object] = {}
    lines = path.read_text().splitlines()
    data_start = 0
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if parts[0] == "NB":
            meta["NB"] = int(parts[1])
        elif parts[0] in {"I_LO", "I_HI", "TH_LO", "TH_HI"}:
            meta[parts[0]] = float(parts[1])
        elif parts[0] in {"TOTAL", "OVERFLOW"}:
            meta[parts[0]] = int(parts[1])
        elif len(parts) == 4:
            data_start = i
            break
    nb = int(meta["NB"])
    counts = np.zeros((nb, nb, nb), dtype=np.int64)
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 4:
            continue
        ia, ib, it, count = map(int, parts)
        counts[ia, ib, it] = count
    meta["counts"] = counts
    return meta


def hist_coordinates(hist: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    nb = int(hist["NB"])
    i_lo = float(hist["I_LO"])
    i_hi = float(hist["I_HI"])
    th_lo = float(hist["TH_LO"])
    th_hi = float(hist["TH_HI"])
    d_i = (i_hi - i_lo) / nb
    d_th = (th_hi - th_lo) / nb
    i_centers = i_lo + d_i * (np.arange(nb) + 0.5)
    th_centers = th_lo + d_th * (np.arange(nb) + 0.5)
    return i_centers, i_centers.copy(), th_centers, d_i * d_i * d_th


def fit_rescaled_equilibrium(q_hist: dict[str, object], p_hist: dict[str, object], min_count: int = 50, i_max: float = 2.5) -> dict[str, object]:
    q_counts = q_hist["counts"]
    p_counts = p_hist["counts"]
    assert isinstance(q_counts, np.ndarray)
    assert isinstance(p_counts, np.ndarray)
    ia, ib, th, bin_vol = hist_coordinates(q_hist)
    iam, ibm, _ = np.meshgrid(ia, ib, th, indexing="ij")

    q_total = float(q_counts.sum())
    p_total = float(p_counts.sum())
    q_density = q_counts / q_total / bin_vol
    p_density = p_counts / p_total / bin_vol

    mask = (q_counts >= min_count) & (p_counts >= min_count) & (iam < i_max) & (ibm < i_max)
    y = np.log(q_density[mask])
    x = np.log(p_density[mask])
    w = q_counts[mask].astype(float)

    design = np.c_[np.ones_like(x), x]
    xtw = design.T * w
    beta = np.linalg.solve(xtw @ design, xtw @ y)
    fitted = design @ beta
    ybar = np.sum(w * y) / np.sum(w)
    r2 = 1.0 - np.sum(w * (y - fitted) ** 2) / np.sum(w * (y - ybar) ** 2)
    residual = np.full_like(q_density, np.nan, dtype=float)
    valid = (q_counts > 0) & (p_counts > 0)
    residual[valid] = np.log(q_density[valid]) - (beta[0] + beta[1] * np.log(p_density[valid]))

    return {
        "intercept": float(beta[0]),
        "slope": float(beta[1]),
        "weighted_r2": float(r2),
        "bins_used": int(mask.sum()),
        "residual": residual,
        "q_counts": q_counts,
        "p_counts": p_counts,
        "i_centers": ia,
        "theta_centers": th,
    }


def make_lte_residual_figure() -> dict[str, object]:
    cases = [
        (
            "n=50, $a/n=0.48$",
            ROOT / "lte/n50 data/simd_n50_j24.hist",
            ROOT / "lte/n50 data/simd_n50_eq_j24.hist",
        ),
        (
            "n=100, $a/n=0.48$",
            ROOT / "lte/n100 data/n100_j48.hist",
            ROOT / "lte/n100 data/n100_eq_j48.hist",
        ),
    ]

    results = []
    for label, q_path, p_path in cases:
        q_hist = load_hist(q_path)
        p_hist = load_hist(p_path)
        fit = fit_rescaled_equilibrium(q_hist, p_hist)
        fit["label"] = label
        fit["q_source"] = str(q_path.relative_to(ROOT))
        fit["p_source"] = str(p_path.relative_to(ROOT))
        results.append(fit)

    theta_index = int(np.argmin(np.abs(results[0]["theta_centers"])))
    theta_value = float(results[0]["theta_centers"][theta_index])

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), sharex=True, sharey=True)
    vmax = 0.45
    image = None
    metrics: dict[str, object] = {"theta_slice": theta_value}

    for ax, fit in zip(axes, results):
        residual = fit["residual"][:, :, theta_index]
        q_counts = fit["q_counts"][:, :, theta_index]
        p_counts = fit["p_counts"][:, :, theta_index]
        i_centers = fit["i_centers"]
        masked = np.ma.array(residual, mask=(q_counts < 200) | (p_counts < 200) | ~np.isfinite(residual))
        image = ax.imshow(
            masked.T,
            origin="lower",
            extent=[i_centers[0], i_centers[-1], i_centers[0], i_centers[-1]],
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_title(f"{fit['label']}\n$x={fit['slope']:.4f}$, $R_w^2={fit['weighted_r2']:.4f}$")
        ax.set_xlabel(r"$I_a$")
        ax.grid(alpha=0.12)
        metrics[str(fit["label"])] = {
            "q_source": fit["q_source"],
            "p_source": fit["p_source"],
            "slope": fit["slope"],
            "intercept": fit["intercept"],
            "weighted_r2": fit["weighted_r2"],
            "bins_used": fit["bins_used"],
            "displayed_bins": int(np.ma.count(masked)),
            "masked_threshold": "q_count>=200 and p_count>=200 at theta slice",
        }
    axes[0].set_ylabel(r"$I_b$")
    assert image is not None
    cbar = fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.82)
    cbar.set_label(r"$\log q - (x\log p_6+c)$")
    fig.suptitle(r"Mid-chain LTE residual at fixed $\theta\approx 0$", y=0.98)
    fig.subplots_adjust(top=0.78, right=0.88, wspace=0.12)

    out_pdf = REVISION / "lte_residual_midchain.pdf"
    out_png = REVISION / "lte_residual_midchain.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {"files": [str(out_pdf.relative_to(ROOT)), str(out_png.relative_to(ROOT))], "metrics": metrics}


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 10.5,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "figure.dpi": 140,
            "savefig.dpi": 220,
        }
    )
    outputs = {
        "action_profiles": make_action_profile_figure(),
        "cascade_embedding": make_embedding_figure(),
        "lte_residual_midchain": make_lte_residual_figure(),
    }
    metrics_path = REVISION / "manuscript_figure_metrics.json"
    metrics_path.write_text(json.dumps(outputs, indent=2) + "\n")
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
