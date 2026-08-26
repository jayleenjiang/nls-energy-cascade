#!/usr/bin/env python3
"""Supplementary two-tail and time-dependence plots for entropy FT runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from analyze_entropy_ft import aggregate_blocks, load_run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blocks", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--analysis-dir", type=Path)
    parser.add_argument(
        "--taus", default="20,40,60,80,100,120,140,160,180,200"
    )
    parser.add_argument("--threshold-step", type=float, default=0.01)
    parser.add_argument("--minimum-raw-count", type=int, default=5)
    return parser.parse_args()


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def survival_rows(
    values: np.ndarray,
    n: int,
    tau: float,
    step: float,
) -> list[dict[str, object]]:
    flattened = values.ravel()
    sample_count = flattened.size
    maximum = float(np.max(np.abs(flattened)))
    thresholds = np.arange(step, maximum + 0.5 * step, step)
    mu = float(np.mean(flattened))
    sigma = float(np.std(flattened, ddof=0))
    rows = []
    for threshold in thresholds:
        plus_count = int(np.count_nonzero(flattened >= threshold))
        minus_count = int(np.count_nonzero(flattened <= -threshold))
        plus_raw = plus_count / sample_count
        minus_raw = minus_count / sample_count
        rows.append(
            {
                "n": n,
                "tau": tau,
                "A": float(threshold),
                "n_samples": sample_count,
                "mean": mu,
                "std": sigma,
                "plus_count": plus_count,
                "minus_count": minus_count,
                "p_plus_raw": plus_raw,
                "p_minus_raw": minus_raw,
                "p_plus_plus_four": (plus_count + 2.0) / (sample_count + 4.0),
                "p_minus_plus_four": (minus_count + 2.0) / (sample_count + 4.0),
                "p_plus_normal": float(norm.sf(threshold, loc=mu, scale=sigma)),
                "p_minus_normal": float(norm.cdf(-threshold, loc=mu, scale=sigma)),
            }
        )
    return rows


def tail_fit_metrics(
    rows: list[dict[str, object]], minimum_raw_count: int
) -> dict[str, float]:
    usable_plus = [
        row
        for row in rows
        if int(row["plus_count"]) >= minimum_raw_count
        and float(row["p_plus_normal"]) > 0.0
    ]
    usable_minus = [
        row
        for row in rows
        if int(row["minus_count"]) >= minimum_raw_count
        and float(row["p_minus_normal"]) > 0.0
    ]

    def rmse(items: list[dict[str, object]], empirical: str, reference: str) -> float:
        if not items:
            return float("nan")
        differences = [
            math.log(float(row[empirical])) - math.log(float(row[reference]))
            for row in items
            if float(row[empirical]) > 0.0
        ]
        return float(math.sqrt(np.mean(np.square(differences)))) if differences else float("nan")

    return {
        "plus_log_survival_rmse_vs_normal": rmse(
            usable_plus, "p_plus_raw", "p_plus_normal"
        ),
        "minus_log_survival_rmse_vs_normal": rmse(
            usable_minus, "p_minus_raw", "p_minus_normal"
        ),
        "plus_thresholds_used": len(usable_plus),
        "minus_thresholds_used": len(usable_minus),
    }


def plot_two_tail_survival(
    output_dir: Path,
    n: int,
    rows: list[dict[str, object]],
    taus: list[float],
    minimum_raw_count: int,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    for tau in taus:
        selected = [row for row in rows if row["n"] == n and row["tau"] == tau]
        plus = [row for row in selected if int(row["plus_count"]) >= minimum_raw_count]
        minus = [row for row in selected if int(row["minus_count"]) >= minimum_raw_count]
        if plus:
            axes[0].plot(
                [row["A"] for row in plus],
                [math.log(float(row["p_plus_raw"])) for row in plus],
                label=rf"$t={tau:g}$",
            )
        if minus:
            axes[1].plot(
                [row["A"] for row in minus],
                [math.log(float(row["p_minus_raw"])) for row in minus],
                label=rf"$t={tau:g}$",
            )
    axes[0].set_xlabel(r"$A$")
    axes[0].set_ylabel(r"$\log P(J_t\geq A)$")
    axes[1].set_xlabel(r"$A$")
    axes[1].set_ylabel(r"$\log P(J_t\leq-A)$")
    for axis in axes:
        axis.grid(alpha=0.2)
        handles, _ = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / f"action_two_tail_logprob_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"action_two_tail_logprob_n{n}.pdf")
    plt.close(fig)


def plot_negative_probability(
    output_dir: Path,
    n: int,
    rows: list[dict[str, object]],
) -> None:
    selected = [row for row in rows if row["n"] == n]
    if not selected:
        return
    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    for observable, label in [
        ("entropy_rate", r"$\Sigma_t^{\rm m}/t$"),
        ("heat_current", r"$J_E(t)$"),
        ("action_current", r"$J_M(t)$"),
    ]:
        subset = sorted(
            (row for row in selected if row["observable"] == observable),
            key=lambda row: float(row["tau"]),
        )
        if subset:
            probability = np.asarray(
                [float(row["negative_probability"]) for row in subset]
            )
            probability[probability <= 0.0] = np.nan
            ax.semilogy(
                [float(row["tau"]) for row in subset],
                probability,
                "o-",
                label=label,
            )
    ax.set_xlabel(r"averaging time $t$")
    ax.set_ylabel(r"raw probability $P(X_t<0)$")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / f"negative_probability_vs_time_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"negative_probability_vs_time_n{n}.pdf")
    plt.close(fig)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def plot_symmetry_slopes(
    output_dir: Path,
    n: int,
    entropy_rows: list[dict[str, str]],
    heat_rows: list[dict[str, str]],
    action_rows: list[dict[str, str]],
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    specifications = [
        (axes[0], entropy_rows, 1.0, "medium entropy"),
        (axes[1], heat_rows, 0.4, "bath heat current"),
        (axes[2], action_rows, None, "action current"),
    ]
    for axis, rows, target, title in specifications:
        selected = sorted(
            (
                row
                for row in rows
                if int(row["n"]) == n and math.isfinite(float(row["ft_slope"]))
            ),
            key=lambda row: float(row["tau"]),
        )
        if selected:
            tau = np.asarray([float(row["tau"]) for row in selected])
            slope = np.asarray([float(row["ft_slope"]) for row in selected])
            low = np.asarray([float(row["ft_slope_ci_low"]) for row in selected])
            high = np.asarray([float(row["ft_slope_ci_high"]) for row in selected])
            errors = np.vstack([slope - low, high - slope])
            errors[~np.isfinite(errors)] = 0.0
            axis.errorbar(1.0 / tau, slope, yerr=errors, fmt="o-", capsize=3)
        if target is not None:
            axis.axhline(target, color="k", linestyle="--", label=f"reference {target:g}")
            axis.legend(fontsize=8)
        axis.set_title(title)
        axis.set_xlabel(r"$1/t$")
        axis.set_ylabel("raw symmetry slope")
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / f"all_symmetry_slopes_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"all_symmetry_slopes_n{n}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    requested_taus = [float(value) for value in args.taus.split(",") if value]
    all_survival_rows: list[dict[str, object]] = []
    all_fit_rows: list[dict[str, object]] = []
    negative_rows: list[dict[str, object]] = []
    chain_lengths: list[int] = []

    for path in args.blocks:
        metadata, raw = load_run(path)
        n = int(metadata["n"])
        chain_lengths.append(n)
        base_tau = float(metadata["block_time"])
        run_rows: list[dict[str, object]] = []
        available_taus: list[float] = []
        for tau in requested_taus:
            factor_float = tau / base_tau
            factor = int(round(factor_float))
            if factor <= 0 or not math.isclose(factor_float, factor, abs_tol=1.0e-12):
                continue
            if factor > raw.shape[1]:
                continue
            aggregated = aggregate_blocks(raw, factor, base_tau)
            available_taus.append(tau)
            rows = survival_rows(
                aggregated["action_current"], n, tau, args.threshold_step
            )
            run_rows.extend(rows)
            all_survival_rows.extend(rows)
            all_fit_rows.append(
                {
                    "source": str(path),
                    "n": n,
                    "tau": tau,
                    **tail_fit_metrics(rows, args.minimum_raw_count),
                }
            )
            for observable in ["entropy_rate", "heat_current", "action_current"]:
                flattened = aggregated[observable].ravel()
                negative_rows.append(
                    {
                        "source": str(path),
                        "n": n,
                        "tau": tau,
                        "observable": observable,
                        "n_samples": flattened.size,
                        "negative_count": int(np.count_nonzero(flattened < 0.0)),
                        "negative_probability": float(np.mean(flattened < 0.0)),
                    }
                )
        plot_two_tail_survival(
            args.output_dir,
            n,
            run_rows,
            available_taus,
            args.minimum_raw_count,
        )

    write_rows(args.output_dir / "action_two_tail_survival.csv", all_survival_rows)
    write_rows(args.output_dir / "action_normal_tail_fit_metrics.csv", all_fit_rows)
    write_rows(args.output_dir / "negative_probability_vs_time.csv", negative_rows)

    analysis_dir = args.analysis_dir or (args.blocks[0].parent / "analysis")
    if analysis_dir.exists():
        entropy_path = analysis_dir / "ft_summary.csv"
        heat_path = analysis_dir / "heat_symmetry_summary.csv"
        action_path = analysis_dir / "action_symmetry_summary.csv"
        if entropy_path.exists() and heat_path.exists() and action_path.exists():
            entropy_rows = read_csv(entropy_path)
            heat_rows = read_csv(heat_path)
            action_rows = read_csv(action_path)
            for n in chain_lengths:
                plot_symmetry_slopes(
                    args.output_dir, n, entropy_rows, heat_rows, action_rows
                )
    for n in chain_lengths:
        plot_negative_probability(args.output_dir, n, negative_rows)
    print(f"wrote supplementary analysis to {args.output_dir}")


if __name__ == "__main__":
    main()
