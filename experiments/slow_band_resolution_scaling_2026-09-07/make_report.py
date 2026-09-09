#!/usr/bin/env python3
"""Create figures, verdict files, and the Part-1 LaTeX report."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
A = HERE / "part1_analysis"
F = HERE / "figures"
R = HERE / "report"
F.mkdir(exist_ok=True)
R.mkdir(exist_ok=True)

CASES = ("driven_dt2p5e-4", "equilibrium_dt2p5e-4")
CASE_LABEL = {
    "driven_dt2p5e-4": "Driven (2,8)",
    "equilibrium_dt2p5e-4": "Equilibrium (5,5)",
}
STREAM_COUNTS = (64, 128, 256, 512, 1024, 2048)
RULES = (
    ("fixed", "signed", "fixed / signed", "#0072B2", "o"),
    ("fixed", "nonnegative", "fixed / nonnegative", "#D55E00", "s"),
    ("extended", "signed", "extended / signed", "#009E73", "^"),
    ("extended", "nonnegative", "extended / nonnegative", "#CC79A7", "D"),
)


def read(name: str) -> list[dict[str, str]]:
    with (A / name).open(newline="") as f:
        return list(csv.DictReader(f))


def fmt(value: object, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return r"$\infty$" if number > 0 else "N/A"
    return f"{number:.{digits}f}"


def make_figures() -> None:
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "figure.dpi": 150,
    })
    curves = read("resolution_curve.csv")
    summary = read("synthetic_scaling_summary.csv")

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.25), sharex=True, sharey=True)
    for ax, case in zip(axes, CASES):
        for window, constraint, label, color, marker in RULES:
            rows = sorted(
                (row for row in curves if row["case"] == case
                 and row["window"] == window and row["constraint"] == constraint),
                key=lambda row: int(row["N_streams"]),
            )
            y = [
                0.68 if row["smallest_resolved_separation"].startswith(">")
                else float(row["smallest_resolved_separation"])
                for row in rows
            ]
            ax.plot(STREAM_COUNTS, y, marker=marker, ms=4, lw=1.2,
                    color=color, label=label)
        ax.axhline(0.20, color="black", lw=1, ls="--", label=r"target $\Delta=0.20$")
        ax.set_xscale("log", base=2)
        ax.set_xticks(STREAM_COUNTS, [str(n) for n in STREAM_COUNTS], rotation=35)
        ax.set_ylim(0.15, 0.71)
        ax.set_yticks((0.2, 0.3, 0.4, 0.5, 0.6, 0.68),
                      ("0.2", "0.3", "0.4", "0.5", "0.6", ">0.6"))
        ax.set_title(CASE_LABEL[case])
        ax.set_xlabel("emulated independent streams $N$")
    axes[0].set_ylabel("smallest tested resolved separation")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.08), fontsize=8)
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(F / f"resolution_curves.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.25), sharex=True, sharey=True)
    for ax, case in zip(axes, CASES):
        for window, constraint, label, color, marker in RULES:
            rows = sorted(
                (row for row in summary if row["case"] == case
                 and abs(float(row["separation"]) - 0.20) < 1e-12
                 and row["window"] == window and row["constraint"] == constraint),
                key=lambda row: int(row["N_streams"]),
            )
            errors = [max(float(row["median_abs_error_slow"]),
                          float(row["median_abs_error_fast"])) for row in rows]
            ax.plot(STREAM_COUNTS, errors, marker=marker, ms=4, lw=1.2,
                    color=color, label=label)
        benchmark = 0.20 * np.sqrt(256 / np.asarray(STREAM_COUNTS, float))
        ax.plot(STREAM_COUNTS, benchmark, color="0.35", ls=":", lw=1.1,
                label=r"$N^{-1/2}$ reference")
        ax.axhline(0.05, color="black", lw=1, ls="--", label="frozen tolerance")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks(STREAM_COUNTS, [str(n) for n in STREAM_COUNTS], rotation=35)
        ax.set_title(CASE_LABEL[case])
        ax.set_xlabel("emulated independent streams $N$")
    axes[0].set_ylabel(r"max. median rate error at $\Delta=0.20$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.08), fontsize=8)
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(F / f"delta020_error_scaling.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.8), sharex=True)
    for column, case in enumerate(CASES):
        for window, constraint, label, color, marker in RULES:
            rows = sorted(
                (row for row in summary if row["case"] == case
                 and abs(float(row["separation"]) - 0.20) < 1e-12
                 and row["window"] == window and row["constraint"] == constraint),
                key=lambda row: int(row["N_streams"]),
            )
            axes[0, column].plot(
                STREAM_COUNTS, [float(row["fraction_select_m_ge_2"]) for row in rows],
                marker=marker, ms=3.5, lw=1.1, color=color, label=label,
            )
            axes[1, column].plot(
                STREAM_COUNTS,
                [float(row["fraction_numerical_identifiability_pass"]) for row in rows],
                marker=marker, ms=3.5, lw=1.1, color=color,
            )
        axes[0, column].axhline(0.90, color="black", lw=1, ls="--")
        axes[0, column].set_title(CASE_LABEL[case])
        for row in range(2):
            axes[row, column].set_xscale("log", base=2)
            axes[row, column].set_xticks(STREAM_COUNTS,
                                         [str(n) for n in STREAM_COUNTS], rotation=35)
            axes[row, column].set_ylim(-0.03, 1.03)
        axes[1, column].set_xlabel("emulated independent streams $N$")
    axes[0, 0].set_ylabel(r"fraction selecting $m\geq2$")
    axes[1, 0].set_ylabel("numerical-identifiability fraction")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 1.02), fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for suffix in ("pdf", "png"):
        fig.savefig(F / f"delta020_selection_identifiability.{suffix}",
                    dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_markdown() -> None:
    summary = read("synthetic_scaling_summary.csv")
    extrap = read("delta020_extrapolation.csv")
    gate = json.loads((HERE / "PART1_GATE.json").read_text())
    passes = [row for row in summary if row["resolution_pass"] == "1"]
    target = [row for row in summary if abs(float(row["separation"]) - 0.20) < 1e-12]
    saturated = sum(row["saturated_by_frozen_rule"] == "1" for row in extrap)
    verdict = f"""# Final verdict

Part 1 completed all **336 frozen cells** and **168,000 raw fit replicates**.
The requested separation, `delta=0.20`, failed the frozen resolution rule in
all {len(target)} tested cells, including both bath cases, all six stream
counts through `N=2048`, both windows, and both amplitude constraints.

The Part-2 gate is therefore **closed**: `{gate['reason']}`. No real trajectory
production was launched. This is the protocol-mandated result, not an
interrupted or incomplete run.

There were {len(passes)} passing cells elsewhere on the grid. They occur almost
entirely at the much wider separation `delta=0.60`; the driven fixed,
non-negative rule reaches `delta=0.40` only at `N=2048`. Thus the current
inverse problem can resolve coarse spectral splitting, but not the physically
relevant approximately 0.20 slow-band spacing.

For `delta=0.20`, {saturated} of 8 rule/case error curves meet the frozen
saturation definition. The fitted error exponents are much shallower than the
independent-stream `N^(-1/2)` benchmark and the extrapolated stream counts are
unstable, ranging from about `3.5e4` to effectively infinite. These are
descriptive extrapolations only and cannot authorize Part 2.

The defensible claim is a numerical resolution bound: with the existing
empirical covariance, lag support, and frozen multi-exponential pipeline,
`delta=0.20` is not recoverable with up to 2048 emulated independent streams.
This does not prove that the physical slow band is absent.
"""
    (HERE / "FINAL_VERDICT.md").write_text(verdict)

    audit = json.loads((HERE / "provenance" / "integrity_audit.json").read_text())
    validation = "# Validation report\n\n"
    validation += f"Overall audit: **{'PASS' if audit['overall_pass'] else 'FAIL'}**.\n\n"
    validation += "## Frozen-output checks\n\n"
    for item in audit["checks"]:
        validation += f"- **{item['check']}**: {'PASS' if item['passed'] else 'FAIL'}; {item['detail']}.\n"
    validation += "\n## Claim boundary\n\n"
    validation += (
        "Part 2 was not run because the predeclared delta=0.20 gate failed. "
        "Power-law stream-count estimates are reported only as descriptive "
        "extrapolations and are not observed resolution results.\n"
    )
    (HERE / "VALIDATION_REPORT.md").write_text(validation)


def write_tex() -> None:
    summary = read("synthetic_scaling_summary.csv")
    curves = read("resolution_curve.csv")
    extrap = read("delta020_extrapolation.csv")
    manifest = json.loads((A / "run_manifest.json").read_text())
    audit = json.loads((HERE / "provenance" / "integrity_audit.json").read_text())
    target = [row for row in summary if abs(float(row["separation"]) - 0.20) < 1e-12]

    lines = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=0.72in]{geometry}",
        r"\usepackage{amsmath,booktabs,longtable,array,graphicx,xcolor,hyperref,seqsplit}",
        r"\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}",
        r"\newcommand{\pass}{\textcolor{green!45!black}{PASS}}",
        r"\newcommand{\fail}{\textcolor{red!70!black}{FAIL}}",
        r"\title{Data-volume calibration for resolving the three-mode slow relaxation band}",
        r"\author{Frozen synthetic-covariance audit}",
        r"\date{9 September 2026}",
        r"\begin{document}",
        r"\maketitle",
        r"\begin{abstract}",
        r"We calibrate how many independent trajectories are required to separate two nearby real relaxation rates using the empirical cross-lag covariance of the controlled three-mode data. The complete frozen grid comprises 336 cells and 168,000 nonlinear-fit replicates. The target spacing $\Delta=0.20$ fails for both driven and equilibrium covariance at every tested stream count through $N=2048$, for both fit windows and both amplitude constraints. The predeclared Part-2 gate is therefore closed and no additional physical trajectory simulation is launched. The result is a quantitative non-resolution bound, not evidence that the physical slow band is absent.",
        r"\end{abstract}",
        r"\section{Question and frozen protocol}",
        r"The synthetic truth is",
        r"\[C(t)=0.6e^{-0.90t}+0.4e^{(-0.90-\Delta)t},\]",
        r"with $\Delta\in\{0.05,0.10,0.15,0.20,0.30,0.40,0.60\}$ and $N\in\{64,128,256,512,1024,2048\}$. Noise preserves the full empirical cross-lag covariance of the saved $I_2$ correlation curves and scales as $\sqrt{64/N}$. We use 500 deterministic replicates per case--$N$--$\Delta$--window--constraint cell.",
        r"Real sums of one to four exponentials are selected sequentially by $\Delta\mathrm{AIC}\leq-10$. A cell resolves the two rates only if at least 90\% of successful fits select $m\geq2$ and both median matched-rate errors are at most $\max(0.05,\Delta/4)$. The fixed window is $[0.20,4.00]$; the extended window follows the frozen contiguous $C/\mathrm{SE}\geq3$ rule. No rule was changed after output existed.",
        r"\section{Primary result: the Part-2 gate is closed}",
        r"None of the 48 cells at the target spacing $\Delta=0.20$ passes. Consequently there is no common tested $N\leq2048$, window, and amplitude rule that resolves the target in both bath cases. Under the frozen protocol, Part 2 is forbidden.",
        r"\begin{figure}[ht]\centering\includegraphics[width=.98\linewidth]{../figures/resolution_curves.pdf}\caption{Smallest tested separation satisfying the complete resolution rule. Values marked $>0.6$ mean no tested separation passed. The dashed line is the target $\Delta=0.20$.}\end{figure}",
        r"\begin{table}[ht]\centering\small\begin{tabular}{lllrr}\toprule Case & Window & amplitudes & smallest resolved at $N=2048$ & resolved cells\\\midrule",
    ]
    for case in CASES:
        for window, constraint, label, _, _ in RULES:
            row = next(row for row in curves if row["case"] == case
                       and int(row["N_streams"]) == 2048
                       and row["window"] == window and row["constraint"] == constraint)
            value = row["smallest_resolved_separation"]
            value_tex = r"$>0.60$" if value.startswith(">") else fmt(value, 2)
            lines.append(
                f"{CASE_LABEL[case]} & {window} & {constraint} & {value_tex} & {row['number_resolved']} " + r"\\"
            )
    lines += [
        r"\bottomrule\end{tabular}\caption{Observed resolution at the largest tested stream count.}\end{table}\clearpage",
        r"\section{Why increasing \texorpdfstring{$N$}{N} did not resolve \texorpdfstring{$\Delta=0.20$}{Delta=0.20}}",
        r"The probability of selecting multiple exponentials approaches one, but the matched rate errors remain well above the 0.05 tolerance. Simultaneously, the original numerical-identifiability fraction generally decreases. Model-order selection alone is therefore not evidence that the two physical rates have been recovered.",
        r"\begin{figure}[ht]\centering\includegraphics[width=.98\linewidth]{../figures/delta020_error_scaling.pdf}\caption{Maximum of the two median matched-rate errors at $\Delta=0.20$. The dotted $N^{-1/2}$ curve is a slope reference, not a fit; the dashed line is the frozen tolerance.}\end{figure}",
        r"\begin{figure}[ht]\centering\includegraphics[width=.98\linewidth]{../figures/delta020_selection_identifiability.pdf}\caption{At $\Delta=0.20$, selection of $m\geq2$ becomes frequent while numerical identifiability deteriorates.}\end{figure}\clearpage",
        r"\section{All target-spacing summaries}",
        r"\scriptsize\begin{longtable}{llrlrrrrr}\toprule Case & Window & amplitudes & $N$ & $P(m\geq2)$ & error slow & error fast & num.-ID & pass\\\midrule\endhead",
    ]
    for row in sorted(target, key=lambda x: (
            x["case"], x["window"], x["constraint"], int(x["N_streams"]))):
        lines.append(
            f"{CASE_LABEL[row['case']]} & {row['window']} & {row['constraint']} & "
            f"{row['N_streams']} & {fmt(row['fraction_select_m_ge_2'])} & "
            f"{fmt(row['median_abs_error_slow'])} & {fmt(row['median_abs_error_fast'])} & "
            f"{fmt(row['fraction_numerical_identifiability_pass'])} & "
            f"{'yes' if row['resolution_pass']=='1' else 'no'} " + r"\\"
        )
    lines += [
        r"\bottomrule\end{longtable}\normalsize",
        r"\section{Scaling, saturation, and descriptive extrapolation}",
        r"For each case and rule, a descriptive power law was fit to the maximum median error over $N\geq256$. Seven of eight curves meet the frozen saturation rule: the target remains unresolved at $N=2048$ and the error decreases by less than 10\% from $N=1024$ to 2048. Exponents are much shallower than the independent-stream benchmark $-1/2$. Extrapolated $N$ values do not include the separate 90\% model-selection requirement and cannot authorize production.",
        r"\small\begin{longtable}{lllrrrr}\toprule Case & Window & amplitudes & exponent & implied $N$ & reduction & saturated\\\midrule\endhead",
    ]
    for row in extrap:
        implied = fmt(row["implied_N_from_error_power_law"], 0)
        lines.append(
            f"{CASE_LABEL[row['case']]} & {row['window']} & {row['constraint']} & "
            f"{fmt(row['loglog_error_exponent'])} & {implied} & "
            f"{100*float(row['relative_error_reduction_1024_to_2048']):.1f}\\% & "
            f"{'yes' if row['saturated_by_frozen_rule']=='1' else 'no'} " + r"\\"
        )
    lines += [
        r"\bottomrule\end{longtable}\normalsize",
        r"\section{Audit and reproducibility}",
        f"The integrity audit is {'\\pass' if audit['overall_pass'] else '\\fail'}. It reconstructs all 336 summaries from the 168,000 raw rows, verifies exactly 500 replicate indices in every cell, reproduces the Part-2 gate, and verifies the source and both input hashes.",
        r"\begin{verbatim}",
        r"python3 experiments/slow_band_resolution_scaling_2026-09-07/run_scaling.py",
        r"  --output experiments/slow_band_resolution_scaling_2026-09-07/part1_analysis",
        r"  --replicates 500 --workers 8",
        r"\end{verbatim}",
        f"Master seed: {manifest['master_seed']}. Runtime: {manifest['elapsed_seconds']/3600:.2f} hours. Analysis source SHA-256: \\texttt{{\\seqsplit{{{manifest['analysis_source_sha256']}}}}}.",
        r"\begin{itemize}",
        f"\\item Driven input SHA-256: \\texttt{{\\seqsplit{{{manifest['input_npz_sha256']['driven_dt2p5e-4']}}}}}.",
        f"\\item Equilibrium input SHA-256: \\texttt{{\\seqsplit{{{manifest['input_npz_sha256']['equilibrium_dt2p5e-4']}}}}}.",
        f"\\item Previous fitting pipeline SHA-256: \\texttt{{\\seqsplit{{{manifest['previous_pipeline_sha256']}}}}}.",
        r"\end{itemize}",
        r"All raw replicates, cell summaries, task seeds, resolution curves, extrapolation diagnostics, hashes, and the frozen protocol accompany this report.",
        r"\section{Claim boundary}",
        r"The result supports only the following statement: under the empirical covariance and frozen fitting protocol, a separation of 0.20 is not recoverable with up to 2048 emulated independent streams. It does not establish that the physical slow band is absent, and it does not justify quoting a unique spectral gap. No real-data Part-2 simulation was run because its predeclared gate failed.",
        r"\end{document}",
    ]
    (R / "slow_band_resolution_scaling_report.tex").write_text("\n".join(lines) + "\n")


def main() -> None:
    make_figures()
    write_markdown()
    write_tex()


if __name__ == "__main__":
    main()
