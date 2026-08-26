#!/usr/bin/env python3
"""Build a compact LaTeX report from an audited entropy/action FT run."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--supplement-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--status-label", default="production")
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def fmt(value: str | float, digits: int = 4) -> str:
    number = float(value)
    if not math.isfinite(number):
        return "--"
    if number == 0.0:
        return "0"
    if abs(number) < 1.0e-3 or abs(number) >= 1.0e4:
        return f"{number:.2e}"
    return f"{number:.{digits}f}"


def relative_path(path: Path, output_dir: Path) -> str:
    import os

    return Path(os.path.relpath(path.resolve(), output_dir.resolve())).as_posix()


def reference_count(
    rows: list[dict[str, str]], target: float
) -> tuple[int, int]:
    usable = [
        row
        for row in rows
        if finite(row["ft_slope"])
        and finite(row["ft_slope_ci_low"])
        and finite(row["ft_slope_ci_high"])
        and int(row["ft_bins_used"]) >= 3
    ]
    covering = [
        row
        for row in usable
        if float(row["ft_slope_ci_low"]) <= target <= float(row["ft_slope_ci_high"])
    ]
    return len(usable), len(covering)


def figure_block(path: Path, output_dir: Path, caption: str) -> str:
    if not path.exists():
        return ""
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=0.98\textwidth]{{{relative_path(path, output_dir)}}}
\caption{{{caption}}}
\end{{figure}}
"""


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for n in [10, 20, 30, 40]:
        path = args.run_dir / f"n{n}_summary.csv"
        if path.exists():
            rows = read_rows(path)
            if len(rows) == 1:
                summaries.append(rows[0])
    if len(summaries) != 4:
        raise FileNotFoundError("all four n=10,20,30,40 summary files are required")

    entropy = read_rows(args.analysis_dir / "ft_summary.csv")
    heat = read_rows(args.analysis_dir / "heat_symmetry_summary.csv")
    action = read_rows(args.analysis_dir / "action_symmetry_summary.csv")
    coupling = read_rows(args.analysis_dir / "coupling_summary.csv")
    tails = read_rows(args.supplement_dir / "action_normal_tail_fit_metrics.csv")
    negative = read_rows(args.supplement_dir / "negative_probability_vs_time.csv")

    entropy_usable, entropy_covering = reference_count(entropy, 1.0)
    heat_usable, heat_covering = reference_count(heat, 0.4)

    conclusion = (
        "The direct finite-time medium-entropy diagnostic is consistent with its "
        "unit-slope reference in "
        f"{entropy_covering} of {entropy_usable} raw-count-qualified fits; the bath-heat "
        f"diagnostic contains its $\\Delta\\beta=0.4$ reference in {heat_covering} of "
        f"{heat_usable} qualified fits.  These counts are a resolution summary, not a "
        "binary theorem test.  The scientific interpretation must also use the time trend, "
        "the fitted range, and the omitted system-entropy boundary term."
    )

    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=1in]{geometry}",
        r"\usepackage{amsmath,amssymb,booktabs,graphicx,float}",
        r"\usepackage[hidelinks]{hyperref}",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{0.55em}",
        rf"\title{{Bath-heat and action-current fluctuation study ({args.status_label})}}",
        r"\author{}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Question and scope}",
        r"We test two distinct questions using the same projection-free Cartesian trajectory. "
        r"The thermodynamic diagnostic uses the bath heats and medium entropy, while the action "
        r"current is analyzed as a separate finite-time observable.  Correlation between the two "
        r"does not by itself identify action current with entropy production.",
        r"The recorded quantities are",
        r"\[J_E(t)=\frac{Q_L-Q_R}{2t},\qquad "
        r"\Sigma_t^{\rm m}=-\frac{Q_L}{T_L}-\frac{Q_R}{T_R},\qquad "
        r"J_M(t)=t^{-1}\int_0^t j_M(s)\,ds.\]",
        r"For $T_L=10$ and $T_R=2$, the heat-current affinity is "
        r"$\Delta\beta=T_R^{-1}-T_L^{-1}=0.4$.  At finite time,",
        r"\[\frac{\Sigma_t^{\rm m}}t=\Delta\beta J_E(t)-"
        r"\frac12(T_L^{-1}+T_R^{-1})\frac{\Delta E}{t}+"
        r"O(\text{integration error}).\]",
        r"\section*{Executive result}",
        conclusion,
        r"Plus-four tail probabilities are used only to display zero or very small counts. "
        r"All symmetry fits below require nonzero raw counts on both sides.",
        r"\section*{Simulation and numerical audit}",
        rf"The {args.status_label} sampler uses $c_j=x_j+i y_j$, $\gamma=0.1$, "
        r"$dt=5\times10^{-4}$, burn-in $B=500$, and non-overlapping base blocks of "
        r"length $t=20$.  There is no action-floor projection.  Each block records "
        r"$Q_L,Q_R,\Delta E,\Sigma_t^{\rm m}$, the action current, stream ID, and block ID.",
        r"\begin{table}[H]",
        r"\centering",
        r"\begin{tabular}{rrrrrr}",
        r"\toprule",
        r"$n$ & blocks & $\langle J_M\rangle$ & $\langle J_E\rangle$ & "
        r"$\langle\Sigma^{\rm m}/t\rangle$ & balance RMS \\ ",
        r"\midrule",
    ]
    for row in summaries:
        heat_mean = 0.5 * (
            float(row["mean_q_left_rate"]) - float(row["mean_q_right_rate"])
        )
        lines.append(
            f"{row['n']} & {int(row['n_blocks']):,} & {fmt(row['mean_action_current'], 5)} & "
            f"{fmt(heat_mean, 5)} & {fmt(row['mean_entropy_rate'], 5)} & "
            f"{fmt(row['rms_energy_balance_error_rate'], 3)} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Production means and discrete first-law residual.}",
            r"\end{table}",
            r"\section*{Medium-entropy and heat-current symmetry}",
            r"For the entropy-production rate $a=\Sigma_t^{\rm m}/t$, the plotted "
            r"symmetry function is",
            r"\[R_t(a)=t^{-1}\log\!\left[\frac{p_t(a)}{p_t(-a)}\right].\]",
            r"The simple asymptotic reference has slope one.  For $J_E$, the corresponding "
            r"reference slope is $\Delta\beta=0.4$.",
            r"\par",
            r"\begin{table}[H]",
            r"\centering\small",
            r"\resizebox{\textwidth}{!}{%",
            r"\begin{tabular}{rrrrrrrr}",
            r"\toprule",
            r"$n$ & $t$ & $N$ & negative & slope & 95\% CI & $R^2$ & bins \\ ",
            r"\midrule",
        ]
    )
    selected_taus = {20.0, 40.0, 80.0, 120.0, 160.0, 200.0}
    for row in entropy:
        if float(row["tau"]) not in selected_taus or not finite(row["ft_slope"]):
            continue
        lines.append(
            f"{row['n']} & {fmt(row['tau'], 0)} & {int(row['n_samples']):,} & "
            f"{int(row['negative_count']):,} & {fmt(row['ft_slope'])} & "
            f"[{fmt(row['ft_slope_ci_low'])}, {fmt(row['ft_slope_ci_high'])}] & "
            f"{fmt(row['ft_r_squared'])} & {row['ft_bins_used']} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"}",
            r"\caption{Raw-count-qualified medium-entropy symmetry fits. Missing rows have "
            r"insufficient two-sided counts.}",
            r"\end{table}",
            r"\section*{Action-current tails and symmetry}",
            r"The central 1--99\% histogram is compared with a fitted normal density.  Each "
            r"empirical upper and lower survival tail is also compared with the corresponding "
            r"normal CDF.  The action-current symmetry slope is reported empirically, without "
            r"assigning it the thermodynamic reference value one.",
            r"\par",
            r"\begin{table}[H]",
            r"\centering\small",
            r"\resizebox{\textwidth}{!}{%",
            r"\begin{tabular}{rrrrrrr}",
            r"\toprule",
            r"$n$ & $t$ & $N$ & negative & slope & 95\% CI & bins \\ ",
            r"\midrule",
        ]
    )
    for row in action:
        if float(row["tau"]) not in selected_taus or not finite(row["ft_slope"]):
            continue
        lines.append(
            f"{row['n']} & {fmt(row['tau'], 0)} & {int(row['n_samples']):,} & "
            f"{int(row['negative_count']):,} & {fmt(row['ft_slope'])} & "
            f"[{fmt(row['ft_slope_ci_low'])}, {fmt(row['ft_slope_ci_high'])}] & "
            f"{row['ft_bins_used']} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"}",
            r"\caption{Raw-count-qualified action-current symmetry fits.}",
            r"\end{table}",
            r"\section*{Heat--action coupling}",
            r"Tight coupling would require more than a large Pearson correlation: the "
            r"heat-on-action regression must be stable with $t$, and the residual variance "
            r"fraction must become small.  The joint plots and tables report all three quantities.",
            r"\par",
            r"\begin{table}[H]",
            r"\centering\small",
            r"\begin{tabular}{rrrrr}",
            r"\toprule",
            r"$n$ & $t$ & Pearson $r$ & heat-on-action slope & residual fraction \\ ",
            r"\midrule",
        ]
    )

    for row in coupling:
        if float(row["tau"]) not in {20.0, 200.0}:
            continue
        lines.append(
            f"{row['n']} & {fmt(row['tau'], 0)} & {fmt(row['pearson_correlation'])} & "
            f"{fmt(row['heat_on_action_slope'])} & "
            f"{fmt(row['residual_variance_fraction'])} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Heat--action coupling at the shortest and longest requested windows.}",
            r"\end{table}",
        ]
    )

    for n in [10, 20, 30, 40]:
        lines.append(
            figure_block(
                args.analysis_dir / f"entropy_pdf_symmetry_n{n}.pdf",
                args.output_dir,
                f"Medium-entropy PDF and raw/plus-four symmetry diagnostic for $n={n}$.",
            )
        )
        lines.append(
            figure_block(
                args.analysis_dir / f"action_tail_normal_fit_n{n}_t20.pdf",
                args.output_dir,
                f"Central action-current histogram and both normal-tail benchmarks for $n={n}$ at $t=20$.",
            )
        )
        lines.append(
            figure_block(
                args.supplement_dir / f"action_two_tail_logprob_n{n}.pdf",
                args.output_dir,
                f"Raw two-tail log probabilities versus threshold and averaging time for $n={n}$.",
            )
        )
        lines.append(
            figure_block(
                args.analysis_dir / f"heat_action_coupling_n{n}.pdf",
                args.output_dir,
                f"Joint bath-heat/action-current sample for $n={n}$.",
            )
        )

    lines.extend(
        [
            r"\clearpage",
            r"\section*{Interpretation guardrail}",
            r"A finite-time detailed FT is a statement about total trajectory entropy, "
            r"$\Delta s_{\rm tot}=\Sigma_t^{\rm m}+\Delta s_{\rm sys}$.  The unknown "
            r"NESS system-entropy endpoint term is not reconstructed here.  Consequently, a "
            r"deviation of the medium-entropy slope from one is reported as ``not verified in "
            r"the sampled window,'' not as a universal violation.  Conversely, visual linearity "
            r"or a Gaussian-looking center is not by itself evidence for the FT.",
            r"\end{document}",
        ]
    )

    output = args.output_dir / "entropy_ft_report.tex"
    output.write_text("\n".join(lines) + "\n")
    print(output)


if __name__ == "__main__":
    main()
