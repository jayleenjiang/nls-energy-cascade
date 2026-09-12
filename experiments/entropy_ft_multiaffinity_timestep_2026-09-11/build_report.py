#!/usr/bin/env python3
"""Build the frozen multi-affinity timestep audit report."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_key_values(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def number(value, digits: int = 6) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value).replace("_", r"\_")
    if not math.isfinite(numeric):
        return "--"
    if numeric == 0:
        return "0"
    if abs(numeric) < 1e-3 or abs(numeric) >= 1e4:
        return f"{numeric:.{digits}e}"
    return f"{numeric:.{digits}f}"


def yes(value) -> str:
    return "yes" if str(value).lower() in {"1", "true"} else "no"


def label(value: str) -> str:
    return value.replace("_", r"\_")


def rows_tex(rows: list[list[str]]) -> str:
    return "\n".join(" & ".join(row) + r" \\" for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_dir.resolve()
    analysis = root / "analysis"
    report = root / "report"
    report.mkdir(parents=True, exist_ok=True)

    per_tau = read_csv(analysis / "per_tau_results.csv")
    extrap = read_csv(analysis / "extrapolation_summary.csv")
    contrasts = read_csv(analysis / "timestep_contrasts.csv")
    dt_zero = read_csv(analysis / "dt_zero_extrapolation.csv")
    eq_mean = read_csv(analysis / "equilibrium_weighted_baseline.csv")
    eq_diff = read_csv(analysis / "equilibrium_baseline_difference.csv")[0]
    three = read_csv(analysis / "three_affinity_finest_summary.csv")
    hashes = read_csv(analysis / "input_hashes.csv")
    verdict = json.loads((analysis / "VERDICT.json").read_text())
    manifest = json.loads((analysis / "analysis_manifest.json").read_text())
    production = read_key_values(root / "production_manifest.txt")

    names = {"weak": r"$(6.5,5.5)$", "moderate": r"$(8,4)$", "mid": r"$(7,5)$"}
    per_tables = {}
    extrap_tables = {}
    contrast_tables = {}
    for affinity in ("weak", "moderate"):
        selected = [row for row in per_tau if row["affinity"] == affinity]
        selected.sort(key=lambda row: (-float(row["dt"]), int(row["tau"])))
        per_tables[affinity] = rows_tex([[
            number(row["dt"], 3), row["tau"], row["N_windows"], row["n_negative"],
            row["n_symmetric_pairs"], number(row["a_fit"]),
            f"[{number(row['a_ci_low'])}, {number(row['a_ci_high'])}]",
            number(row["R2"], 4), yes(row["admissible"]),
        ] for row in selected])

        selected_extrap = [row for row in extrap if row["affinity"] == affinity]
        selected_extrap.sort(key=lambda row: (-float(row["dt"]), row["selection"]))
        extrap_tables[affinity] = rows_tex([[
            number(row["dt"], 3), label(row["selection"]), row["times_used"],
            number(row["a_inf"]),
            f"[{number(row['a_inf_ci_low'])}, {number(row['a_inf_ci_high'])}]",
            number(row["a_inf_over_delta_beta"]),
            f"[{number(row['ratio_ci_low'])}, {number(row['ratio_ci_high'])}]",
            number(row["linear_reduced_chi2"], 3),
            number(row["max_abs_standardized_residual"], 3),
            yes(row["linear_adequate"]),
        ] for row in selected_extrap])

        selected_contrasts = [
            row for row in contrasts
            if row["affinity"] == affinity and row["quantity"] == "per_tau_slope"
        ]
        selected_contrasts.sort(key=lambda row: (-float(row["finer_dt"]), int(row["tau"])))
        contrast_tables[affinity] = rows_tex([[
            number(row["finer_dt"], 3), row["tau"],
            number(row["finer_minus_baseline"]),
            f"[{number(row['difference_ci_low'])}, {number(row['difference_ci_high'])}]",
            yes(row["difference_ci_contains_zero"]),
            number(row["mean_Q_relative_change"], 6),
        ] for row in selected_contrasts])

    dt_table = rows_tex([[
        label(row["affinity"]), label(row["selection"]),
        number(row["dt_zero_ratio"]),
        f"[{number(row['dt_zero_ci_low'])}, {number(row['dt_zero_ci_high'])}]",
        number(row["dt_slope"], 3),
        f"[{number(row['dt_slope_ci_low'], 3)}, {number(row['dt_slope_ci_high'], 3)}]",
        number(row["linear_reduced_chi2"], 3), yes(row["linear_Odt_compatible"]),
        yes(row["all_finite_time_fits_adequate"]), yes(row["same_times_used"]),
    ] for row in dt_zero])

    eq_tau = [row for row in per_tau if row["affinity"] == "equilibrium"]
    eq_tau.sort(key=lambda row: (-float(row["dt"]), int(row["tau"])))
    eq_tau_table = rows_tex([[
        number(row["dt"], 3), row["tau"], number(row["a_fit"]),
        f"[{number(row['a_ci_low'])}, {number(row['a_ci_high'])}]",
        yes(row["ci_contains_reference"]), row["n_negative"],
        row["n_symmetric_pairs"], number(row["R2"], 4),
    ] for row in eq_tau])
    eq_mean_table = rows_tex([[
        number(row["dt"], 3), number(row["weighted_mean_slope"]),
        number(row["bootstrap_se"]),
        f"[{number(row['ci_low'])}, {number(row['ci_high'])}]",
        number(row["z_from_zero"], 3),
        f"{row['positive_windows']}/{row['negative_windows']}/{row['zero_windows']}",
        row["sign_pattern_tau_order"],
        number(row["fraction_of_dbeta_0p057143"], 4),
    ] for row in eq_mean])

    finest_table = rows_tex([[
        names[row["affinity"]], label(row["selection"]), number(row["dt"], 3),
        number(row["a_inf_over_delta_beta"]),
        f"[{number(row['ratio_ci_low'])}, {number(row['ratio_ci_high'])}]",
        number(row["linear_reduced_chi2"], 3),
        number(row["max_abs_standardized_residual"], 3),
        yes(row["adequate"]), yes(row["contains_one"]), yes(row["finest_step_pass"]),
    ] for row in sorted(three, key=lambda row: (row["selection"], row["affinity"]))])

    provenance_table = rows_tex([[
        label(row["case"]), label(row["kind"]), row["rows"],
        r"\texttt{" + row["sha256"][:16] + r"\ldots}",
    ] for row in hashes])

    primary = [row for row in three if row["selection"] == "all_admissible"]
    weak = next(row for row in primary if row["affinity"] == "weak")
    moderate = next(row for row in primary if row["affinity"] == "moderate")
    primary_text = (
        rf"The predeclared verdict is \texttt{{{label(verdict['verdict'])}}}. "
        rf"At the finest timestep, $(6.5,5.5)$ has "
        rf"$a_\infty/\Delta\beta={number(weak['a_inf_over_delta_beta'])}$ "
        rf"with 95\% CI [{number(weak['ratio_ci_low'])}, {number(weak['ratio_ci_high'])}]. "
        rf"For $(8,4)$ the finest primary fit "
        rf"{'passes' if yes(moderate['adequate']) == 'yes' else 'fails'} the frozen adequacy gate. "
        rf"The nominal equilibrium estimator floor is "
        rf"{number(verdict['nominal_estimator_floor_absolute'])} in slope units "
        rf"({number(verdict['nominal_estimator_floor_fraction_of_dbeta_0p057143'], 4)} "
        rf"of $\Delta\beta=0.057143$)."
    )

    tex = rf"""\documentclass[10pt]{{article}}
\usepackage[margin=0.68in]{{geometry}}
\usepackage{{amsmath,amssymb,booktabs,longtable,array,graphicx,xcolor,hyperref,pdflscape}}
\hypersetup{{colorlinks=true,linkcolor=blue,urlcolor=blue}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.4em}}
\newcommand{{\Db}}{{\Delta\beta}}
\title{{Multi-affinity Timestep and Estimator-floor Audit}}
\author{{Numerical audit report}}
\date{{11 September 2026}}
\begin{{document}}
\sloppy
\maketitle

\section{{Primary answer}}
{primary_text}

Every new case uses the same validated sampler, 128 independent streams,
31,252 duration-5 blocks per stream, and 4,000,256 raw rows.  Both requested
extra affinities were run at $dt=2.5\times10^{{-4}}$ and
$1.25\times10^{{-4}}$; the $dt=5\times10^{{-4}}$ trajectories are immutable
baselines.

\begin{{figure}}[ht]
\centering
\includegraphics[width=0.98\linewidth]{{../figures/timestep_summary.pdf}}
\caption{{Left: finite-time intercept ratios at three timesteps. Right:
relative mean-heat shifts at $\tau=5$ and 25; colours encode timestep.}}
\end{{figure}}

\section{{Finest-step three-affinity summary}}
\scriptsize
\begin{{longtable}}{{l l r r l r r c c c}}
\toprule
affinity & selection & $dt$ & $a_\infty/\Db$ & 95\% CI & $\chi^2_\nu$ & max $|r|$ & adequate & CI has 1 & pass\\
\midrule
\endfirsthead
{finest_table}
\bottomrule
\end{{longtable}}
\normalsize

The direct finest-step criterion requires both an adequate $1/\tau$ model and
a ratio CI containing one.  An inadequate $(8,4)$ fit is reported as
unresolved rather than converted into evidence against the fluctuation
relation.

\section{{Per-window results at $(6.5,5.5)$}}
\begin{{figure}}[ht]
\centering
\includegraphics[width=0.72\linewidth]{{../figures/finite_time_weak.pdf}}
\caption{{Weak-affinity finite-time slopes and frozen $1/\tau$ fits.}}
\end{{figure}}
\scriptsize
\begin{{longtable}}{{r r r r r r l r c}}
\toprule
$dt$ & $\tau$ & $N$ & $n_-$ & pairs & $a(\tau)$ & bootstrap 95\% CI & $R^2$ & adm.\\
\midrule
\endfirsthead
{per_tables['weak']}
\bottomrule
\end{{longtable}}
\normalsize

\section{{Per-window results at $(8,4)$}}
\begin{{figure}}[ht]
\centering
\includegraphics[width=0.72\linewidth]{{../figures/finite_time_moderate.pdf}}
\caption{{Moderate-affinity finite-time slopes and frozen $1/\tau$ fits.}}
\end{{figure}}
\scriptsize
\begin{{longtable}}{{r r r r r r l r c}}
\toprule
$dt$ & $\tau$ & $N$ & $n_-$ & pairs & $a(\tau)$ & bootstrap 95\% CI & $R^2$ & adm.\\
\midrule
\endfirsthead
{per_tables['moderate']}
\bottomrule
\end{{longtable}}
\normalsize

\section{{Finite-time extrapolations}}
\begin{{landscape}}
\scriptsize
\subsection*{{$(6.5,5.5)$}}
\begin{{longtable}}{{r l l r l r l r r c}}
\toprule
$dt$ & selection & $\tau$ used & $a_\infty$ & 95\% CI & ratio & ratio CI & $\chi^2_\nu$ & max $|r|$ & adequate\\
\midrule
\endfirsthead
{extrap_tables['weak']}
\bottomrule
\end{{longtable}}
\subsection*{{$(8,4)$}}
\begin{{longtable}}{{r l l r l r l r r c}}
\toprule
$dt$ & selection & $\tau$ used & $a_\infty$ & 95\% CI & ratio & ratio CI & $\chi^2_\nu$ & max $|r|$ & adequate\\
\midrule
\endfirsthead
{extrap_tables['moderate']}
\bottomrule
\end{{longtable}}

\subsection*{{$dt\to0$ diagnostics}}
\begin{{longtable}}{{l l r l r l r c c c}}
\toprule
affinity & selection & $r_0$ & 95\% CI & slope & slope CI & $\chi^2_\nu$ & $O(dt)$ & all adequate & same $\tau$\\
\midrule
\endfirsthead
{dt_table}
\bottomrule
\end{{longtable}}
\normalsize
\end{{landscape}}

The three-point $dt\to0$ lines are diagnostics.  The direct finest-step CI
and the frozen adequacy gate take priority.

\section{{Mean-heat and slope shifts from $dt=5\times10^{{-4}}$}}
\scriptsize
\subsection*{{$(6.5,5.5)$}}
\begin{{longtable}}{{r r r l c r}}
\toprule
finer $dt$ & $\tau$ & $\Delta a$ & bootstrap 95\% CI & contains 0 & relative $\Delta\langle Q\rangle$\\
\midrule
\endfirsthead
{contrast_tables['weak']}
\bottomrule
\end{{longtable}}
\subsection*{{$(8,4)$}}
\begin{{longtable}}{{r r r l c r}}
\toprule
finer $dt$ & $\tau$ & $\Delta a$ & bootstrap 95\% CI & contains 0 & relative $\Delta\langle Q\rangle$\\
\midrule
\endfirsthead
{contrast_tables['moderate']}
\bottomrule
\end{{longtable}}
\normalsize

\section{{Equilibrium estimator baseline}}
\begin{{figure}}[ht]
\centering
\includegraphics[width=0.72\linewidth]{{../figures/equilibrium_estimator_floor.pdf}}
\caption{{Equilibrium per-window asymmetry slopes and their frozen weighted
means.  The exact physical reference is zero.}}
\end{{figure}}

\scriptsize
\begin{{longtable}}{{r r r l c r r r}}
\toprule
$dt$ & $\tau$ & $a(\tau)$ & bootstrap 95\% CI & contains 0 & $n_-$ & pairs & $R^2$\\
\midrule
\endfirsthead
{eq_tau_table}
\bottomrule
\end{{longtable}}

\begin{{longtable}}{{r r r l r l l r}}
\toprule
$dt$ & weighted mean & boot. SE & 95\% CI & $z$ & $+/-/0$ & signs ($\tau$ order) & fraction of 0.057143\\
\midrule
\endfirsthead
{eq_mean_table}
\bottomrule
\end{{longtable}}
\normalsize

The fine-minus-coarse equilibrium weighted baseline is
{number(eq_diff['fine_minus_coarse'])}, with bootstrap SE
{number(eq_diff['bootstrap_se'])} and 95\% CI
[{number(eq_diff['ci_low'])}, {number(eq_diff['ci_high'])}].  The two baselines
are {'statistically consistent' if yes(eq_diff['ci_contains_zero']) == 'yes' else 'statistically different'} under the frozen criterion.
The nominal floor is {number(verdict['nominal_estimator_floor_absolute'])}; its
ratio to the residual $(7,5)$ finest-step deficit is
{number(verdict['floor_to_mid_residual_ratio'], 3)}, which is
{'within' if verdict['same_order_by_predeclared_factor_three_rule'] else 'outside'}
the predeclared factor-three ``same order'' range.

\section{{Protocol and provenance}}
The freeze/launch commit is
\path{{{production.get('repository_commit_at_launch', 'not-recorded')}}}.\par
Production source SHA-256:\par
\path{{98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3}}.\par
Binary SHA-256:\par
\path{{4c4880d721733897d200f2601690da873f5b226aaa1013df23df2351c1dfe7d1}}.\par
Analysis source SHA-256:\par
\path{{{manifest['analysis_source_sha256']}}}.\par

\begin{{tabular}}{{l l r l}}
\toprule
case & kind & rows & SHA-256\\
\midrule
{provenance_table}
\bottomrule
\end{{tabular}}

Exact commands and measured runtimes are in \texttt{{COMMANDS.tsv}} and
\texttt{{RUNTIMES.tsv}}.  Complete per-window moments, raw symmetric-bin
counts, all bootstrap draws, residuals, and hashes accompany this report.  No
scientific setting or gate was changed after production began.

\end{{document}}
"""
    (report / "multiaffinity_timestep_report.tex").write_text(tex)


if __name__ == "__main__":
    main()

