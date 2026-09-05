#!/usr/bin/env python3
"""Build a concise auditable LaTeX report from the frozen profile analysis."""

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"
REPORT = ROOT / "report"


def rows(name):
    with (ANALYSIS / name).open() as handle:
        return list(csv.DictReader(handle))


def val(x, digits=4):
    try:
        number = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not math.isfinite(number):
        return "--"
    return f"{number:.{digits}g}"


def esc(text):
    return str(text).replace("_", r"\_")


def main():
    REPORT.mkdir(parents=True, exist_ok=True)
    run = rows("logical_run_summary.csv")
    eq = rows("equilibrium_checks.csv")
    scaling = rows("midpoint_scaling.csv")
    gates = [json.loads((ANALYSIS / f"bc1_gate_n{n}.json").read_text()) for n in (25, 50, 100)]
    nonfinite_total = sum(int(r["nonfinite_trajectories"]) for r in run)
    discarded_total = sum(int(r["discarded_trajectories"]) for r in run)
    profile_paths = list((ROOT / "curated_results" / "profiles").glob("*_profile.csv"))

    gate_lines = []
    for gate in gates:
        for item in gate["endpoint_checks"]:
            gate_lines.append(
                f'{gate["n"]} & {item["side"]} & {val(item["reference"],7)} & '
                f'{val(item["new_mean"],7)} $\\pm$ {val(item["new_se"],3)} & '
                f'{val(item["signed_discrepancy"],3)} & {val(item["tolerance"],3)} & '
                f'{"PASS" if item["pass"] else "FAIL"} \\\\'
            )

    scaling_lines = []
    for r in scaling:
        scaling_lines.append(
            f'{esc(r["bc_label"])} & {esc(r["temp_label"])} & '
            f'{val(r["loglog_slope_mid_I"],4)} $\\pm$ {val(r["slope_se_three_point"],2)} & '
            f'{val(r["R2_mid_I"],3)} & {val(r["reference_slope"],2)} & '
            f'{esc(r["mid_I_values"])} \\\\'
        )

    eq_lines = []
    for r in eq:
        eq_lines.append(
            f'{esc(r["bc_label"])} & {r["n"]} & {val(r["mean_I_spatial"],4)} & '
            f'{val(r["relative_I_range"],3)} & '
            f'{r["sine_pointwise_95_failures"]} & {val(r["sine_max_abs_z"],3)} & '
            f'{"PASS" if r["sine_simultaneous_pass"] == "True" else "FAIL"} \\\\'
        )

    stability_lines = []
    for r in run:
        stability_lines.append(
            f'{esc(r["bc_label"])} & {esc(r["temp_label"])} & {r["n"]} & '
            f'{r["nonfinite_trajectories"]} & {r["discarded_trajectories"]} & '
            f'{r["projection_count"]} & {val(r["max_seed_z"],3)} & '
            f'{val(r["left_last_quarter_relative_change"],3)} & '
            f'{val(r["right_last_quarter_relative_change"],3)} \\\\'
        )

    tex = rf"""\documentclass[10pt]{{article}}
\usepackage[margin=0.7in]{{geometry}}
\usepackage{{amsmath,booktabs,graphicx,longtable,array}}
\usepackage[T1]{{fontenc}}
\newcommand{{\code}}[1]{{\texttt{{#1}}}}
\title{{Stationary profiles under four boundary couplings}}
\author{{Numerical audit report}}
\date{{2026-09-05}}
\begin{{document}}
\maketitle

\section*{{Frozen design}}
The corrected SIMD implementation evolves 16 trajectories per batch with
$\gamma=0.1$, maximum $dt=5\times10^{{-4}}$, and corrected noises
$\sigma_I=2\sqrt{{2\gamma T I}}$ and
$\sigma_\phi=\sqrt{{2\gamma T/I}}$.  The four couplings are
BC1 $(b_I,b_\phi)=(2\gamma(2T-F),\gamma P)$,
BC2 $(2\gamma(2T-F),0)$,
BC3 $(2\gamma(2T-I^2),0)$, and
BC3b $(2\gamma(2T-I^2),\gamma P)$.
Here $F=2MI_1-I_1^2+2I_1I_2\cos\delta$ and
$P=2I_2\sin\delta$, with the reflected right-boundary expression.
The run matrix contains 36 logical conditions, each formed from two fixed seeds
and 512 trajectories. Burn-ins are 2000, 8000, and 32000 for
$n=25,50,100$; all measurement durations are 2000.

\section*{{BC1 reproduction gate}}
\begin{{center}}\small
\begin{{tabular}}{{rrrllll}}
\toprule
$n$ & side & reference & new mean $\pm$ SE & difference & tolerance & verdict\\
\midrule
{chr(10).join(gate_lines)}
\bottomrule
\end{{tabular}}
\end{{center}}

\section*{{Profiles}}
Shaded bands are pointwise 95\% intervals across independent trajectory time
averages.  Every per-site numerical value is in the accompanying merged CSV.
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/action_profiles_T10_T2.pdf}}
\caption{{Action profiles for $(T_1,T_n)=(10,2)$.}}\end{{figure}}
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/action_profiles_T4_T6.pdf}}
\caption{{Action profiles for $(T_1,T_n)=(4,6)$.}}\end{{figure}}
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/action_profiles_T6_T6.pdf}}
\caption{{Equal-temperature action profiles.}}\end{{figure}}
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/sine_profiles_T10_T2.pdf}}
\caption{{Bond-sine profiles for $(T_1,T_n)=(10,2)$.}}\end{{figure}}
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/sine_profiles_T4_T6.pdf}}
\caption{{Bond-sine profiles for $(T_1,T_n)=(4,6)$.}}\end{{figure}}
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/sine_profiles_T6_T6.pdf}}
\caption{{Equal-temperature bond-sine profiles.}}\end{{figure}}

\clearpage
\section*{{Equal-temperature control}}
The action columns report the spatial mean and relative range.  ``Sine 95\%
failures'' counts individual bonds outside a pointwise interval; the simultaneous
verdict uses a Bonferroni 95\% threshold across all bonds and is not tuned by
case.
\begin{{center}}\small
\begin{{tabular}}{{lrrrrl}}
\toprule BC & $n$ & mean $I$ & relative range & sine 95\% failures & simultaneous sine\\
\midrule
{chr(10).join(eq_lines)}
\bottomrule\end{{tabular}}\end{{center}}

\section*{{Central-profile scaling}}
The exponent is the unconstrained three-point OLS slope of
$\log\langle I\rangle_{{\rm mid}}$ against $\log n$.  The quoted standard error
has one residual degree of freedom and must not be read as an asymptotic
confidence interval.
\begin{{figure}}[ht]\centering
\includegraphics[width=0.98\textwidth]{{../figures/midpoint_scaling.pdf}}
\caption{{Central action versus chain length.}}\end{{figure}}
\begin{{center}}\scriptsize
\begin{{tabular}}{{llrrrr}}
\toprule BC & temperatures & slope $\pm$ SE & $R^2$ & reference & values at 25;50;100\\
\midrule
{chr(10).join(scaling_lines)}
\bottomrule\end{{tabular}}\end{{center}}

\clearpage
\section*{{Stability and reproducibility}}
Endpoint changes compare the cumulative profile at 75\% and 100\% of the
fixed measurement interval.  Seed $z_{{\max}}$ is the largest pointwise
replicate difference divided by its independent-replicate SE over all sites and
bonds; all individual values remain available in the CSVs.
\begin{{longtable}}{{lllrrrrrr}}
\toprule BC & temperatures & $n$ & nonfinite & discarded & projections & seed $z_{{\max}}$ & left change & right change\\
\midrule\endfirsthead
\toprule BC & temperatures & $n$ & nonfinite & discarded & projections & seed $z_{{\max}}$ & left change & right change\\
\midrule\endhead
{chr(10).join(stability_lines)}
\bottomrule
\end{{longtable}}

\section*{{Audit and provenance}}
There are {len(profile_paths)} merged logical-run profile CSVs. Across the run
matrix the analysis found {nonfinite_total} non-finite trajectories and
{discarded_total} discarded trajectories.  The immutable base source SHA-256 is
\code{{3919ab963e9d94bcb25ae5ef1c30c2bb032636525db7f6df4d6a318dd41f0656}}.
The instrumented source SHA-256 is
\code{{40d13c7547d9c147ca8241ed701e1a1e85540fc0c21db48c4a63a943eccfe3ec}},
the binary SHA-256 is
\code{{8d3438d3ea893e59216e9784d2104a33ce00d051d4df54c8c19f814fe2c0e6f4}},
and the integrated protocol commit is
\code{{3e6d8e817ee44c8c686f2d80f7c333f3e9d3a973}}.
Exact commands and seeds are in \code{{COMMANDS.csv}}; file-level hashes are in
\code{{analysis/FILE\_HASHES.csv}}.

\paragraph{{Numerical scope.}}
The requested inherited implementation uses single-precision SIMD dynamics,
an adaptive Euler--Maruyama step with a $10^{{-5}}$ floor, a positive-action
projection, and a polynomial trigonometric approximation. Projection counts
are therefore reported rather than hidden. Claims from unstable BC3/BC3b runs
are limited to the observed instability.
\end{{document}}
"""
    (REPORT / "boundary_profile_report.tex").write_text(tex)


if __name__ == "__main__":
    main()
