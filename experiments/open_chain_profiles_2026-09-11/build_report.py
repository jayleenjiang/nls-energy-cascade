#!/usr/bin/env python3
"""Build the open-chain profile diagnostic LaTeX report."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_key_values(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def num(value, digits: int = 5) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return str(value).replace("_", r"\_")
    if not math.isfinite(value):
        return "--"
    if value == 0:
        return "0"
    if abs(value) < 1e-3 or abs(value) >= 1e4:
        return f"{value:.{digits}e}"
    return f"{value:.{digits}f}"


def yn(value) -> str:
    return "yes" if str(value).lower() in {"1", "true"} else "no"


def tex_rows(rows: list[list[str]]) -> str:
    return "\n".join(" & ".join(row) + r" \\" for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_dir.resolve()
    analysis = root / "analysis"
    report = root / "report"
    report.mkdir(parents=True, exist_ok=True)

    stationarity = read_csv(analysis / "stationarity_gates.csv")
    scaling_points = read_csv(analysis / "scaling_points.csv")
    scaling = read_csv(analysis / "scaling_fits.csv")
    penetration = read_csv(analysis / "penetration_depth.csv")
    spread = read_csv(analysis / "collapse_spread.csv")
    runs = read_csv(analysis / "run_summary.csv")
    verdict = json.loads((analysis / "VERDICT.json").read_text())
    manifest = json.loads((analysis / "analysis_manifest.json").read_text())
    production = read_key_values(root / "production_manifest.txt")

    stationarity.sort(key=lambda row: (row["temperature"], int(row["n"])))
    stationarity_table = tex_rows([[
        row["temperature"], row["n"],
        num(row["delta_mass_q4_minus_q3"]), num(row["se_delta_mass"]),
        num(row["mass_relative_change"], 4), yn(row["mass_ci_contains_zero"]),
        num(row["delta_free_end_q4_minus_q3"]), num(row["se_delta_free_end"]),
        num(row["free_end_relative_change"], 4), yn(row["free_end_ci_contains_zero"]),
        num(row["interior_median_site_relative_change"], 4),
        num(row["interior_max_site_relative_change"], 4),
        yn(row["stationary_gate_pass"]),
    ] for row in stationarity])

    points_table = tex_rows([[
        row["temperature"], row["metric"].replace("_", r"\_"), row["n"], row["j"],
        num(row["mean_I"]), num(row["se_mean_I"]), yn(row["stationary"]),
    ] for row in sorted(scaling_points,
                        key=lambda row: (row["temperature"], row["metric"], int(row["n"])))])

    scaling_table = tex_rows([[
        row["temperature"], row["metric"].replace("_", r"\_"),
        num(row["exponent"]), num(row["exponent_se_ols"]),
        f"[{num(row['ci_low_df1'])}, {num(row['ci_high_df1'])}]",
        num(row["log_sse"]), yn(row["all_stationary"]),
    ] for row in scaling])

    penetration_table = tex_rows([[
        row["temperature"], row["n"], num(row["left_action"]), num(row["threshold"]),
        yn(row["threshold_reached"]), row["first_site_at_or_below"] or "--",
        num(row["interpolated_site"]), num(row["interpolated_fraction_of_n"], 4),
        yn(row["stationary"]),
    ] for row in penetration])

    spread_table = tex_rows([[
        row["temperature"], row["scale"].replace("_", r"\_"),
        row["region"].replace("_", r"\_"), row["grid_points"],
        num(row["median_relative_spread"], 4), num(row["max_relative_spread"], 4),
    ] for row in spread])

    run_table = tex_rows([[
        row["T1"], row["n"], row["seed"], row["threads"],
        row["valid_trajectories"], row["nonfinite_trajectories"],
        row["projection_count"], num(row["elapsed_seconds"], 2),
    ] for row in sorted(runs, key=lambda row: (float(row["T1"]), int(row["n"]), int(row["seed"])))])

    if verdict["all_six_logical_conditions_stationary"]:
        primary = (
            "All six logical conditions pass the prospectively frozen stationarity gate. "
            "The profile, scaling, penetration, and collapse summaries may therefore be "
            "interpreted as stationary within the resolution of this protocol."
        )
    else:
        failed = ", ".join(verdict["nonstationary_conditions"]).replace("_", r"\_")
        primary = (
            "The open-chain stationarity gate does not pass for every condition. "
            f"Failed conditions: {failed}. Scaling and collapse numbers involving these "
            "conditions are reported descriptively and are not stationary-profile claims."
        )

    figures = "\n".join(
        rf"""\begin{{figure}}[p]
\centering
\includegraphics[width=0.88\linewidth]{{../figures/{stem}.pdf}}
\caption{{{caption}}}
\end{{figure}}"""
        for stem, caption in (
            ("open_stationarity_T10", "Stationarity diagnostics at left-bath temperature 10: instantaneous total action and free-end action."),
            ("open_stationarity_T6", "Stationarity diagnostics at left-bath temperature 6: instantaneous total action and free-end action."),
            ("open_quarter_profiles_T10", "Separate quarter-averaged action profiles at left-bath temperature 10."),
            ("open_quarter_profiles_T6", "Separate quarter-averaged action profiles at left-bath temperature 6."),
            ("open_action_profiles_T10", "Raw open-chain action profiles at left-bath temperature 10; bands are pointwise 95 percent intervals."),
            ("open_action_profiles_T6", "Raw open-chain action profiles at left-bath temperature 6; bands are pointwise 95 percent intervals."),
            ("open_action_profiles_rescaled_T10", r"Fixed-exponent collapse test $\sqrt{n}\langle I_j\rangle$ at left-bath temperature 10."),
            ("open_action_profiles_rescaled_T6", r"Fixed-exponent collapse test $\sqrt{n}\langle I_j\rangle$ at left-bath temperature 6."),
        )
    )

    tex = rf"""\documentclass[10pt]{{article}}
\usepackage[margin=0.68in]{{geometry}}
\usepackage{{amsmath,booktabs,longtable,graphicx,pdflscape,hyperref}}
\hypersetup{{colorlinks=true,linkcolor=blue,urlcolor=blue}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.4em}}
\title{{One-bath Open-chain Profile Diagnostic}}
\author{{Frozen numerical audit}}
\date{{September 2026}}
\begin{{document}}
\sloppy
\maketitle

\section{{Boundary definition}}
The left endpoint carries the canonical BC1 bath,
\[
b_{{I,L}}=2\gamma\left[2T_L-(2MI_1-I_1^2+2I_1I_2\cos\delta_1)\right],\qquad
b_{{\phi,L}}=2\gamma I_2\sin\delta_1,
\]
with $\sigma_{{I,L}}=2\sqrt{{2\gamma T_LI_1}}$ and
$\sigma_{{\phi,L}}=\sqrt{{2\gamma T_L/I_1}}$.  At the free endpoint,
$b_{{I,R}}=b_{{\phi,R}}=\sigma_{{I,R}}=\sigma_{{\phi,R}}=0$ exactly; only
Hamiltonian drift remains.  Thus this is the genuinely open definition from
\texttt{{flux\_V2.cpp}}, not the non-open historical \texttt{{main\_fixed.cpp}}
cases.

The port deliberately retains the corrected-noise SIMD profile integrator so
the action profiles can be compared directly with the existing 72-run matrix.
It also retains single precision, adaptive Euler--Maruyama stepping,
polynomial trigonometry, and positive-action projection.  The saved bond-sine
values are therefore archival only: the independent controlled-integrator
audit showed that this scheme produces a spurious sine signal.

\section{{Stationarity result (reported first)}}
{primary}

The gate compares the separate third and fourth measurement quarters.  It
requires total-action and free-end paired differences to contain zero and be
at most 5 percent in relative magnitude; the interior sitewise median and
maximum changes must be at most 5 and 10 percent, respectively.

\begin{{landscape}}
\scriptsize
\begin{{longtable}}{{l r r r r c r r r c r r c}}
\toprule
$T_L$ & $n$ & $\Delta M$ & SE & rel. & CI0 & $\Delta I_n$ & SE & rel. & CI0 & bulk med. & bulk max & pass\\
\midrule
\endfirsthead
{stationarity_table}
\bottomrule
\end{{longtable}}
\end{{landscape}}

{figures}

\clearpage
\section{{Mid-chain and free-end scaling}}
\scriptsize
\begin{{longtable}}{{l l r r r r c}}
\toprule
$T_L$ & metric & $n$ & site & mean & SE & stationary\\
\midrule
\endfirsthead
{points_table}
\bottomrule
\end{{longtable}}

\begin{{longtable}}{{l l r r l r c}}
\toprule
$T_L$ & metric & exponent & OLS SE & 95\% CI ($df=1$) & log-SSE & all stationary\\
\midrule
\endfirsthead
{scaling_table}
\bottomrule
\end{{longtable}}
\normalsize
Each exponent uses only three chain lengths and therefore has one residual
degree of freedom.  No exponent is adjusted to improve collapse.

\section{{Fixed $\sqrt{{n}}$ collapse}}
Profiles are linearly interpolated to the predeclared 101-point grid
$x=0,0.01,\ldots,1$.  Relative spread is
$(\max_n y-\min_n y)/\operatorname{{mean}}_n y$ and is reported on the full
grid and fixed interior $0.10\le x\le0.90$.

\scriptsize
\begin{{longtable}}{{l l l r r r}}
\toprule
$T_L$ & scale & region & grid points & median spread & maximum spread\\
\midrule
\endfirsthead
{spread_table}
\bottomrule
\end{{longtable}}
\normalsize

\section{{Penetration depth}}
The prospectively fixed threshold is 50 percent of the driven-end action.  The
first raw crossing and its linearly interpolated location are reported without
tail extrapolation.

\scriptsize
\begin{{longtable}}{{l r r r c r r r c}}
\toprule
$T_L$ & $n$ & $I_1$ & threshold & reached & first site & interpolated site & site/$n$ & stationary\\
\midrule
\endfirsthead
{penetration_table}
\bottomrule
\end{{longtable}}
\normalsize

\clearpage
\section{{Run-level integrity and provenance}}
\scriptsize
\begin{{longtable}}{{r r r r r r r r}}
\toprule
$T_L$ & $n$ & seed & threads & valid & nonfinite & projections & seconds\\
\midrule
\endfirsthead
{run_table}
\bottomrule
\end{{longtable}}
\normalsize

Parent profile-matrix source SHA-256:
\path{{40d13c7547d9c147ca8241ed701e1a1e85540fc0c21db48c4a63a943eccfe3ec}}.

Remote freeze commit:
\path{{{production['remote_freeze_commit']}}}.

Open-chain source SHA-256:
\path{{{manifest['source_sha256']}}}.

Binary SHA-256:
\path{{{manifest['binary_sha256']}}}.

Analysis SHA-256:
\path{{{manifest['analysis_sha256']}}}.

Exact commands, seeds, runtimes, complete raw files, merged profiles, separate
quarter profiles, snapshot trajectories, spread curves, and all SHA-256 values
accompany this report.

\end{{document}}
"""
    (report / "open_chain_profiles_report.tex").write_text(tex)


if __name__ == "__main__":
    main()
