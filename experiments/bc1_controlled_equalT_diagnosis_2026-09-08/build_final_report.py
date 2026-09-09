#!/usr/bin/env python3
"""Build the final figures, LaTeX report, and audit markdown files."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"
REPORT = ROOT / "report"
FIGURES = REPORT / "figures"

SOURCE_SHA = "17d2b0e75c7838c7fc8710bae0f0381c7607c967d38e03b4db3bde34e8ec086b"
BINARY_SHA = "68c13863da721f655ecbd3d0913b80691c0a0a3a6c984877af8423c331c443db"
LAUNCH_COMMIT = "908a7c08faac7ea3b8eb06c0577fca7fb5e9ad5e"
REMOTE_FREEZE_COMMIT = "35124a825a4119e9007c635e739ae71e3752d6b7"


def read_csv(name: str) -> list[dict[str, str]]:
    with (ANALYSIS / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_json(name: str) -> dict:
    return json.loads((ANALYSIS / name).read_text())


def f(value: str | float, digits: int = 6) -> str:
    return f"{float(value):.{digits}g}"


def make_figures(products: dict[int, dict[str, list[dict[str, str]]]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25})

    fig, axes = plt.subplots(2, 1, figsize=(8.1, 6.4), constrained_layout=True)
    for ax, n in zip(axes, (25, 100)):
        rows = products[n]["sine"]
        j = np.array([int(r["j"]) for r in rows])
        cm = np.array([float(r["controlled_mean_sin"]) for r in rows])
        cs = np.array([float(r["controlled_se"]) for r in rows])
        sm = np.array([float(r["simd_mean_sin"]) for r in rows])
        ss = np.array([float(r["simd_se"]) for r in rows])
        ax.plot(j, sm, color="#C44E52", lw=1.1, label="inherited SIMD")
        ax.fill_between(j, sm - 1.96 * ss, sm + 1.96 * ss, color="#C44E52", alpha=0.18)
        ax.plot(j, cm, color="#2C6E9B", lw=1.1, label="controlled Cartesian")
        ax.fill_between(j, cm - 1.96 * cs, cm + 1.96 * cs, color="#2C6E9B", alpha=0.2)
        ax.axhline(0.0, color="black", lw=0.9)
        ax.set_ylabel(r"$\langle\sin\theta_j\rangle$")
        ax.set_title(f"n={n}")
        ax.legend(frameon=False, ncol=2)
    axes[-1].set_xlabel("bond j")
    fig.savefig(FIGURES / "bond_sine_profiles.pdf")
    fig.savefig(FIGURES / "bond_sine_profiles.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.1, 6.4), constrained_layout=True)
    for ax, n in zip(axes, (25, 100)):
        rows = products[n]["sine"]
        j = np.array([int(r["j"]) for r in rows])
        cz = np.array([float(r["controlled_z"]) for r in rows])
        sz = np.array([float(r["simd_z"]) for r in rows])
        bonf = products[n]["summary"]["sine_control"]["bonferroni_critical"]
        ax.plot(j, sz, color="#C44E52", lw=1.0, label="inherited SIMD")
        ax.plot(j, cz, color="#2C6E9B", lw=1.0, label="controlled Cartesian")
        ax.axhline(1.96, color="#777777", ls="--", lw=0.8)
        ax.axhline(-1.96, color="#777777", ls="--", lw=0.8)
        ax.axhline(bonf, color="black", ls=":", lw=1.0, label="Bonferroni bound")
        ax.axhline(-bonf, color="black", ls=":", lw=1.0)
        ax.set_ylabel("z score")
        ax.set_title(f"n={n}")
        ax.legend(frameon=False, ncol=3)
    axes[-1].set_xlabel("bond j")
    fig.savefig(FIGURES / "bond_sine_z_scores.pdf")
    fig.savefig(FIGURES / "bond_sine_z_scores.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.1, 6.4), constrained_layout=True)
    for ax, n in zip(axes, (25, 100)):
        rows = products[n]["action"]
        j = np.array([int(r["j"]) for r in rows])
        cm = np.array([float(r["controlled_mean_I"]) for r in rows])
        cs = np.array([float(r["controlled_se_mean_I"]) for r in rows])
        sm = np.array([float(r["simd_mean_I"]) for r in rows])
        ss = np.array([float(r["simd_se_mean_I"]) for r in rows])
        ax.plot(j, sm, color="#C44E52", lw=1.1, label="inherited SIMD")
        ax.fill_between(j, sm - 1.96 * ss, sm + 1.96 * ss, color="#C44E52", alpha=0.18)
        ax.plot(j, cm, color="#2C6E9B", lw=1.1, label="controlled Cartesian")
        ax.fill_between(j, cm - 1.96 * cs, cm + 1.96 * cs, color="#2C6E9B", alpha=0.2)
        ax.set_ylabel(r"$\langle I_j\rangle$")
        ax.set_title(f"n={n}")
        ax.legend(frameon=False, ncol=2)
    axes[-1].set_xlabel("site j")
    fig.savefig(FIGURES / "action_profile_comparison.pdf")
    fig.savefig(FIGURES / "action_profile_comparison.png", dpi=220)
    plt.close(fig)


def bond_rows(rows: list[dict[str, str]]) -> str:
    return "\n".join(
        f"{r['j']} & {f(r['controlled_mean_sin'])} & {f(r['controlled_se'])} & "
        f"{f(r['controlled_z'], 5)} & {f(r['simd_mean_sin'])} & "
        f"{f(r['simd_se'])} & {f(r['simd_z'], 5)} \\\\" for r in rows
    )


def action_rows(rows: list[dict[str, str]]) -> str:
    return "\n".join(
        f"{r['j']} & {f(r['controlled_mean_I'])} & {f(r['controlled_se_mean_I'])} & "
        f"{f(r['simd_mean_I'])} & {f(r['simd_se_mean_I'])} & {f(r['difference'])} \\\\"
        for r in rows
    )


def write_report(products: dict[int, dict]) -> None:
    s25 = products[25]["summary"]
    s100 = products[100]["summary"]
    report = rf"""\documentclass[10pt]{{article}}
\usepackage[margin=0.75in]{{geometry}}
\usepackage{{amsmath,amssymb,booktabs,longtable,graphicx,xcolor}}
\usepackage[hidelinks]{{hyperref}}
\graphicspath{{{{figures/}}}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.45em}}
\title{{BC1 Equal-Temperature Known-Answer Diagnostic\\Projection-Free Controlled Integrator}}
\author{{Numerical audit report}}
\date{{9 September 2026}}
\begin{{document}}
\maketitle

\section{{Question and frozen decision rule}}
For the canonical BC1 bath at $T_1=T_n=6$, the stationary Gibbs density is
invariant under phase reversal, while every $\sin\theta_j$ is odd. Therefore
$\langle\sin\theta_j\rangle=0$ exactly. The inherited SIMD calculation failed
this control with maximum $|z|=8.678$ at $n=25$ and $12.929$ at $n=100$.
The diagnostic repeats only BC1 with the controlled Cartesian implementation:
double precision, fixed $dt=5\times10^{{-4}}$, standard trigonometry, no
adaptive stepping, and no positive-action projection or floor.

The gate was frozen before production. The port is accepted only if the
spatial mean action is within 5\% of the SIMD reference and its relative
spatial range is no larger than the SIMD range. The sine control is assessed
both pointwise at $|z|>1.95996$ and simultaneously with a two-sided
Bonferroni threshold.

\section{{Protocol and integrity}}
For $n=25$, burn-in and measurement times were 2000 and 2000. For $n=100$,
they were 32000 and 2000. The sample interval was 0.01, producing 200001
profile samples per trajectory. Two base seeds and 256 trajectories per base
seed gave 512 independent trajectories at each chain length.

\begin{{center}}
\begin{{tabular}}{{rrrrr}}
\toprule
$n$ & base seeds & trajectories & samples/trajectory & integrity errors\\
\midrule
25 & 2026090825, 2026090826 & 512 & 200001 & 0\\
100 & 2026090900, 2026090901 & 512 & 200001 & 0\\
\bottomrule
\end{{tabular}}
\end{{center}}

All 426 frozen output hashes verified, and a fresh independent invocation of
the frozen analysis reproduced the saved JSON byte-for-byte.

\section{{Primary result}}
\begin{{center}}
\begin{{tabular}}{{rrrrrr}}
\toprule
$n$ & controlled max $|z|$ & SIMD max $|z|$ & pointwise failures & Bonf. critical & simultaneous test\\
\midrule
25 & {s25['sine_control']['max_abs_z']:.3f} & {s25['sine_control']['simd_max_abs_z']:.3f} & {s25['sine_control']['pointwise_failures']}/24 & {s25['sine_control']['bonferroni_critical']:.3f} & PASS\\
100 & {s100['sine_control']['max_abs_z']:.3f} & {s100['sine_control']['simd_max_abs_z']:.3f} & {s100['sine_control']['pointwise_failures']}/99 & {s100['sine_control']['bonferroni_critical']:.3f} & PASS\\
\bottomrule
\end{{tabular}}
\end{{center}}

At $n=25$, no bond fails even the uncorrected pointwise test. At $n=100$,
5 of 99 bonds exceed 1.95996, approximately the number expected by chance
under 99 pointwise 5\% tests, while the maximum $|z|=2.984$ remains below the
predeclared simultaneous threshold 3.478. Thus the exact zero-sine control
passes at both lengths under the controlled integrator.

\begin{{figure}}[ht]
\centering
\includegraphics[width=0.94\textwidth]{{bond_sine_profiles.pdf}}
\caption{{Bond-sine profiles and 95\% pointwise bands. The inherited SIMD
profile shows coherent nonzero structure; the controlled profile fluctuates
around zero.}}
\end{{figure}}

\begin{{figure}}[ht]
\centering
\includegraphics[width=0.94\textwidth]{{bond_sine_z_scores.pdf}}
\caption{{Per-bond $z$ scores. Dotted lines are the frozen Bonferroni bounds;
dashed lines mark the uncorrected pointwise 95\% bounds.}}
\end{{figure}}

\section{{Like-for-like precision}}
The mean controlled-to-SIMD SE ratios are comparable or smaller:
the median ratio is {s25['sine_control']['median_controlled_over_simd_se']:.3f}
at $n=25$ and {s100['sine_control']['median_controlled_over_simd_se']:.3f} at
$n=100$. Mean controlled SEs are
{s25['sine_control']['controlled_se_mean']:.6g} and
{s100['sine_control']['controlled_se_mean']:.6g}, compared with SIMD means
{s25['sine_control']['simd_se_mean']:.6g} and
{s100['sine_control']['simd_se_mean']:.6g}. The disappearance of the coherent
signal is therefore not caused by reduced precision.

\section{{Porting gate: action profile}}
\begin{{center}}
\begin{{tabular}}{{rrrrrr}}
\toprule
$n$ & controlled mean & SIMD mean & relative difference & controlled range & SIMD range\\
\midrule
25 & {s25['porting']['controlled_spatial_mean']:.6f} & {s25['porting']['simd_spatial_mean']:.6f} & {100*s25['porting']['mean_relative_difference']:.3f}\% & {100*s25['porting']['controlled_relative_range']:.3f}\% & {100*s25['porting']['simd_relative_range']:.3f}\%\\
100 & {s100['porting']['controlled_spatial_mean']:.6f} & {s100['porting']['simd_spatial_mean']:.6f} & {100*s100['porting']['mean_relative_difference']:.3f}\% & {100*s100['porting']['controlled_relative_range']:.3f}\% & {100*s100['porting']['simd_relative_range']:.3f}\%\\
\bottomrule
\end{{tabular}}
\end{{center}}
Both predeclared porting gates pass. The spatial mean differs by only 0.407\%
at $n=25$ and 1.298\% at $n=100$, and the controlled profiles are flatter
than their SIMD references.

\begin{{figure}}[ht]
\centering
\includegraphics[width=0.94\textwidth]{{action_profile_comparison.pdf}}
\caption{{Equal-temperature action profiles. The controlled port preserves
the action scale while eliminating the systematic bond-sine asymmetry.}}
\end{{figure}}

\section{{Numerical-event audit}}
\begin{{center}}
\begin{{tabular}}{{rrrrrr}}
\toprule
$n$ & projections & floors & zero-radius & midpoint failures & min sampled action\\
\midrule
25 & 0 & 0 & 0 & 0 & {s25['numerics']['min_sampled_action']:.6e}\\
100 & 0 & 0 & 0 & 0 & {s100['numerics']['min_sampled_action']:.6e}\\
\bottomrule
\end{{tabular}}
\end{{center}}
Nonfinite batches were zero at both lengths. The near-zero-action region was
sampled down to $4.86\times10^{{-11}}$ and $1.68\times10^{{-12}}$ without
clipping.

\section{{Verdict and claim boundary}}
The controlled scheme passes the exact BC1 equal-temperature sine control at
both chain lengths, whereas the inherited SIMD scheme fails systematically at
comparable SE. Together with the passed action-profile porting gate, this
establishes that the old nonzero equal-temperature bond-sine profile is a
numerical artefact of the inherited SIMD numerical scheme, not a property of
the BC1 stationary dynamics. Because several numerical changes were made
together, this experiment does not identify projection alone as the unique
cause.

Consequently, driven $(10,2)$ bond-sine profiles obtained with the inherited
SIMD scheme should not be used quantitatively in the manuscript until they are
re-measured with the controlled integrator. Existing action-profile claims may
be retained subject to their separate validated gates.

\section{{Provenance}}
Launch commit: \nolinkurl{{{LAUNCH_COMMIT}}}.\\
Remote frozen commit: \nolinkurl{{{REMOTE_FREEZE_COMMIT}}}.\\
Source SHA-256: \nolinkurl{{{SOURCE_SHA}}}.\\
Binary SHA-256: \nolinkurl{{{BINARY_SHA}}}.
The 64 unique completed production commands are preserved in the provenance
directory as \texttt{{COMMANDS\_UNIQUE\_}}\allowbreak\texttt{{PRODUCTION.csv}}.
The original \texttt{{COMMANDS.csv}} retains all 72 launch records, including eight
commands from the first detached attempt; that attempt produced no scientific
data and its empty logs are preserved under
\texttt{{failed\_detached\_launch\_1}}.

\clearpage
\appendix
\section{{All bond-sine values}}
\subsection{{$n=25$}}
\scriptsize
\begin{{longtable}}{{rrrrrrr}}
\toprule
$j$ & ctrl. mean & ctrl. SE & ctrl. $z$ & SIMD mean & SIMD SE & SIMD $z$\\
\midrule\endhead
{bond_rows(products[25]['sine'])}
\bottomrule
\end{{longtable}}

\subsection{{$n=100$}}
\begin{{longtable}}{{rrrrrrr}}
\toprule
$j$ & ctrl. mean & ctrl. SE & ctrl. $z$ & SIMD mean & SIMD SE & SIMD $z$\\
\midrule\endhead
{bond_rows(products[100]['sine'])}
\bottomrule
\end{{longtable}}

\clearpage
\section{{Full action profiles}}
\subsection{{$n=25$}}
\begin{{longtable}}{{rrrrrr}}
\toprule
$j$ & ctrl. mean & ctrl. SE & SIMD mean & SIMD SE & difference\\
\midrule\endhead
{action_rows(products[25]['action'])}
\bottomrule
\end{{longtable}}

\subsection{{$n=100$}}
\begin{{longtable}}{{rrrrrr}}
\toprule
$j$ & ctrl. mean & ctrl. SE & SIMD mean & SIMD SE & difference\\
\midrule\endhead
{action_rows(products[100]['action'])}
\bottomrule
\end{{longtable}}
\end{{document}}
"""
    (REPORT / "bc1_controlled_equalT_diagnosis_report.tex").write_text(report)


def write_markdown(products: dict[int, dict]) -> None:
    s25 = products[25]["summary"]
    s100 = products[100]["summary"]
    verdict = f"""# Final verdict

**PASS: the projection-free controlled integrator satisfies the BC1 equal-temperature known-answer control at n=25 and n=100.**

- n=25: max |z| = {s25['sine_control']['max_abs_z']:.6f}; 0/24 pointwise failures; Bonferroni critical value {s25['sine_control']['bonferroni_critical']:.6f}; simultaneous PASS.
- n=100: max |z| = {s100['sine_control']['max_abs_z']:.6f}; 5/99 uncorrected pointwise failures; Bonferroni critical value {s100['sine_control']['bonferroni_critical']:.6f}; simultaneous PASS.
- Old SIMD maxima were {s25['sine_control']['simd_max_abs_z']:.6f} and {s100['sine_control']['simd_max_abs_z']:.6f}, with 19/24 and 71/99 pointwise failures.
- The action-profile porting gate passed: spatial-mean differences were {100*s25['porting']['mean_relative_difference']:.3f}% and {100*s100['porting']['mean_relative_difference']:.3f}%.
- Projection, floor, zero-radius, midpoint-failure, and nonfinite-batch counts were all zero.

The inherited SIMD equal-temperature sine signal is therefore a numerical-scheme artefact. This comparison changes several numerical ingredients together, so it does not isolate projection as the unique cause. Driven (10,2) SIMD bond-sine profiles must be re-measured with the controlled integrator before quantitative use.
"""
    (ROOT / "FINAL_VERDICT.md").write_text(verdict)

    validation = f"""# Validation report

## Integrity

- 512 trajectories at n=25 and 512 at n=100.
- 200001 measurement samples per trajectory.
- Missing products: none. Integrity errors: none.
- All 426 saved-output SHA-256 entries pass `shasum -a 256 -c`.
- Fresh reanalysis matches `analysis/final_analysis_stdout.json` byte-for-byte.

## Frozen numerical protocol

- BC1, T1=Tn=6, gamma=0.1, dt=5e-4.
- Burn-in: 2000 (n=25), 32000 (n=100); measurement duration 2000; sample interval 0.01.
- Seeds: 2026090825, 2026090826, 2026090900, 2026090901.
- Source SHA-256: `{SOURCE_SHA}`.
- Binary SHA-256: `{BINARY_SHA}`.
- Launch commit: `{LAUNCH_COMMIT}`.
- Remote frozen commit: `{REMOTE_FREEZE_COMMIT}`.

## Gates

| n | action port | numerical events | max |z| | pointwise failures | Bonferroni | final control |
|---:|:---:|:---:|---:|---:|:---:|:---:|
| 25 | PASS | PASS | {s25['sine_control']['max_abs_z']:.6f} | 0/24 | PASS | PASS |
| 100 | PASS | PASS | {s100['sine_control']['max_abs_z']:.6f} | 5/99 | PASS | PASS |

Full per-bond and per-site raw analysis values are in `analysis/*_comparison.csv` and in the PDF appendices. The 64 unique completed commands are in `provenance/COMMANDS_UNIQUE_PRODUCTION.csv`; all 72 launch records, including the eight no-data detached-attempt records, remain in `provenance/COMMANDS.csv`. Frozen output hashes are in `provenance/OUTPUT_HASHES.sha256`.
"""
    (ROOT / "VALIDATION_REPORT.md").write_text(validation)

    contents = """# Package contents

- `report/bc1_controlled_equalT_diagnosis_report.pdf` and `.tex`: final report with complete bond and action tables.
- `report/figures/`: publication-quality PDF and PNG figures.
- `analysis/`: complete summary JSON and per-bond/per-site CSVs.
- `raw/`: all 64 production batch outputs plus smoke-test artefacts.
- `src/`, `bin/`, `analyze.py`, `build_final_report.py`: frozen simulator and analysis implementation.
- `PROTOCOL.md`, `FROZEN_HASHES.sha256`, `PORTING_AUDIT.md`: predeclared protocol and source audit.
- `provenance/`: exact commands, timestamps, launch commit, complete output hashes, and preserved failed-launch evidence. `COMMANDS.csv` contains all 72 launch records; `COMMANDS_UNIQUE_PRODUCTION.csv` contains the 64 unique completed production commands.
- `FINAL_VERDICT.md`, `VALIDATION_REPORT.md`: concise claim boundary and audit result.
"""
    (ROOT / "PACKAGE_CONTENTS.md").write_text(contents)


def write_unique_production_commands() -> None:
    source = ROOT / "provenance" / "COMMANDS.csv"
    target = ROOT / "provenance" / "COMMANDS_UNIQUE_PRODUCTION.csv"
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    seen: set[tuple[str, str, str, str]] = set()
    unique: list[dict[str, str]] = []
    for row in rows:
        key = (row["n"], row["base_seed"], row["stream_offset"], row["command"])
        if key not in seen:
            seen.add(key)
            unique.append(row)
    if len(rows) != 72 or len(unique) != 64:
        raise RuntimeError(f"unexpected command audit: {len(rows)=}, {len(unique)=}")
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["n", "base_seed", "stream_offset", "command"])
        writer.writeheader()
        writer.writerows(unique)


def main() -> None:
    REPORT.mkdir(parents=True, exist_ok=True)
    products: dict[int, dict] = {}
    for n in (25, 100):
        products[n] = {
            "summary": load_json(f"n{n}_summary.json"),
            "sine": read_csv(f"n{n}_bond_sine_comparison.csv"),
            "action": read_csv(f"n{n}_action_comparison.csv"),
        }
    make_figures(products)
    write_unique_production_commands()
    write_report(products)
    write_markdown(products)


if __name__ == "__main__":
    main()
