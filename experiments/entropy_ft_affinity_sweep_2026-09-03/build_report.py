#!/usr/bin/env python3
"""Build a concise LaTeX report from the frozen affinity-sweep CSV outputs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: str | float, digits: int = 4) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "--"
    if not math.isfinite(numeric):
        return "--"
    if numeric == 0.0:
        return "0"
    if abs(numeric) < 1.0e-3 or abs(numeric) >= 1.0e4:
        return f"{numeric:.{digits}e}"
    return f"{numeric:.{digits}f}"


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "_": r"\_",
        "%": r"\%",
        "&": r"\&",
        "#": r"\#",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(char, char) for char in text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.experiment_dir
    analysis = root / "analysis"
    window = rows(analysis / "window_summary.csv")
    extrapolation = rows(analysis / "infinite_time_extrapolation.csv")
    crossover = rows(analysis / "crossover_summary.csv")
    hashes = rows(analysis / "input_hashes.csv")
    manifest_text = (root / "production_manifest.txt").read_text()
    manifest = {}
    for line in manifest_text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            manifest[key] = value

    equilibrium_manifest = {}
    equilibrium_manifest_path = root / "EQUILIBRIUM_PRODUCTION_MANIFEST.txt"
    if equilibrium_manifest_path.is_file():
        for line in equilibrium_manifest_path.read_text().splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                equilibrium_manifest[key] = value

    commands = (root / "COMMANDS.tsv").read_text().splitlines()[1:]
    equilibrium_commands_path = root / "EQUILIBRIUM_COMMAND.tsv"
    if equilibrium_commands_path.is_file():
        commands.extend(equilibrium_commands_path.read_text().splitlines()[1:])
    args.output.parent.mkdir(parents=True, exist_ok=True)

    parts: list[str] = []
    parts.append(r"""\documentclass[10pt]{article}
\usepackage[margin=0.72in]{geometry}
\usepackage{booktabs,longtable,graphicx,amsmath,xcolor,hyperref,fvextra}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}
\newcommand{\NA}{\textemdash}
\begin{document}
\begin{center}
{\Large Heat-flux affinity sweep and Gaussian crossover}\par
\vspace{0.3em}
{\small Frozen $n=10$ Cartesian experiment with auditable bootstrap-gate amendment, 2026-09-04}
\end{center}

\section*{Scope and frozen estimator}
The transported heat in a window is
\[
Q=\frac{Q_L-Q_R}{2},\qquad
\Delta\beta=\frac{1}{T_R}-\frac{1}{T_L}.
\]
Non-overlapping windows are formed within each of 128 independent streams.
The directly sampled fluctuation-relation slope is obtained from raw matched
bins of $\log[p(Q)/p(-Q)]$ against $Q$.  The fixed rule is
$\Delta Q=\operatorname{std}(Q)/20$, at least 10 raw observations on each
side of every accepted bin pair, and the largest contiguous reliable block
that straddles zero.  No plus-four points or extrapolated densities enter the
fit.  Stream bootstrap intervals resample the 128 complete streams.

The Gaussian bulk diagnostics are
\[
a_{\rm Gauss}=\frac{2\langle Q\rangle}{\operatorname{Var}(Q)},\qquad
G_{\rm FT}=\frac{\operatorname{Var}(Q)}{2\langle Q\rangle/\Delta\beta}.
\]
They characterize the CLT crossover but are not substitutes for a directly
sampled negative tail.
""")

    parts.append(r"\section*{Provenance}")
    parts.append(
        "Experiment repository commit: \\nolinkurl{%s}.\\\\\n"
        "Production source commit: \\nolinkurl{%s}.\\\\\n"
        "Frozen source SHA-256: \\nolinkurl{%s}.\\\\\n"
        "Driven binary SHA-256: \\nolinkurl{%s}.\\\\\n"
        "Equilibrium binary SHA-256: \\nolinkurl{%s}.\\\\\n"
        % (
            latex_escape(manifest.get("experiment_repository_commit", "missing")),
            latex_escape(manifest.get("production_source_commit", "missing")),
            latex_escape(manifest.get("production_source_sha256", "missing")),
            latex_escape(manifest.get("binary_sha256", "missing")),
            latex_escape(equilibrium_manifest.get("binary_sha256", "missing")),
        )
    )
    parts.append(
        r"The four driven runs and the equilibrium amendment use $n=10$, "
        r"$\gamma=0.1$, $dt=5\times10^{-4}$, "
        r"burn-in 500, base block duration 20, 8 SIMD batches, 16 lanes, "
        r"7,813 blocks per stream, and 1,000,064 blocks per case."
    )

    parts.append(r"\paragraph{Exact production commands.}")
    parts.append(
        r"\begin{Verbatim}[fontsize=\scriptsize,breaklines=true,"
        r"breakanywhere=true]"
    )
    for command in commands:
        parts.append(command.replace("\t", "  "))
    parts.append(r"\end{Verbatim}")

    parts.append(r"\paragraph{Input hashes.}")
    parts.append(r"\begin{longtable}{lrp{0.62\textwidth}}")
    parts.append(r"\toprule Case & rows & SHA-256 \\")
    parts.append(r"\midrule\endhead")
    for row in hashes:
        parts.append(
            r"\texttt{%s} & %s & \scriptsize\nolinkurl{%s} \\"
            % (
                latex_escape(row["case"]),
                row["rows"],
                latex_escape(row["sha256"]),
            )
        )
    parts.append(r"\bottomrule\end{longtable}")
    parts.append(
        r"\noindent Full absolute input paths are retained in "
        r"\texttt{analysis/input\_hashes.csv}."
    )

    parts.append(r"\section*{Bootstrap-gate amendment and audit trail}")
    parts.append(
        r"The original all-resolved-time gate reported 0/1000 accepted "
        r"intercepts for every case because sub-threshold raw acceptance was "
        r"written as zero.  The amended analysis leaves the production data, "
        r"$\Delta Q$, minimum count, minimum three pairs, and straddle-zero "
        r"rule unchanged.  Each $t$ is bootstrapped independently on its "
        r"frozen full-sample bin grid.  Driven-case intercepts use only "
        r"full-sample fits with $n_-\ge500$ and $R^2\ge0.98$."
    )
    parts.append(r"\begin{longtable}{lrrrlr}")
    parts.append(
        r"\toprule Case & old reported & old raw valid & old $t$ set & "
        r"amended $t$ set & amended valid \\"
    )
    parts.append(r"\midrule\endhead")
    for row in extrapolation:
        amended_valid = (
            f'{row["reliable_joint_bootstrap_resolved"]}/'
            f'{row["reliable_joint_bootstrap_total"]}'
            if int(row["n_reliable_times_used"]) >= 3
            else "no intercept"
        )
        parts.append(
            r"$\Delta\beta=%s$ & %s/%s & %s/%s & %s & %s & %s \\"
            % (
                number(row["delta_beta"], 6),
                row["old_joint_bootstrap_reported"],
                row["old_joint_bootstrap_total"],
                row["old_joint_bootstrap_resolved"],
                row["old_joint_bootstrap_total"],
                latex_escape(row["old_resolved_times"]),
                latex_escape(row["reliable_times_used"]),
                amended_valid,
            )
        )
    parts.append(r"\bottomrule\end{longtable}")
    parts.append(
        r"At equilibrium the literal $R^2\ge0.98$ set is empty because the "
        r"exact reference is flat.  The displayed equilibrium set therefore "
        r"uses $n_-\ge500$ without the inapplicable nonzero-signal $R^2$ gate."
    )

    parts.append(r"\section*{Per-time stream-bootstrap intervals}")
    parts.append(
        r"The final column tests the finite-time interval against "
        r"$\Delta\beta$ (zero at equilibrium), not the long-time verdict."
    )
    parts.append(r"\scriptsize\begin{longtable}{lrrrrrrrl}")
    parts.append(
        r"\toprule $\Delta\beta$ & $t$ & $n_-$ & $R^2$ & $a_{\rm fit}$ & "
        r"CI low & CI high & accepted & in CI \\"
    )
    parts.append(r"\midrule\endhead")
    for row in sorted(
        (item for item in window if item["resolved"] == "1"),
        key=lambda item: (float(item["delta_beta"]), int(item["t"])),
    ):
        parts.append(
            "%s & %s & %s & %s & %s & %s & %s & %s/%s & %s \\\\"
            % (
                number(row["delta_beta"], 6),
                row["t"],
                row["n_negative"],
                number(row["a_fit_R2"], 4),
                number(row["a_fit"], 7),
                number(row["per_t_bootstrap_ci_low"], 7),
                number(row["per_t_bootstrap_ci_high"], 7),
                row["per_t_bootstrap_resolved"],
                row["per_t_bootstrap_total"],
                "yes" if row["per_t_ci_contains_reference"] == "1" else "no",
            )
        )
    parts.append(r"\bottomrule\end{longtable}\normalsize")
    parts.append(
        r"Finite-time corrections remain visible.  Rows failing the fixed "
        r"negative-count or $R^2$ gate are reported here but excluded from "
        r"the amended intercept."
    )

    parts.append(r"\section*{Original per-window outputs retained for audit}")
    parts.append(
        r"The original bootstrap CI column below uses the fixed-full-window "
        r"conditional bootstrap.  All full-sample moments and slopes are unchanged."
    )
    by_case: dict[str, list[dict[str, str]]] = {}
    for row in window:
        by_case.setdefault(row["case"], []).append(row)
    for case_name, case_window in by_case.items():
        beta = float(case_window[0]["delta_beta"])
        tl = case_window[0]["T_left"]
        tr = case_window[0]["T_right"]
        parts.append(
            r"\subsection*{%s: $(T_L,T_R)=(%s,%s)$, $\Delta\beta=%s$}"
            % (latex_escape(case_name), tl, tr, number(beta, 6))
        )
        parts.append(r"\scriptsize\begin{longtable}{rrrrrrrrrrr}")
        parts.append(
            r"\toprule $t$ & $N$ & mean & std & skew & ex.kurt & $n_-$ & "
            r"$a_{fit}$ & boot.95\% CI & $R^2$ & $a_G$ / $G_{FT}$ \\"
        )
        parts.append(r"\midrule\endhead")
        for row in sorted(case_window, key=lambda item: int(item["t"])):
            if row["resolved"] == "1":
                slope = number(row["a_fit"], 5)
                interval = "[%s,%s]" % (
                    number(row["a_bootstrap_ci_low"], 5),
                    number(row["a_bootstrap_ci_high"], 5),
                )
                r2 = number(row["a_fit_R2"], 4)
            else:
                slope = r"\textsc{unres.}"
                interval = "--"
                r2 = "--"
            parts.append(
                "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s / %s \\\\"
                % (
                    row["t"], row["N_windows"], number(row["mean_Q"]),
                    number(row["std_Q"]), number(row["skew_Q"]),
                    number(row["excess_kurtosis_Q"]), row["n_negative"],
                    slope, interval, r2, number(row["a_Gauss"], 5),
                    number(row["gaussFT"], 4),
                )
            )
        parts.append(r"\bottomrule\end{longtable}\normalsize")

    parts.append(r"\section*{Amended reliable-time extrapolation}")
    parts.append(r"\begin{center}\resizebox{\textwidth}{!}{%")
    parts.append(r"\begin{tabular}{lrrrrrrl}")
    parts.append(
        r"\toprule Case & $\Delta\beta$ & reliable $t$ & $a_\infty$ & "
        r"95\% CI & $a_\infty/\Delta\beta$ & ratio CI & status \\"
    )
    parts.append(r"\midrule")
    for row in extrapolation:
        parts.append(
            r"\texttt{%s} & %s & %s & %s & [%s,%s] & %s & [%s,%s] & %s \\"
            % (
                latex_escape(row["case"]), number(row["delta_beta"], 6),
                latex_escape(row["reliable_times_used"]), number(row["a_inf"], 5),
                number(row["a_inf_ci_low"], 5), number(row["a_inf_ci_high"], 5),
                number(row["a_inf_over_delta_beta"], 4),
                number(row["ratio_ci_low"], 4), number(row["ratio_ci_high"], 4),
                latex_escape(row["FT_status"]),
            )
        )
    parts.append(r"\bottomrule\end{tabular}}\end{center}")
    parts.append(
        r"The amended intercept accepts 1000/1000 replicates for the first "
        r"three driven rows and equilibrium.  Regression $R^2$ is retained "
        r"as a diagnostic; the weakest-drive linear $1/t$ extrapolation is noisy."
    )

    parts.append(r"\section*{Crossover summary}")
    parts.append(r"\begin{center}\resizebox{\textwidth}{!}{%")
    parts.append(r"\begin{tabular}{rrrrrrrrl}")
    parts.append(
        r"\toprule $\Delta\beta$ & $a_\infty/\Delta\beta$ & ratio CI & "
        r"$a_G/\Delta\beta$ & $G_{FT}(640)$ & skew(160) & ex.kurt(160) & "
        r"$n_-(160)$ & status \\"
    )
    parts.append(r"\midrule")
    for row in crossover:
        parts.append(
            "%s & %s & [%s,%s] & %s & %s & %s & %s & %s & %s \\\\"
            % (
                number(row["delta_beta"], 6),
                number(row["a_inf_over_delta_beta"], 4),
                number(row["ratio_ci_low"], 4), number(row["ratio_ci_high"], 4),
                number(row["a_Gauss_over_delta_beta"], 4),
                number(row["gaussFT_t640"], 4), number(row["skew_t160"], 4),
                number(row["excess_kurtosis_t160"], 4),
                row["n_negative_t160"] or "--", latex_escape(row["FT_status"]),
            )
        )
    parts.append(r"\bottomrule\end{tabular}}\end{center}")

    parts.append(r"""\section*{Figures}
\begin{figure}[ht]
\centering
\includegraphics[width=0.76\textwidth]{../figures/a_fit_vs_inverse_time.pdf}
\caption{Direct matched-bin slopes against inverse window duration.  The
amended intercept uses only the reliable-time subsets stated above.}
\end{figure}
\begin{figure}[ht]
\centering
\includegraphics[width=0.76\textwidth]{../figures/slope_crossover_vs_affinity.pdf}
\caption{Amended reliable-time direct-tail intercepts and Gaussian bulk slopes
against affinity.}
\end{figure}
\begin{figure}[ht]
\centering
\includegraphics[width=0.76\textwidth]{../figures/gaussFT_vs_affinity.pdf}
\caption{Gaussian FT ratio at $t=640$ across affinity.}
\end{figure}

\section*{Equilibrium row and claim boundary}
The added $(6,6)$ production uses the same $n=10$ dynamics, integration step,
burn-in, stream count, and block count as the driven cases.  At
$\Delta\beta=0$, ratios that divide by $\Delta\beta$ and $G_{\rm FT}$ are
undefined and are reported as dashes.  The direct asymmetry-slope reference is
zero; its extrapolated confidence interval is tested against zero rather than
against a ratio of one.

The complete raw positive/negative counts used by every fit are in
\texttt{analysis/symmetric\_bin\_raw\_counts.csv}; all per-window statistics
are in \texttt{analysis/window\_summary.csv}.  Old and amended gates are side
by side in \nolinkurl{analysis/bootstrap_gate_comparison.csv}; pre-amendment
artifacts are under \nolinkurl{pre_gate_respec_2026-09-04/}.  A case labelled
\texttt{CONSISTENT\_WITH\_FT} means only that the amended reliable-time
extrapolated interval includes the reference.  It is not a proof.  A case
whose interval excludes the reference is labelled \texttt{FAIL};
insufficient raw two-tail support remains \texttt{UNRESOLVED}.  The label
\texttt{CONSISTENT\_WITH\_EQUILIBRIUM} has the corresponding limited meaning
for the zero-affinity control.
\end{document}
""")

    args.output.write_text("\n".join(parts))


if __name__ == "__main__":
    main()
