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
{\small Frozen $n=10$ Cartesian experiment with equilibrium amendment, 2026-09-04}
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

    parts.append(r"\section*{Outcome under the frozen gates}")
    parts.append(
        r"The full-sample $1/t$ intercepts give "
        r"$a_\infty/\Delta\beta=1.0175,\ 1.0153,\ 1.0220$ at "
        r"$\Delta\beta=0.027972,\ 0.057143,\ 0.125$, respectively.  "
        r"These point estimates lie close to the FT reference, but none of "
        r"the cases reaches the predeclared threshold of 800 valid joint "
        r"stream-bootstrap intercepts.  Consequently every nonzero-affinity "
        r"FT verdict is \texttt{UNRESOLVED}; the point estimates alone are "
        r"not promoted to a pass."
    )
    parts.append(
        r"For the equilibrium control, the full-sample intercept is "
        r"$7.80\times10^{-5}$.  Its $t=20$ direct slope is "
        r"$1.63\times10^{-4}$ with stream-bootstrap 95\% interval "
        r"$[-8.98\times10^{-5},\,4.18\times10^{-4}]$, which contains zero.  "
        r"The joint long-time interval again fails the 800-replicate gate, "
        r"so the formal equilibrium extrapolation remains "
        r"\texttt{UNRESOLVED}."
    )
    parts.append(
        r"The Gaussian bulk approaches the FT value only at weak drive: "
        r"$a_{\rm Gauss}/\Delta\beta=0.960,\ 0.931,\ 0.745,\ 0.515,\ 0.301$ "
        r"as the affinity increases.  At $t=160$, the negative counts fall "
        r"from 20,136 and 2,997 in the two weakest driven cases to 2, 0, and "
        r"0 in the three stronger cases.  This directly records the loss of "
        r"two-tail support rather than replacing it by a Gaussian estimate."
    )

    parts.append(r"\section*{Raw window statistics and two-tail fits}")
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

    parts.append(r"\section*{Long-time extrapolation}")
    parts.append(r"\begin{longtable}{lrrrrrrl}")
    parts.append(
        r"\toprule Case & $\Delta\beta$ & resolved $t$ & $a_\infty$ & "
        r"95\% CI & $a_\infty/\Delta\beta$ & ratio CI & status \\"
    )
    parts.append(r"\midrule\endhead")
    for row in extrapolation:
        parts.append(
            r"\texttt{%s} & %s & %s & %s & [%s,%s] & %s & [%s,%s] & %s \\"
            % (
                latex_escape(row["case"]), number(row["delta_beta"], 6),
                latex_escape(row["resolved_times"]), number(row["a_inf"], 5),
                number(row["a_inf_ci_low"], 5), number(row["a_inf_ci_high"], 5),
                number(row["a_inf_over_delta_beta"], 4),
                number(row["ratio_ci_low"], 4), number(row["ratio_ci_high"], 4),
                latex_escape(row["FT_status"]),
            )
        )
    parts.append(r"\bottomrule\end{longtable}")

    parts.append(r"\section*{Crossover summary}")
    parts.append(r"\scriptsize\begin{longtable}{rrrrrrrrl}")
    parts.append(
        r"\toprule $\Delta\beta$ & $a_\infty/\Delta\beta$ & ratio CI & "
        r"$a_G/\Delta\beta$ & $G_{FT}(640)$ & skew(160) & ex.kurt(160) & "
        r"$n_-(160)$ & status \\"
    )
    parts.append(r"\midrule\endhead")
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
    parts.append(r"\bottomrule\end{longtable}\normalsize")

    parts.append(r"""\section*{Figures}
\begin{figure}[ht]
\centering
\includegraphics[width=0.76\textwidth]{../figures/a_fit_vs_inverse_time.pdf}
\caption{Direct matched-bin slopes against inverse window duration.}
\end{figure}
\begin{figure}[ht]
\centering
\includegraphics[width=0.76\textwidth]{../figures/slope_crossover_vs_affinity.pdf}
\caption{Extrapolated direct-tail slopes and Gaussian bulk slopes against affinity.}
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
are in \texttt{analysis/window\_summary.csv}.  A case labelled
\texttt{CONSISTENT\_WITH\_FT} means only that the predeclared extrapolated
confidence interval includes the reference.  It is not a proof.  A case whose
interval excludes the reference is labelled \texttt{FAIL}; insufficient raw
two-tail support remains \texttt{UNRESOLVED}.  The label
\texttt{CONSISTENT\_WITH\_EQUILIBRIUM} has the corresponding limited meaning
for the zero-affinity control.
\end{document}
""")

    args.output.write_text("\n".join(parts))


if __name__ == "__main__":
    main()
