#!/usr/bin/env python3
"""Generate auditable LaTeX tables directly from frozen CSV/JSON outputs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ANA = ROOT / "analysis"
OUT = ROOT / "report" / "generated_tables.tex"
CASES = [
    ("driven_dt1e-3", r"driven, $10^{-3}$"),
    ("driven_dt2p5e-4", r"driven, $2.5\!\times\!10^{-4}$"),
    ("equilibrium_dt1e-3", r"equal $T$, $10^{-3}$"),
    ("equilibrium_dt2p5e-4", r"equal $T$, $2.5\!\times\!10^{-4}$"),
]
OBS = {
    "cos_theta1": r"$\cos\theta_1$",
    "sin_theta1": r"$\sin\theta_1$",
    "cos_theta3": r"$\cos\theta_3$",
    "sin_theta3": r"$\sin\theta_3$",
    "cos_theta1_minus_theta3": r"$\cos(\theta_1-\theta_3)$",
    "I2": r"$I_2$",
    "I1_plus_I3": r"$I_1+I_3$",
    "I1_minus_I3": r"$I_1-I_3$",
    "cos_even_sum": r"$\cos\theta_1+\cos\theta_3$",
    "cos_odd_diff": r"$\cos\theta_1-\cos\theta_3$",
    "sin_even_sum": r"$\sin\theta_1+\sin\theta_3$",
    "sin_odd_diff": r"$\sin\theta_1-\sin\theta_3$",
}


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def f4(x):
    return f"{float(x):.4f}"


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append(r"\begin{longtable}{lllrr}")
    lines.append(r"\caption{Frozen-rule plateau results for every saved observable.}\label{tab:allrates}\\")
    lines.append(r"\toprule Case & Observable & Window & $\widehat\lambda_R$ [95\% CI] & $N_{\rm eff}$ \\")
    lines.append(r"\midrule\endfirsthead")
    lines.append(r"\toprule Case & Observable & Window & $\widehat\lambda_R$ [95\% CI] & $N_{\rm eff}$ \\")
    lines.append(r"\midrule\endhead")
    for case, label in CASES:
        for row in read(ANA / f"{case}_plateaus.csv"):
            if row["lambda_R"]:
                window = f"[{float(row['plateau_start']):.2f},{float(row['plateau_end']):.2f}]"
                rate = f"{f4(row['lambda_R'])} [{f4(row['ci_low'])},{f4(row['ci_high'])}]"
            else:
                window = "--"
                rate = r"\multicolumn{1}{c}{UNRESOLVED}"
            lines.append(
                f"{label} & {OBS[row['observable']]} & {window} & {rate} & "
                f"{float(row['effective_count']):.0f} \\\\"
            )
        lines.append(r"\addlinespace")
    lines.append(r"\bottomrule\end{longtable}")

    summary = json.load((ANA / "controlled_summary.json").open())
    lines.append(r"\begin{table}[t]\centering\small")
    lines.append(r"\caption{Pooled common-plateau gate.  The eligible column lists observables that individually passed the frozen plateau rule.}\label{tab:common}")
    lines.append(r"\begin{tabular}{lp{4.1cm}ccc}\toprule Case & Eligible & Spread pass & CI intersection & Verdict\\\midrule")
    for item in summary["case_summaries"]:
        common = item["common_plateau"]
        eligible = ", ".join(OBS[x] for x in common.get("eligible", [])) or "none"
        lines.append(
            f"{dict(CASES)[item['case']]} & {eligible} & "
            f"{'yes' if common.get('rate_spread_ok', False) else 'no'} & "
            f"{'yes' if common.get('ci_overlap', False) else 'no'} & NO COMMON PLATEAU \\\\"
        )
    lines.append(r"\bottomrule\end{tabular}\end{table}")

    lines.append(r"\begin{longtable}{llll}")
    lines.append(r"\caption{Historical versus controlled rates. Blank numerical entries mean the observable was unresolved; $\dagger$ means it was not saved historically.}\label{tab:comparison}\\")
    lines.append(r"\toprule Bath & Observable & Historical & Controlled: $10^{-3}$; $2.5\times10^{-4}$\\\midrule\endfirsthead")
    lines.append(r"\toprule Bath & Observable & Historical & Controlled: $10^{-3}$; $2.5\times10^{-4}$\\\midrule\endhead")
    for row in read(ANA / "previous_vs_controlled_rates.csv"):
        hist = r"UNRESOLVED" if row["previous_status"] == "UNRESOLVED" else (
            r"not saved$^\dagger$" if row["previous_status"] == "NOT_SAVED"
            else f"{f4(row['previous_lambda_R'])} [{f4(row['previous_ci_low'])},{f4(row['previous_ci_high'])}]"
        )
        def new(prefix):
            return (
                f"{f4(row[prefix + '_lambda_R'])} [{f4(row[prefix + '_ci_low'])},{f4(row[prefix + '_ci_high'])}]"
                if row[prefix + "_lambda_R"] else "UNRESOLVED"
            )
        lines.append(
            f"{row['bath_case']} & {OBS[row['observable']]} & {hist} & "
            f"{new('new_dt1e-3')}; {new('new_dt2p5e-4')} \\\\"
        )
    lines.append(r"\bottomrule\end{longtable}")

    lines.append(r"\begin{longtable}{lllrrrl}")
    lines.append(r"\caption{Predeclared early-window damped fits, all using $t\in[0.05,1.00]$.}\label{tab:odd}\\")
    lines.append(r"\toprule Case & Observable & $\lambda_R$ [95\% CI] & $\lambda_I$ & 95\% CI & Period & $\Delta$AIC\\\midrule\endfirsthead")
    lines.append(r"\toprule Case & Observable & $\lambda_R$ [95\% CI] & $\lambda_I$ & 95\% CI & Period & $\Delta$AIC\\\midrule\endhead")
    for row in read(ANA / "early_odd_damped_fits.csv"):
        lines.append(
            f"{dict(CASES)[row['case']]} & {OBS[row['observable']]} & "
            f"{f4(row['lambda_R'])} [{f4(row['lambda_R_ci_low'])},{f4(row['lambda_R_ci_high'])}] & "
            f"{f4(row['lambda_I'])} & [{f4(row['lambda_I_ci_low'])},{f4(row['lambda_I_ci_high'])}] & "
            f"{float(row['period']):.3f} & {float(row['delta_AIC_damped_minus_exponential']):.1f} \\\\"
        )
    lines.append(r"\bottomrule\end{longtable}")

    lines.append(r"\begin{table}[t]\centering\small")
    lines.append(r"\caption{Run-level numerical and stationarity audit.  $N_{\rm eff,last}$ is the smallest effective count among accepted plateaus after accounting for the terminal plateau lag.}\label{tab:audit}")
    lines.append(r"\begin{tabular}{lrrrrr}\toprule Case & Proj./floor & Midpoint fail & Nonfinite & Worst burn-in metric & $N_{\rm eff,last}$\\\midrule")
    audit = {x["case"]: x for x in json.load((ROOT / "provenance" / "integrity_audit.json").open())["cases"]}
    for case, label in CASES:
        rows = read(ANA / f"{case}_plateaus.csv")
        worst = max(float(r["stationarity_rms_over_C0"]) for r in rows)
        neff = []
        for r in rows:
            if r["plateau_end"]:
                neff.append(64 * (1000 - float(r["plateau_end"])) / (2 * float(r["tau_int"])))
        a = audit[case]
        lines.append(
            f"{label} & {a['projection_count']}/{a['floor_count']} & "
            f"{a['midpoint_failure_count']} & {a['nonfinite']} & {worst:.5f} & {min(neff):.0f} \\\\"
        )
    lines.append(r"\bottomrule\end{tabular}\end{table}")

    OUT.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
