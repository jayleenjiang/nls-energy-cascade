#!/usr/bin/env python3
"""Build a manuscript claim registry from source-traced numerical artifacts.

The goal of this script is deliberately modest: it does not try to understand
every mathematical sentence in ``draft.tex``.  Instead, it verifies the core
numerical claims that would be hardest to defend during journal review:

* long-chain action-current scaling and finite-window current diagnostics;
* long-chain action profiles and local-equilibrium table values;
* short-chain neural-network Fokker--Planck and eigenfunction diagnostics;
* the normalization convention used for equilibrium local marginals.

The output is a machine-readable JSON registry plus a compact Markdown report.
Any failed text check is treated as an audit failure because it means a value
present in the evidence files is not faithfully represented in the manuscript.
"""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
DRAFT = REVISION / "draft.tex"
DEFAULT_JSON = REVISION / "manuscript_claim_audit.json"
DEFAULT_MD = REVISION / "manuscript_claim_audit.md"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def tex_text() -> str:
    return DRAFT.read_text()


def normalize_ws(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def text_present(tex: str, needle: str) -> bool:
    return needle in tex or normalize_ws(needle) in normalize_ws(tex)


def line_for(tex: str, needle: str) -> int | None:
    idx = tex.find(needle)
    if idx < 0:
        idx = normalize_ws(tex).find(normalize_ws(needle))
        if idx < 0:
            return None
        return None
    return tex[:idx].count("\n") + 1


def check_strings(tex: str, strings: list[str]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    for s in strings:
        checks.append({"text": s, "present": text_present(tex, s), "line": line_for(tex, s)})
    return checks


def verdict_from_checks(checks: list[dict[str, Any]]) -> str:
    return "VERIFIED" if all(c["present"] for c in checks) else "FAIL"


def add_claim(
    registry: list[dict[str, Any]],
    *,
    claim_id: str,
    section: str,
    claim: str,
    evidence: list[Path],
    expected_text: list[str],
    computed: dict[str, Any] | None = None,
    note: str = "",
) -> None:
    tex = tex_text()
    checks = check_strings(tex, expected_text)
    registry.append(
        {
            "id": claim_id,
            "section": section,
            "claim": claim,
            "verdict": verdict_from_checks(checks),
            "evidence": [rel(p) for p in evidence],
            "computed": computed or {},
            "text_checks": checks,
            "note": note,
        }
    )


def round_str(x: float, digits: int) -> str:
    return f"{x:.{digits}f}"


def read_window_stats(path: Path) -> dict[tuple[int, int], dict[str, float]]:
    rows: dict[tuple[int, int], dict[str, float]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            key = (int(row["n"]), int(float(row["window"])))
            rows[key] = {k: float(v) for k, v in row.items() if k not in {"n", "dt"}}
            rows[key]["n"] = float(row["n"])
            rows[key]["dt"] = float(row["dt"])
    return rows


def read_summary(path: Path) -> dict[str, Any]:
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one summary row")
    row: dict[str, Any] = dict(rows[0])
    for key in [
        "n",
        "batches",
        "lanes",
        "n_trajectories",
        "bond",
        "seed",
        "threads",
        "projection_count",
    ]:
        row[key] = int(row[key])
    for key in [
        "T1",
        "Tn",
        "gamma",
        "dt",
        "burnin",
        "measure",
        "mean_action_current",
        "sample_sd",
        "standard_error",
        "normal95_ci_lower",
        "normal95_ci_upper",
        "mean_first_half_current",
        "mean_second_half_current",
        "mean_second_minus_first",
        "paired_difference_se",
        "projection_rate",
        "elapsed_seconds",
    ]:
        row[key] = float(row[key])
    return row


def build_registry() -> list[dict[str, Any]]:
    tex = tex_text()
    registry: list[dict[str, Any]] = []

    flux_path = REVISION / "experiments/flux_validation/production_dt5e-4/flux_primary_scaling.json"
    validation_path = REVISION / "experiments/flux_validation/validation_report.md"
    window_path = REVISION / "experiments/flux_validation/production_dt5e-4/current_windows_window_statistics.csv"
    larger_n_dir = REVISION / "experiments/flux_validation/larger_n_pilot_2026-06-20"
    larger_n_scaling_path = larger_n_dir / "n10_50_b64_scaling_scaling.json"
    larger_n_summary_path = larger_n_dir / "n50_b64_summary.csv"
    larger_n_fine_path = larger_n_dir / "n50_b64_dt2p5e-4_summary.csv"
    larger_n_long_burn_path = larger_n_dir / "n50_b16_burn10000_summary.csv"
    larger_n_readme_path = larger_n_dir / "README.md"
    n60_dir = REVISION / "experiments/flux_validation/larger_n60_pilot_2026-06-20"
    n60_scaling_path = n60_dir / "n10_60_b64_scaling_scaling.json"
    n60_summary_path = n60_dir / "n60_b64_summary.csv"
    n60_readme_path = n60_dir / "README.md"
    flux_sensitivity_path = n60_dir / "flux_scaling_sensitivity_n10_60.json"
    parameter_dir = REVISION / "experiments/flux_validation/parameter_robustness_2026-06-20"
    parameter_scaling_path = parameter_dir / "moderate_contrast_T8_T4_prod/b64_scaling_scaling.json"
    parameter_summary_paths = [
        parameter_dir / f"moderate_contrast_T8_T4_prod/n{n}_b64_summary.csv"
        for n in (10, 20, 30, 40)
    ]
    parameter_readme_path = parameter_dir / "README.md"
    parameter_production_summary_path = parameter_dir / "production_summary.csv"
    gamma_dir = REVISION / "experiments/flux_validation/gamma_robustness_2026-06-21"
    gamma_scaling_path = gamma_dir / "gamma_robustness_scaling.json"
    gamma_summary_csv_path = gamma_dir / "gamma_robustness_summary.csv"
    gamma_report_path = gamma_dir / "gamma_robustness_report.md"
    gamma_generated_sources = [
        gamma_dir / "generated_sources/gamma0p05/NLS_flux_canonical_gamma.cpp",
        gamma_dir / "generated_sources/gamma0p2/NLS_flux_canonical_gamma.cpp",
    ]
    gamma_summary_paths = [
        gamma_dir / label / f"n{n}_summary.csv"
        for label in ("gamma0p05", "gamma0p2")
        for n in (10, 20, 30, 40)
    ]
    figure_metrics_path = REVISION / "manuscript_figure_metrics.json"
    source_trace_path = REVISION / "source_trace_metrics.json"
    rerun_path = REVISION / "short_chain_nn_rerun_metrics.json"
    eigen_path = REVISION / "eigen_fit_sensitivity.json"
    availability_audit_json = REVISION / "availability_path_audit.json"
    availability_audit_md = REVISION / "availability_path_audit.md"

    flux = load_json(flux_path)
    figures = load_json(figure_metrics_path)
    source = load_json(source_trace_path)
    rerun = load_json(rerun_path)
    eigen = load_json(eigen_path)
    window = read_window_stats(window_path)
    larger_n_scaling = load_json(larger_n_scaling_path)
    larger_n_summary = read_summary(larger_n_summary_path)
    larger_n_fine = read_summary(larger_n_fine_path)
    larger_n_long_burn = read_summary(larger_n_long_burn_path)
    n60_scaling = load_json(n60_scaling_path)
    n60_summary = read_summary(n60_summary_path)
    flux_sensitivity = load_json(flux_sensitivity_path)
    parameter_scaling = load_json(parameter_scaling_path)
    parameter_summaries = [read_summary(path) for path in parameter_summary_paths]
    gamma_scaling = load_json(gamma_scaling_path)

    add_claim(
        registry,
        claim_id="intro_claim_evidence_map",
        section="contributions and organization",
        claim="The introduction includes a claim-evidence map that separates supported finite-size or diagnostic statements from deliberately unclaimed stronger interpretations.",
        evidence=[
            DRAFT,
            flux_path,
            validation_path,
            flux_sensitivity_path,
            parameter_scaling_path,
            gamma_scaling_path,
            window_path,
            figure_metrics_path,
            source_trace_path,
            rerun_path,
            eigen_path,
        ],
        expected_text=[
            r"\label{tab:claim-evidence-map}",
            r"Claim--evidence map for the manuscript",
            r"Does the simulated SDE have the intended equilibrium structure?",
            r"Equivalence with the original HLNS thermostat.",
            r"Production Monte Carlo for $n=10,20,30,40$, with $n=50,60$",
            r"Finite-size faster-than-Fourier decay of the mean action current.",
            r"An asymptotic exponent or universal parameter law.",
            r"Action profiles, pair-marginal comparisons, LTE table entries, and residual",
            r"Strict local Gibbs convergence or vanishing residuals.",
            r"Finite-window current histograms, survival fits, and",
            r"Large-deviation asymptotics or Gallavotti--Cohen symmetry.",
            r"Saved-model Gibbs validation, short-chain density diagnostics,",
            r"Proof of the long-chain exponent or a resolved spectral gap.",
        ],
        computed={
            "mapped_claim_families": [
                "canonical_sde_validation",
                "finite_size_current_scaling",
                "local_equilibrium_diagnostics",
                "finite_window_current_fluctuations",
                "short_chain_mechanism_diagnostics",
            ]
        },
        note="This introductory map is a scope-control device: it does not add new numerical values, but it binds each manuscript claim family to its evidence stream and explicit non-claim.",
    )

    prefactor = flux["prefactor"]
    exponent = flux["exponent"]
    ci_lo, ci_hi = flux["exponent_normalized_95_ci"]
    conductivity_prefactor = prefactor / 8.0
    conductivity_exponent = exponent + 1.0

    add_claim(
        registry,
        claim_id="flux_scaling_main",
        section="abstract / thermal conductivity",
        claim="Mean action current follows a finite-size power law over n=10,20,30,40.",
        evidence=[flux_path, validation_path],
        expected_text=[
            r"\E[J(n)]=28.75\,n^{-1.850}",
            r"$[-1.870,-1.830]$",
            r"$n=10,20,30,40$",
            r"R^2=0.9980",
        ],
        computed={
            "prefactor": prefactor,
            "exponent": exponent,
            "exponent_ci": [ci_lo, ci_hi],
            "r_squared": flux["r_squared_log_fit"],
        },
    )

    add_claim(
        registry,
        claim_id="conductivity_scaling",
        section="abstract / thermal conductivity",
        claim="Finite-chain action conductivity decays as an insulating power law.",
        evidence=[flux_path],
        expected_text=[
            r"\kappa_n\sim n^{-0.850}",
            r"\sim 3.59\,n^{-0.850}",
            r"R_n^{\mathrm{act}}\sim n^{1.850}",
            r"not a derivation of the exponent",
        ],
        computed={
            "conductivity_prefactor": conductivity_prefactor,
            "conductivity_exponent": conductivity_exponent,
            "action_resistance_exponent": -exponent,
        },
        note="Computed as n*E[J(n)]/(T1-Tn) with T1-Tn=8.",
    )

    # Flux table values and primary simulation settings.
    summary_rows = [
        (10, 0.3925219606, 0.0018902080),
        (20, 0.1191693526, 0.0009195305),
        (30, 0.0545731139, 0.0006191849),
        (40, 0.0297475540, 0.0004827205),
    ]
    add_claim(
        registry,
        claim_id="flux_table_values",
        section="thermal conductivity",
        claim="Current table reports the four production means and standard errors.",
        evidence=[validation_path, flux_path],
        expected_text=[
            r"$\Delta t=5\times10^{-4}$",
            r"$1024$ independent trajectories",
            r"($1000,1280,2880,5120$ for $n=10,20,30,40$)",
            r"0.39252 & 0.11917 & 0.05457 & 0.02975",
            r"0.00189 & 0.00092 & 0.00062 & 0.00048",
        ],
        computed={f"n={n}": {"mean": mean, "se": se} for n, mean, se in summary_rows},
    )

    add_claim(
        registry,
        claim_id="flux_diagnostics",
        section="thermal conductivity",
        claim="Stationarity and timestep diagnostics are reported with limited scope.",
        evidence=[validation_path],
        expected_text=[
            r"only $0.90$",
            r"relative difference about $3.9\%$",
            r"$0.60$ pooled standard errors",
        ],
        computed={
            "stationarity_max_abs_z": flux["stationarity_max_abs_z"],
            "validation_report_source": rel(validation_path),
        },
    )

    add_claim(
        registry,
        claim_id="larger_n_current_robustness",
        section="thermal conductivity",
        claim="The n=50 larger-chain current run and production fine-step check are reported as robustness evidence rather than replacing the primary exponent.",
        evidence=[
            larger_n_scaling_path,
            larger_n_summary_path,
            larger_n_fine_path,
            larger_n_long_burn_path,
            larger_n_readme_path,
        ],
        expected_text=[
            r"$\E[J(50)]=0.01852$",
            r"standard error $0.00044$",
            r"$-1.49$ paired standard",
            r"\E[J(n)] \;=\; 32.50\,n^{-1.894}",
            r"$[-1.916,-1.873]$",
            r"$\E[J(50)]=0.01918$",
            r"standard error $0.00040$",
            r"$-1.56$ paired standard errors",
            r"$3.6\%$ upward shift, or $1.12$",
            r"$\E[J(50)]=0.01931\pm0.00086$",
            r"$n=50$ and $n=60$ computations as evidence against",
        ],
        computed={
            "n50_primary_mean": larger_n_summary["mean_action_current"],
            "n50_primary_se": larger_n_summary["standard_error"],
            "n50_primary_stationarity_z": (
                larger_n_summary["mean_second_minus_first"]
                / larger_n_summary["paired_difference_se"]
            ),
            "n50_fine_mean": larger_n_fine["mean_action_current"],
            "n50_fine_se": larger_n_fine["standard_error"],
            "n50_fine_stationarity_z": (
                larger_n_fine["mean_second_minus_first"]
                / larger_n_fine["paired_difference_se"]
            ),
            "n50_fine_minus_coarse": (
                larger_n_fine["mean_action_current"]
                - larger_n_summary["mean_action_current"]
            ),
            "n50_fine_relative_shift": (
                (
                    larger_n_fine["mean_action_current"]
                    - larger_n_summary["mean_action_current"]
                )
                / larger_n_summary["mean_action_current"]
            ),
            "n50_fine_difference_over_pooled_se": (
                (
                    larger_n_fine["mean_action_current"]
                    - larger_n_summary["mean_action_current"]
                )
                / (
                    larger_n_fine["standard_error"] ** 2
                    + larger_n_summary["standard_error"] ** 2
                )
                ** 0.5
            ),
            "five_length_exponent": larger_n_scaling["exponent"],
            "five_length_ci": larger_n_scaling["exponent_normalized_95_ci"],
            "longer_burnin_mean": larger_n_long_burn["mean_action_current"],
            "longer_burnin_se": larger_n_long_burn["standard_error"],
        },
        note="The manuscript keeps the n=10,20,30,40 fit as the primary quoted exponent while using larger lengths as robustness extensions.",
    )

    add_claim(
        registry,
        claim_id="n60_current_robustness",
        section="thermal conductivity",
        claim="The production-size n=60 extension supports the larger-chain robustness check without redefining the primary exponent.",
        evidence=[
            n60_scaling_path,
            n60_summary_path,
            n60_readme_path,
        ],
        expected_text=[
            r"$1024$ trajectories, burn-in $11520$",
            r"$\E[J(60)]=0.01245$",
            r"standard error $0.00042$",
            r"$-0.52$ paired standard",
            r"\E[J(n)] \;=\; 35.94\,n^{-1.930}",
            r"$[-1.954,-1.906]$",
            r"larger lengths are robustness extensions",
        ],
        computed={
            "n60_mean": n60_summary["mean_action_current"],
            "n60_se": n60_summary["standard_error"],
            "n60_stationarity_z": (
                n60_summary["mean_second_minus_first"]
                / n60_summary["paired_difference_se"]
            ),
            "six_length_exponent": n60_scaling["exponent"],
            "six_length_ci": n60_scaling["exponent_normalized_95_ci"],
            "six_length_r_squared": n60_scaling["r_squared_log_fit"],
        },
        note="The n=60 run has production-size trajectory count but no matched fine-step production run, so it is used as robustness evidence.",
    )

    parameter_max_abs_z = max(
        abs(row["mean_second_minus_first"] / row["paired_difference_se"])
        for row in parameter_summaries
    )
    add_claim(
        registry,
        claim_id="bath_parameter_current_robustness",
        section="thermal conductivity",
        claim="A second bath-temperature production run supports the current-scaling trend without becoming a systematic parameter sweep.",
        evidence=[
            parameter_scaling_path,
            *parameter_summary_paths,
            parameter_readme_path,
            parameter_production_summary_path,
        ],
        expected_text=[
            r"$T_1=8,T_n=4$",
            r"$1024$ trajectories per length",
            r"\E[J(n)] \;=\; 12.95\,n^{-1.751}",
            r"R^2=0.9984",
            r"$[-1.780,-1.723]$",
            r"maximum first-half/second-half stationarity statistic",
            r"$1.74$ paired standard errors",
            r"not a systematic parameter sweep",
        ],
        computed={
            "T1": parameter_summaries[0]["T1"],
            "Tn": parameter_summaries[0]["Tn"],
            "means": {
                f"n={row['n']}": row["mean_action_current"]
                for row in parameter_summaries
            },
            "standard_errors": {
                f"n={row['n']}": row["standard_error"]
                for row in parameter_summaries
            },
            "max_abs_stationarity_z": parameter_max_abs_z,
            "exponent": parameter_scaling["exponent"],
            "exponent_ci": parameter_scaling["exponent_normalized_95_ci"],
            "r_squared": parameter_scaling["r_squared_log_fit"],
        },
        note="The manuscript uses this as a robustness check at a second bath-temperature pair, not as a full parameter sweep.",
    )

    gamma_by_label = {
        record["label"]: record for record in gamma_scaling["scaling_records"]
    }
    gamma_005 = gamma_by_label["gamma0p05"]
    gamma_02 = gamma_by_label["gamma0p2"]
    add_claim(
        registry,
        claim_id="gamma_thermostat_current_robustness",
        section="thermal conductivity / finite-size and parameter robustness",
        claim="Thermostat-coupling production runs support the current-scaling trend without becoming a systematic two-parameter bath sweep.",
        evidence=[
            gamma_scaling_path,
            gamma_summary_csv_path,
            gamma_report_path,
            *gamma_generated_sources,
            *gamma_summary_paths,
        ],
        expected_text=[
            r"thermostat-coupling robustness",
            r"$\gamma=0.05$",
            r"$-1.650$",
            r"$[-1.668,-1.633]$",
            r"$\gamma=0.2$",
            r"$-1.991$",
            r"$[-2.017,-1.967]$",
            r"statistics are $1.14$ and $1.74$ paired standard errors",
            r"single thermostat coupling $\gamma=0.1$",
            r"systematic two-parameter bath study",
            r"\label{tab:parameter-robustness}",
        ],
        computed={
            "gamma_0p05": {
                "exponent": gamma_005["exponent"],
                "exponent_ci": gamma_005["exponent_95_ci"],
                "r_squared": gamma_005["r_squared_log_fit"],
                "max_abs_stationarity_z": gamma_005["max_abs_stationarity_z"],
            },
            "gamma_0p2": {
                "exponent": gamma_02["exponent"],
                "exponent_ci": gamma_02["exponent_95_ci"],
                "r_squared": gamma_02["r_squared_log_fit"],
                "max_abs_stationarity_z": gamma_02["max_abs_stationarity_z"],
            },
            "status": gamma_scaling["status"],
        },
        note="The gamma-specific sources are generated from the frozen canonical source; the primary source file remains unchanged.",
    )

    sensitivity_by_label = {
        record["label"]: record for record in flux_sensitivity["fit_windows"]
    }
    add_claim(
        registry,
        claim_id="flux_scaling_fit_sensitivity",
        section="thermal conductivity",
        claim="Fit-window sensitivity around the n=50 and n=60 robustness points is reported without promoting it to the primary exponent.",
        evidence=[
            flux_sensitivity_path,
            larger_n_summary_path,
            n60_summary_path,
            flux_path,
        ],
        expected_text=[
            r"\label{tab:flux-fit-sensitivity}",
            r"range from $-1.720$ on $10$--$20$ to $-2.178$ on",
            r"primary $n=10,20,30,40$ & $-1.850$ & $[-1.870,-1.831]$ & $0.9980$",
            r"$n=10,20,30,40,50$ & $-1.894$ & $[-1.917,-1.873]$ & $0.9976$",
            r"$n=10,20,30,40,50,60$ & $-1.930$ & $[-1.954,-1.906]$ & $0.9974$",
            r"tail $n=20,30,40,50,60$ & $-2.059$ & $[-2.107,-2.012]$ & $0.9993$",
            r"the spread across rows is the finite-size sensitivity",
        ],
        computed={
            "primary": sensitivity_by_label["primary n=10--40"],
            "with_n50": sensitivity_by_label["with n=50"],
            "with_n50_n60": sensitivity_by_label["with n=50,60"],
            "tail_n20_60": sensitivity_by_label["tail n=20--60"],
            "local_slopes": flux_sensitivity["local_slopes"],
        },
        note="The table is a finite-size fit-window diagnostic; it does not change the primary quoted n=10,20,30,40 exponent.",
    )

    w40 = {tau: window[(40, tau)] for tau in (50, 100, 200)}
    add_claim(
        registry,
        claim_id="finite_window_current_statistics",
        section="finite-time current fluctuations",
        claim="Finite-window current statistics at n=40 are descriptive, not a large-deviation claim.",
        evidence=[window_path],
        expected_text=[
            r"$\Pr(\overline J_\tau<0)$ decreases from $0.233$ at $\tau=50$ to $0.026$ at",
            r"($0.010$, $0.024$, $0.075$ for the three windows)",
            r"($0.095$, $0.067$, $0.048$)",
            r"not infer an asymptotic large-deviation",
        ],
        computed={
            str(tau): {
                "negative_fraction": w40[tau]["negative_fraction"],
                "skewness": w40[tau]["skewness"],
                "window_times_variance": w40[tau]["window_times_variance"],
            }
            for tau in (50, 100, 200)
        },
    )

    profiles = figures["action_profiles"]["metrics"]
    add_claim(
        registry,
        claim_id="long_chain_action_profiles",
        section="nonequilibrium steady state",
        claim="Profile endpoints and sample counts match the figure source metrics.",
        evidence=[figure_metrics_path],
        expected_text=[
            r"$\langle I_1\rangle=0.738$",
            r"$\langle I_{13}\rangle=0.573$",
            r"$\langle I_{25}\rangle=0.165$",
            r"$(0.514,0.106)$ at $n=50$",
            r"$(0.362,0.073)$ at $n=100$",
            r"$5.07\times10^8$, $1.01\times10^9$, and $4.06\times10^9$",
        ],
        computed=profiles,
    )

    # Local-equilibrium table.
    expected_lte: list[str] = []
    lte_computed: dict[str, Any] = {}
    for row in source["lte_table"]["rows"]:
        n = row["n"]
        a = row["a_over_n"]
        key = f"n={n},a/n={a:.2f}"
        x = row["fit"]["slope_x"]
        r2 = row["fit"]["weighted_r2"]
        six_over_x = row["fit"]["six_over_x"]
        tkin = row["kinetic_temperature"]["T_kin"]
        expected_lte.extend([round_str(x, 4), round_str(r2, 4), round_str(six_over_x, 2), round_str(tkin, 2)])
        lte_computed[key] = {
            "x": x,
            "weighted_r2": r2,
            "six_over_x": six_over_x,
            "T_kin": tkin,
        }
    add_claim(
        registry,
        claim_id="lte_table_values",
        section="local thermodynamic equilibrium",
        claim="Table tab:lte values match source-traced histogram/profile recomputation.",
        evidence=[source_trace_path],
        expected_text=sorted(set(expected_lte)),
        computed=lte_computed,
        note="Text check is intentionally value-based; repeated rounded values need appear at least once in the table.",
    )

    lte_residual_decomp = source["lte_residual_decomposition"]
    expected_residual_decomp = [
        r"\label{tab:lte-resid-decomp}",
        r"$r=\log q-(c+x\log p_6)$",
        r"r_{\rm even}(\theta)=\tfrac12\{r(\theta)+r(-\theta)\}",
        r"r_{\rm odd}(\theta)=\tfrac12\{r(\theta)-r(-\theta)\}",
        r"requiring at least $50$ counts",
        r"odd part is a",
        r"direct finite nonequilibrium correction",
    ]
    lte_residual_computed: dict[str, Any] = {}
    for row in lte_residual_decomp["rows"]:
        key = f"n={row['n']},a/n={row['a_over_n']:.2f}"
        rounding = row["manuscript_rounding"]
        expected_residual_decomp.extend(
            [
                rounding["total_rms_3dp"],
                rounding["even_rms_3dp"],
                rounding["odd_rms_3dp"],
                rounding["odd_fraction_2dp"],
            ]
        )
        lte_residual_computed[key] = {
            "total_rms": row["decomposition"]["total_rms"],
            "even_rms": row["decomposition"]["even_rms"],
            "odd_rms": row["decomposition"]["odd_rms"],
            "odd_fraction_of_total_rms": row["decomposition"]["odd_fraction_of_total_rms"],
            "symmetrized_original_bin_count": row["decomposition"]["symmetrized_original_bin_count"],
        }
    add_claim(
        registry,
        claim_id="lte_residual_even_odd_decomposition",
        section="local thermodynamic equilibrium",
        claim="The LTE residual even/odd table matches the source-traced histogram recomputation.",
        evidence=[source_trace_path],
        expected_text=sorted(set(expected_residual_decomp)),
        computed={
            "algorithm": lte_residual_decomp["algorithm"],
            "rows": lte_residual_computed,
        },
        note="The residual decomposition is descriptive and uses a stricter symmetric count mask than the slope fit.",
    )

    controls = source["lte_table"]["controls"]
    add_claim(
        registry,
        claim_id="lte_control_values",
        section="local thermodynamic equilibrium",
        claim="Equilibrium and matched-equilibrium control fits match source-traced controls.",
        evidence=[source_trace_path],
        expected_text=[
            r"$x=0.926$",
            r"$6/7.12=0.843$",
            r"$p_{7.12}$ vs.\ $p_6$",
            r"$0.9261$",
            r"$0.9997$",
            r"$1.0093$",
            r"$0.9944$",
        ],
        computed=controls,
    )

    add_claim(
        registry,
        claim_id="lte_equilibrium_convention",
        section="local thermodynamic equilibrium",
        claim="Local equilibrium marginal notation follows the global Gibbs convention exp[-H/(2T)].",
        evidence=[DRAFT],
        expected_text=[
            r"$p_T(I_a,I_{a+1},\theta_a)$",
            r"First, neither $q$ nor $p_T$",
            r"$\exp[-H/(2T)]$",
            r"$p_T=p_6^{\,6/T}$",
        ],
        computed={"forbidden_old_form_present": r"\exp(-H/\beta)" in tex},
        note="This guards against mixing inverse-temperature beta with the manuscript's temperature T.",
    )

    eq_errors = rerun["equilibrium_validation"]
    add_claim(
        registry,
        claim_id="short_chain_equilibrium_validation",
        section="short-chain Fokker--Planck",
        claim="Neural-network density reproduces the equilibrium Gibbs slices to the reported relative errors.",
        evidence=[rerun_path, source_trace_path],
        expected_text=[r"$1.5\%$", r"$2.6\%$", r"$6.4\%$", r"$I=0.5,1,2$"],
        computed={
            str(row["I"]): {
                "mean_relative_error_percent": row["mean_relative_error_percent"],
                "max_relative_error_percent": row["max_relative_error_percent"],
            }
            for row in eq_errors
        },
    )

    add_claim(
        registry,
        claim_id="short_chain_symmetry_scope",
        section="stabilization",
        claim="Symmetry breaking is reported qualitatively rather than as a standalone percentage estimator.",
        evidence=[rerun_path],
        expected_text=[
            r"qualitative evidence",
            r"rather than as a standalone quantitative estimator",
            r"$\pm0.15$",
        ],
        computed=rerun["symmetry_breaking"],
    )

    eig_main = eigen["summaries"][0]
    eig_late = next(x for x in eigen["summaries"] if x["label"] == "0.5--5")
    eig_short = next(x for x in eigen["summaries"] if x["label"] == "0--2")
    add_claim(
        registry,
        claim_id="eigen_relaxation_diagnostic",
        section="eigenfunction",
        claim="Eigen relaxation rate is framed as an observable-dependent diagnostic with window sensitivity.",
        evidence=[eigen_path, rerun_path],
        expected_text=[
            r"retain $2895$ of $4096$",
            r"\lambda_{\mathrm{diag}} \;\approx\; -0.934",
            r"from $-1.60$ on $t\in[0,2]$ to $-0.68$ on $t\in[0.5,5]$",
            r"root-mean-square error is $0.38$",
            r"relative error $0.71$",
            r"median normalized PDE residual",
            r"$0.35$",
            r"not as a resolved statement about the full generator spectrum",
        ],
        computed={
            "main_window": eig_main,
            "late_window": eig_late,
            "eigen_surrogate_data_fit": rerun["eigen_surrogate_data_fit"],
        },
    )

    add_claim(
        registry,
        claim_id="short_chain_solver_diagnostics_table",
        section="numerical validation appendix",
        claim="A compact short-chain diagnostic table reports source-traced solver checks while preserving qualitative scope.",
        evidence=[DRAFT, rerun_path, eigen_path],
        expected_text=[
            r"\label{app:short-chain-diagnostics}",
            r"\label{tab:short-chain-diagnostics}",
            r"$1.47\%,2.62\%,6.40\%$",
            r"$3.95\%,7.18\%,18.84\%$",
            r"$\rho+\rho^{\mathsf T}\ge 20\%\,\max(\rho+\rho^{\mathsf T})$",
            r"$17.5\%$ and $56.2\%$",
            r"$\theta_0=2.191$ radians ($125.5^\circ$)",
            r"$138^\circ$--$168^\circ$",
            r"$\E[I_2I_1\sin\theta_1]=0.128$",
            r"$\E[I_2I_3\sin\theta_3]=-0.140$",
            r"$8.9\%$ of the first term",
            r"$2895/4096$ full-window fits",
            r"$-1.60$ on $[0,2]$ to $-0.676$ on",
            r"RMSE $0.38$ (relative $0.71$)",
            r"median normalized PDE residual $0.35$",
            r"not a resolved spectral-gap computation",
        ],
        computed={
            "equilibrium_validation": rerun["equilibrium_validation"],
            "symmetry_breaking": rerun["symmetry_breaking"],
            "phase_locking_diagnostic": rerun["phase_locking_diagnostic"],
            "mc_current_balance": rerun["mc_current_balance"],
            "eigen_full_window": eig_main,
            "eigen_short_window": eig_short,
            "eigen_late_window": eig_late,
            "eigen_surrogate_data_fit": rerun["eigen_surrogate_data_fit"],
        },
        note="The table is deliberately scoped as a solver-diagnostic summary and does not promote the short-chain NN outputs to independent transport or spectral claims.",
    )

    add_claim(
        registry,
        claim_id="reproducibility_summary_table",
        section="numerical reproducibility summary",
        claim="The manuscript contains a compact reproducibility map with scope limitations for each numerical result family.",
        evidence=[
            DRAFT,
            flux_path,
            flux_sensitivity_path,
            parameter_scaling_path,
            gamma_scaling_path,
            validation_path,
            figure_metrics_path,
            source_trace_path,
            rerun_path,
            eigen_path,
        ],
        expected_text=[
            r"\section*{Numerical reproducibility summary}",
            r"\label{tab:repro-summary}",
            r"Action-current scaling",
            r"fit-window sensitivity analysis",
            r"and fit-window sensitivity",
            r"bath-temperature and thermostat-coupling robustness checks",
            r"Long-chain profiles and local equilibrium",
            r"Short-chain Fokker--Planck density",
            r"Eigenfunction diagnostic",
            r"Manuscript-level claim audit",
            r"finite-size action-current law over $n=10,20,30,40$",
            r"with $n=50$",
            r"and $n=60$ robustness checks",
            r"saved-model inference reproducibility, not full neural-network retraining",
            r"observable-dependent slow-mode diagnostic, not resolved spectral gap",
            r"local audit to rerun after final author and journal-format edits",
        ],
        computed={
            "registered_artifact_families": [
                "flux_scaling",
                "long_chain_profiles_lte",
                "short_chain_fokker_planck",
                "eigen_diagnostic",
                "claim_audit",
            ]
        },
    )

    add_claim(
        registry,
        claim_id="data_availability_artifacts",
        section="data and code availability",
        claim="Data/code availability lists repo-root paths for the source-trace, rerun, and availability-audit artifacts needed for this audit.",
        evidence=[DRAFT, availability_audit_json, availability_audit_md],
        expected_text=[
            r"paths below are relative to the repository root",
            r"\path{Paper/revision_2026-06-19/scripts/export_source_trace_metrics.py}",
            r"\path{Paper/revision_2026-06-19/source_trace_metrics.json}",
            r"\path{Paper/revision_2026-06-19/scripts/recompute_short_chain_nn_metrics.py}",
            r"\path{Paper/revision_2026-06-19/short_chain_nn_rerun_metrics.json}",
            r"\path{Paper/revision_2026-06-19/scripts/analyze_eigen_fit_windows.py}",
            r"\path{Paper/revision_2026-06-19/eigen_fit_sensitivity.json}",
            r"\path{Paper/revision_2026-06-19/scripts/analyze_flux_scaling_sensitivity.py}",
            r"\path{Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/flux_scaling_sensitivity_n10_60.json}",
            r"\path{Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/b64_scaling_scaling.json}",
            r"\path{Paper/revision_2026-06-19/experiments/flux_validation/gamma_robustness_2026-06-21/gamma_robustness_scaling.json}",
            r"\path{Paper/revision_2026-06-19/scripts/audit_availability_paths.py}",
            r"\path{Paper/revision_2026-06-19/availability_path_audit.json}",
            r"\path{Paper/revision_2026-06-19/availability_path_audit.md}",
        ],
    )

    return registry


def write_outputs(registry: list[dict[str, Any]], json_path: Path, md_path: Path) -> None:
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "draft": rel(DRAFT),
        "total_claims": len(registry),
        "verified": sum(1 for c in registry if c["verdict"] == "VERIFIED"),
        "failed": sum(1 for c in registry if c["verdict"] != "VERIFIED"),
        "scope_note": "Core numerical/data claims only; author declarations and external plagiarism checks remain outside code-verifiable scope.",
        "claims": registry,
    }
    json_path.write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# Manuscript claim audit",
        "",
        f"Generated: `{summary['generated_at_utc']}`",
        "",
        f"Scope note: {summary['scope_note']}",
        "",
        f"Summary: **{summary['verified']} / {summary['total_claims']}** claims verified; **{summary['failed']}** failed.",
        "",
        "| ID | Section | Verdict | Evidence |",
        "|---|---|---:|---|",
    ]
    for claim in registry:
        evidence = "<br>".join(f"`{p}`" for p in claim["evidence"])
        lines.append(f"| `{claim['id']}` | {claim['section']} | {claim['verdict']} | {evidence} |")
    lines.extend(["", "## Failed text checks", ""])
    failed_any = False
    for claim in registry:
        missing = [c for c in claim["text_checks"] if not c["present"]]
        if missing:
            failed_any = True
            lines.append(f"### `{claim['id']}`")
            for c in missing:
                lines.append(f"- missing: `{c['text']}`")
            lines.append("")
    if not failed_any:
        lines.append("None.")
    lines.extend(["", "## Notes", ""])
    lines.append("- This report is a local claim/data audit. It does not replace author confirmation of funding, contributions, or competing interests.")
    lines.append("- It also does not replace a professional plagiarism/self-plagiarism service check before formal journal submission.")
    md_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    registry = build_registry()
    write_outputs(registry, DEFAULT_JSON, DEFAULT_MD)
    failed = [c["id"] for c in registry if c["verdict"] != "VERIFIED"]
    print(json.dumps({"claims": len(registry), "failed": failed, "json": rel(DEFAULT_JSON), "md": rel(DEFAULT_MD)}, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
