#!/usr/bin/env python3
"""Run production-resolution thermostat-coupling robustness experiments.

The manuscript's primary current-scaling claim uses the frozen canonical
implementation ``flux/NLS_flux_canonical.cpp`` with ``gamma=0.1``.  To test
whether the faster-than-Fourier finite-size decay is tied to this single
thermostat coupling, this helper creates gamma-specific temporary sources
without editing the frozen canonical file, runs the same four primary chain
lengths, verifies the summary files against trajectory samples, and fits a
separate power law for each gamma value.

Default production design:

  * T1=10, Tn=2
  * gamma in {0.05, 0.2}
  * n in {10,20,30,40}
  * 64 batches x 16 lanes = 1024 trajectories per length
  * dt=5e-4, measurement window=200
  * burn-ins {1000,1280,2880,5120}

The script is resumable: existing compatible ``*_summary.csv`` and
``*_samples.csv`` files are verified and reused unless ``--force`` is given.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REVISION_REL = Path("Paper/revision_2026-06-19")
CANONICAL_SOURCE_REL = Path("flux/NLS_flux_canonical.cpp")
EXPECTED_CANONICAL_SHA = (
    "76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e"
)
PRIMARY_REFERENCE_SUMMARIES = [
    REVISION_REL
    / "experiments/flux_validation/production_dt5e-4/n10_summary.csv",
    REVISION_REL
    / "experiments/flux_validation/production_dt5e-4/n20_summary.csv",
    REVISION_REL
    / "experiments/flux_validation/production_dt5e-4/n30_summary.csv",
    REVISION_REL
    / "experiments/flux_validation/production_dt5e-4/n40_summary.csv",
]
DEFAULT_OUTPUT_DIR_REL = (
    REVISION_REL
    / "experiments/flux_validation/gamma_robustness_2026-06-21"
)
DEFAULT_BURNINS = {10: 1000.0, 20: 1280.0, 30: 2880.0, 40: 5120.0}
LANES = 16


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gamma_label(gamma: float) -> str:
    return f"gamma{gamma:g}".replace("-", "m").replace(".", "p")


def format_float_for_cpp(value: float) -> str:
    return format(value, ".17g")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gammas", nargs="+", type=float, default=[0.05, 0.2])
    parser.add_argument("--chains", nargs="+", type=int, default=[10, 20, 30, 40])
    parser.add_argument("--T1", type=float, default=10.0)
    parser.add_argument("--Tn", type=float, default=2.0)
    parser.add_argument("--batches", type=int, default=64)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--measure", type=float, default=200.0)
    parser.add_argument("--dt", type=float, default=5.0e-4)
    parser.add_argument("--seed", type=int, default=20260624)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260624)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR_REL)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--no-primary-reference",
        action="store_true",
        help="Do not include the existing gamma=0.1 production data in the comparison plot/report.",
    )
    return parser.parse_args()


def generate_source(root: Path, gamma: float, output_dir: Path) -> Path:
    source = root / CANONICAL_SOURCE_REL
    text = source.read_text()
    label = gamma_label(gamma)
    text = text.replace(
        "constexpr double GAMMA = 0.1;",
        f"constexpr double GAMMA = {format_float_for_cpp(gamma)};",
    )
    text = text.replace(
        'constexpr const char* MODEL_VERSION = "gibbs-canonical-v1";',
        f'constexpr const char* MODEL_VERSION = "gibbs-canonical-v1-{label}";',
    )
    if f"gibbs-canonical-v1-{label}" not in text:
        raise RuntimeError("failed to replace model version")
    generated = output_dir / "generated_sources" / label / "NLS_flux_canonical_gamma.cpp"
    generated.parent.mkdir(parents=True, exist_ok=True)
    generated.write_text(text)
    return generated


def compile_source(root: Path, generated_source: Path, binary: Path, log_path: Path) -> dict:
    binary.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which("clang++") or shutil.which("g++") or "clang++"
    args = [
        compiler,
        "-O3",
        "-mcpu=native",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Xpreprocessor",
        "-fopenmp",
        "-I/opt/homebrew/include/eigen3",
        "-I/opt/homebrew/opt/libomp/include",
        "-L/opt/homebrew/opt/libomp/lib",
        "-lomp",
        str(generated_source),
        "-o",
        str(binary),
    ]
    started = _dt.datetime.now(_dt.timezone.utc).isoformat()
    with log_path.open("w") as log:
        log.write("COMMAND: " + " ".join(args) + "\n\n")
        log.flush()
        proc = subprocess.run(args, cwd=root, stdout=log, stderr=subprocess.STDOUT)
    return {
        "args": args,
        "returncode": proc.returncode,
        "started_at_utc": started,
        "finished_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "log_path": str(log_path.relative_to(root)),
    }


def output_prefix(output_dir: Path, gamma: float, n: int) -> Path:
    return output_dir / gamma_label(gamma) / f"n{n}"


def summary_path_for_prefix(prefix: Path) -> Path:
    return prefix.with_name(prefix.name + "_summary.csv")


def samples_path_for_prefix(prefix: Path) -> Path:
    return prefix.with_name(prefix.name + "_samples.csv")


def compatible_existing_summary(
    path: Path,
    *,
    gamma: float,
    n: int,
    T1: float,
    Tn: float,
    batches: int,
    dt: float,
    burnin: float,
    measure: float,
) -> bool:
    if not path.is_file() or not samples_path_for_prefix(path.with_name(path.name[:-12])).is_file():
        return False
    try:
        row = read_summary(path)
    except Exception:
        return False
    return (
        row["n"] == n
        and math.isclose(row["gamma"], gamma, rel_tol=0.0, abs_tol=1e-14)
        and math.isclose(row["T1"], T1, rel_tol=0.0, abs_tol=1e-14)
        and math.isclose(row["Tn"], Tn, rel_tol=0.0, abs_tol=1e-14)
        and row["batches"] == batches
        and math.isclose(row["dt"], dt, rel_tol=0.0, abs_tol=1e-15)
        and math.isclose(row["burnin"], burnin, rel_tol=0.0, abs_tol=1e-12)
        and math.isclose(row["measure"], measure, rel_tol=0.0, abs_tol=1e-12)
    )


def run_simulation(
    root: Path,
    binary: Path,
    prefix: Path,
    *,
    T1: float,
    Tn: float,
    n: int,
    batches: int,
    burnin: float,
    measure: float,
    dt: float,
    seed: int,
    threads: int,
    log_path: Path,
) -> dict:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    args = [
        str(binary),
        f"{T1:g}",
        f"{Tn:g}",
        str(n),
        str(batches),
        f"{burnin:g}",
        f"{measure:g}",
        f"{dt:.17g}",
        str(seed),
        str(threads),
        str(prefix),
    ]
    print("RUN", " ".join(args), flush=True)
    started = _dt.datetime.now(_dt.timezone.utc).isoformat()
    with log_path.open("w") as log:
        log.write("COMMAND: " + " ".join(args) + "\n\n")
        log.flush()
        proc = subprocess.run(args, cwd=root, stdout=log, stderr=subprocess.STDOUT)
    return {
        "args": args,
        "returncode": proc.returncode,
        "started_at_utc": started,
        "finished_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "log_path": str(log_path.relative_to(root)),
        "summary_path": str(summary_path_for_prefix(prefix).relative_to(root)),
        "samples_path": str(samples_path_for_prefix(prefix).relative_to(root)),
    }


def read_summary(path: Path) -> dict:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one summary row")
    row = dict(rows[0])
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
    if row["lanes"] != LANES:
        raise ValueError(f"{path}: expected {LANES} lanes")
    suffix = "_summary.csv"
    if not path.name.endswith(suffix):
        raise ValueError(f"{path}: expected suffix {suffix}")
    sample_path = path.with_name(path.name[: -len(suffix)] + "_samples.csv")
    raw = np.genfromtxt(sample_path, delimiter=",", names=True)
    samples = np.atleast_1d(raw["time_averaged_action_current"]).astype(float)
    if samples.size != row["n_trajectories"]:
        raise ValueError(
            f"{sample_path}: {samples.size} samples; summary reports "
            f"{row['n_trajectories']}"
        )
    mean = float(np.mean(samples))
    sd = float(np.std(samples, ddof=1))
    se = sd / math.sqrt(samples.size)
    for key, value in [
        ("mean_action_current", mean),
        ("sample_sd", sd),
        ("standard_error", se),
    ]:
        if not math.isclose(row[key], value, rel_tol=2e-12, abs_tol=2e-14):
            raise ValueError(f"{path}: {key} mismatch against raw samples")
    row["stationarity_z"] = (
        row["mean_second_minus_first"] / row["paired_difference_se"]
        if row["paired_difference_se"] > 0
        else float("nan")
    )
    row["summary_path"] = str(path)
    row["sample_path"] = str(sample_path)
    row["samples"] = samples
    return row


def fit_rows(rows: list[dict]) -> tuple[float, float, float]:
    log_n = np.log([row["n"] for row in rows])
    log_mean = np.log([row["mean_action_current"] for row in rows])
    exponent, log_prefactor = np.polyfit(log_n, log_mean, 1)
    prediction = log_prefactor + exponent * log_n
    residual = float(np.sum((log_mean - prediction) ** 2))
    total = float(np.sum((log_mean - np.mean(log_mean)) ** 2))
    r2 = 1.0 - residual / total if total > 0 else float("nan")
    return float(exponent), float(math.exp(log_prefactor)), r2


def bootstrap_ci(rows: list[dict], count: int, seed: int) -> tuple[float, float, int]:
    rng = np.random.default_rng(seed)
    log_n = np.log([row["n"] for row in rows])
    exponents: list[float] = []
    for _ in range(count):
        means = []
        for row in rows:
            samples = row["samples"]
            means.append(float(np.mean(rng.choice(samples, size=samples.size, replace=True))))
        if min(means) > 0.0:
            exponents.append(float(np.polyfit(log_n, np.log(means), 1)[0]))
    if len(exponents) < 0.99 * count:
        raise RuntimeError("too many invalid bootstrap replicates")
    lo, hi = np.quantile(np.asarray(exponents), [0.025, 0.975])
    return float(lo), float(hi), len(exponents)


def make_scaling_record(label: str, rows: list[dict], bootstrap: int, seed: int) -> dict:
    rows = sorted(rows, key=lambda row: row["n"])
    exponent, prefactor, r2 = fit_rows(rows)
    lo, hi, valid = bootstrap_ci(rows, bootstrap, seed)
    return {
        "label": label,
        "gamma": rows[0]["gamma"],
        "chain_lengths": [row["n"] for row in rows],
        "prefactor": prefactor,
        "exponent": exponent,
        "r_squared_log_fit": r2,
        "bootstrap_replicates_requested": bootstrap,
        "bootstrap_replicates_valid": valid,
        "bootstrap_seed": seed,
        "exponent_95_ci": [lo, hi],
        "max_abs_stationarity_z": float(max(abs(row["stationarity_z"]) for row in rows)),
        "rows": [
            {
                key: row[key]
                for key in [
                    "model_version",
                    "n",
                    "T1",
                    "Tn",
                    "gamma",
                    "dt",
                    "burnin",
                    "measure",
                    "batches",
                    "n_trajectories",
                    "seed",
                    "threads",
                    "mean_action_current",
                    "sample_sd",
                    "standard_error",
                    "stationarity_z",
                    "projection_rate",
                    "elapsed_seconds",
                    "summary_path",
                    "sample_path",
                ]
            }
            for row in rows
        ],
    }


def write_summary_csv(path: Path, scaling_records: list[dict]) -> None:
    fieldnames = [
        "label",
        "gamma",
        "n",
        "T1",
        "Tn",
        "batches",
        "trajectories",
        "dt",
        "burnin",
        "measure",
        "mean_action_current",
        "standard_error",
        "stationarity_z",
        "scaling_exponent",
        "scaling_ci_low",
        "scaling_ci_high",
        "scaling_r2",
        "summary_path",
        "sample_path",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for record in scaling_records:
            lo, hi = record["exponent_95_ci"]
            for row in record["rows"]:
                writer.writerow(
                    {
                        "label": record["label"],
                        "gamma": record["gamma"],
                        "n": row["n"],
                        "T1": row["T1"],
                        "Tn": row["Tn"],
                        "batches": row["batches"],
                        "trajectories": row["n_trajectories"],
                        "dt": row["dt"],
                        "burnin": row["burnin"],
                        "measure": row["measure"],
                        "mean_action_current": row["mean_action_current"],
                        "standard_error": row["standard_error"],
                        "stationarity_z": row["stationarity_z"],
                        "scaling_exponent": record["exponent"],
                        "scaling_ci_low": lo,
                        "scaling_ci_high": hi,
                        "scaling_r2": record["r_squared_log_fit"],
                        "summary_path": row["summary_path"],
                        "sample_path": row["sample_path"],
                    }
                )


def write_markdown(path: Path, payload: dict) -> None:
    lines = [
        "# Thermostat-coupling robustness report",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "## Scope",
        "",
        "This report tests whether the faster-than-Fourier finite-size decay of the",
        "action current is tied to the single thermostat coupling used in the primary",
        "production run.  Gamma-specific sources are generated from the frozen",
        "canonical source; the frozen source itself is not edited.",
        "",
        "## Frozen-source check",
        "",
        f"- Source: `{payload['canonical_source']}`",
        f"- Expected SHA-256: `{payload['expected_canonical_sha256']}`",
        f"- Observed SHA-256: `{payload['observed_canonical_sha256']}`",
        f"- Match: `{payload['canonical_sha_matches']}`",
        "",
        "## Scaling results",
        "",
        "| dataset | gamma | n values | exponent | 95% bootstrap CI | R^2 | max |z| |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for record in payload["scaling_records"]:
        lo, hi = record["exponent_95_ci"]
        ns = ",".join(str(n) for n in record["chain_lengths"])
        lines.append(
            f"| {record['label']} | {record['gamma']:.6g} | {ns} | "
            f"{record['exponent']:.5f} | [{lo:.5f}, {hi:.5f}] | "
            f"{record['r_squared_log_fit']:.5f} | {record['max_abs_stationarity_z']:.2f} |"
        )
    lines.extend(
        [
            "",
            "## Per-length means",
            "",
            "| dataset | gamma | n | mean action current | SE | stationarity z |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for record in payload["scaling_records"]:
        for row in record["rows"]:
            lines.append(
                f"| {record['label']} | {record['gamma']:.6g} | {row['n']} | "
                f"{row['mean_action_current']:.10g} | {row['standard_error']:.3g} | "
                f"{row['stationarity_z']:.2f} |"
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "All production-resolution thermostat-coupling datasets in this report are",
            "interpreted only as finite-size robustness checks.  They support the narrower",
            "claim that the observed faster-than-Fourier action-current decay over",
            "`n=10,20,30,40` is not an artifact of the single `gamma=0.1` coupling.",
            "They do not constitute an asymptotic transport theorem or a systematic",
            "two-parameter bath sweep.",
            "",
            "## Output files",
            "",
            f"- JSON: `{payload['json_path']}`",
            f"- CSV: `{payload['summary_csv']}`",
            f"- Plot PDF: `{payload['plot_pdf']}`",
            f"- Plot PNG: `{payload['plot_png']}`",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def make_plot(path_pdf: Path, path_png: Path, scaling_records: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    colors = {
        0.05: "#3b7ddd",
        0.1: "#333333",
        0.2: "#d95f02",
    }
    for record in sorted(scaling_records, key=lambda item: item["gamma"]):
        rows = sorted(record["rows"], key=lambda row: row["n"])
        ns = np.asarray([row["n"] for row in rows], dtype=float)
        means = np.asarray([row["mean_action_current"] for row in rows], dtype=float)
        ses = np.asarray([row["standard_error"] for row in rows], dtype=float)
        gamma = record["gamma"]
        color = colors.get(round(gamma, 8), None)
        label = rf"$\gamma={gamma:g}$, $\alpha={-record['exponent']:.3f}$"
        ax.errorbar(
            ns,
            means,
            yerr=1.959963984540054 * ses,
            fmt="o",
            capsize=3,
            color=color,
            label=label,
        )
        grid = np.linspace(ns.min(), ns.max(), 200)
        ax.plot(
            grid,
            record["prefactor"] * grid ** record["exponent"],
            "--",
            color=color,
            linewidth=1.5,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"chain length $n$")
    ax.set_ylabel(r"mean action current $\langle J\rangle$")
    ax.grid(True, which="both", alpha=0.18)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = repo_root()
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    canonical = root / CANONICAL_SOURCE_REL
    observed_sha = sha256_file(canonical)
    canonical_match = observed_sha == EXPECTED_CANONICAL_SHA
    if not canonical_match:
        raise RuntimeError(
            f"canonical source SHA mismatch: expected {EXPECTED_CANONICAL_SHA}, got {observed_sha}"
        )

    run_records: list[dict] = []
    gamma_rows: dict[float, list[dict]] = {}

    for gamma in args.gammas:
        label = gamma_label(gamma)
        generated_source = generate_source(root, gamma, output_dir)
        binary = output_dir / "bin" / label / "flux_canonical"
        compile_log = output_dir / "logs" / label / "compile.log"
        compile_record = compile_source(root, generated_source, binary, compile_log)
        if compile_record["returncode"] != 0:
            raise RuntimeError(f"compile failed for gamma={gamma}; see {compile_log}")
        print(f"COMPILED gamma={gamma:g}: {binary}", flush=True)

        rows: list[dict] = []
        for n in args.chains:
            if n not in DEFAULT_BURNINS:
                raise ValueError(f"no default burn-in for n={n}; update DEFAULT_BURNINS")
            burnin = DEFAULT_BURNINS[n]
            prefix = output_prefix(output_dir, gamma, n)
            summary_path = summary_path_for_prefix(prefix)
            samples_path = samples_path_for_prefix(prefix)
            log_path = output_dir / "logs" / label / f"n{n}.log"
            if (
                not args.force
                and summary_path.is_file()
                and samples_path.is_file()
                and compatible_existing_summary(
                    summary_path,
                    gamma=gamma,
                    n=n,
                    T1=args.T1,
                    Tn=args.Tn,
                    batches=args.batches,
                    dt=args.dt,
                    burnin=burnin,
                    measure=args.measure,
                )
            ):
                run_record = {
                    "args": [],
                    "returncode": 0,
                    "started_at_utc": None,
                    "finished_at_utc": None,
                    "log_path": str(log_path.relative_to(root)),
                    "summary_path": str(summary_path.relative_to(root)),
                    "samples_path": str(samples_path.relative_to(root)),
                    "reused_existing": True,
                }
                print(f"REUSE gamma={gamma:g} n={n}: {summary_path}", flush=True)
            else:
                run_record = run_simulation(
                    root,
                    binary,
                    prefix,
                    T1=args.T1,
                    Tn=args.Tn,
                    n=n,
                    batches=args.batches,
                    burnin=burnin,
                    measure=args.measure,
                    dt=args.dt,
                    seed=args.seed,
                    threads=args.threads,
                    log_path=log_path,
                )
                run_record["reused_existing"] = False
                if run_record["returncode"] != 0:
                    raise RuntimeError(f"run failed for gamma={gamma} n={n}; see {log_path}")
            row = read_summary(summary_path)
            rows.append(row)
            run_records.append(
                {
                    "gamma": gamma,
                    "n": n,
                    "generated_source": str(generated_source.relative_to(root)),
                    "generated_source_sha256": sha256_file(generated_source),
                    "binary": str(binary.relative_to(root)),
                    "binary_sha256": sha256_file(binary),
                    "compile": compile_record,
                    "run": run_record,
                }
            )
        gamma_rows[gamma] = rows

    scaling_records: list[dict] = []
    if not args.no_primary_reference:
        reference_rows = [read_summary(root / path) for path in PRIMARY_REFERENCE_SUMMARIES]
        scaling_records.append(
            make_scaling_record(
                "primary_reference_gamma0p1",
                reference_rows,
                args.bootstrap,
                args.bootstrap_seed + 10,
            )
        )
    for gamma, rows in sorted(gamma_rows.items()):
        scaling_records.append(
            make_scaling_record(
                gamma_label(gamma),
                rows,
                args.bootstrap,
                args.bootstrap_seed + int(round(gamma * 1_000_000)),
            )
        )

    status = "PASS"
    failed_reasons: list[str] = []
    for record in scaling_records:
        if record["exponent"] >= -1.0:
            status = "FAIL"
            failed_reasons.append(
                f"{record['label']}: exponent {record['exponent']:.5f} is not faster than Fourier"
            )
        if record["max_abs_stationarity_z"] >= 2.0:
            if status != "FAIL":
                status = "WARN"
            failed_reasons.append(
                f"{record['label']}: max |stationarity z|={record['max_abs_stationarity_z']:.2f}"
            )

    plot_pdf = output_dir / "gamma_robustness_scaling.pdf"
    plot_png = output_dir / "gamma_robustness_scaling.png"
    summary_csv = output_dir / "gamma_robustness_summary.csv"
    json_path = output_dir / "gamma_robustness_scaling.json"
    report_md = output_dir / "gamma_robustness_report.md"

    write_summary_csv(summary_csv, scaling_records)
    make_plot(plot_pdf, plot_png, scaling_records)

    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "status": status,
        "failed_reasons": failed_reasons,
        "scope": "THERMOSTAT_COUPLING_FINITE_SIZE_ROBUSTNESS",
        "canonical_source": str(CANONICAL_SOURCE_REL),
        "expected_canonical_sha256": EXPECTED_CANONICAL_SHA,
        "observed_canonical_sha256": observed_sha,
        "canonical_sha_matches": canonical_match,
        "output_dir": str(args.output_dir),
        "json_path": str(json_path.relative_to(root)),
        "summary_csv": str(summary_csv.relative_to(root)),
        "plot_pdf": str(plot_pdf.relative_to(root)),
        "plot_png": str(plot_png.relative_to(root)),
        "run_records": run_records,
        "scaling_records": scaling_records,
    }
    json_safe = json.loads(json.dumps(payload, default=lambda _: "<non-json>"))
    json_path.write_text(json.dumps(json_safe, indent=2) + "\n")
    payload_for_md = dict(json_safe)
    write_markdown(report_md, payload_for_md)
    print(json.dumps({"status": status, "report": str(report_md.relative_to(root))}, indent=2))
    return 0 if status in {"PASS", "WARN"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
