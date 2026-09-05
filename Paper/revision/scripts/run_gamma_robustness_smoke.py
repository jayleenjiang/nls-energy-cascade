#!/usr/bin/env python3
"""Build and smoke-test gamma-variant canonical current simulators.

The production current data intentionally keep ``flux/NLS_flux_canonical.cpp``
frozen: its SHA-256 is recorded in the validation manifest and supports the
manuscript's primary action-current claims.  This helper therefore does not
edit that file.  Instead it verifies the frozen source hash, generates a
temporary gamma-specific source under ``tmp/``, compiles it, runs a tiny smoke
simulation, and writes a local report.

The output is a readiness check for a future gamma-robustness production run,
not manuscript evidence for a scaling exponent.
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


REVISION_REL = Path("Paper/revision_2026-06-19")
REPORT_JSON_REL = REVISION_REL / "gamma_robustness_smoke_report.json"
REPORT_MD_REL = REVISION_REL / "gamma_robustness_smoke_report.md"
DEFAULT_RUN_DIR_REL = Path("tmp/gamma_robustness_smoke")
CANONICAL_SOURCE_REL = Path("flux/NLS_flux_canonical.cpp")
EXPECTED_CANONICAL_SHA = "76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gamma_label(gamma: float) -> str:
    text = f"{gamma:g}".replace("-", "m").replace(".", "p")
    return f"gamma{text}"


def command_record(args: list[str]) -> str:
    return " ".join(args)


def run_command(args: list[str], cwd: Path) -> dict:
    proc = subprocess.run(
        args,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return {
        "args": args,
        "returncode": proc.returncode,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def generate_source(root: Path, gamma: float, run_dir: Path) -> Path:
    source = root / CANONICAL_SOURCE_REL
    text = source.read_text()
    label = gamma_label(gamma)
    gamma_literal = format(gamma, ".17g")
    text = text.replace(
        "constexpr double GAMMA = 0.1;",
        f"constexpr double GAMMA = {gamma_literal};",
    )
    text = text.replace(
        'constexpr const char* MODEL_VERSION = "gibbs-canonical-v1";',
        f'constexpr const char* MODEL_VERSION = "gibbs-canonical-v1-{label}-smoke";',
    )
    if "gibbs-canonical-v1-" not in text:
        raise RuntimeError("failed to replace model version in generated source")
    generated = run_dir / label / "NLS_flux_canonical_gamma.cpp"
    generated.parent.mkdir(parents=True, exist_ok=True)
    generated.write_text(text)
    return generated


def compile_source(root: Path, generated_source: Path, binary: Path) -> dict:
    binary.parent.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which("clang++") or shutil.which("g++") or "clang++"
    args = [
        compiler,
        "-O2",
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
    return run_command(args, root)


def run_smoke(root: Path, binary: Path, gamma: float, run_dir: Path) -> dict:
    label = gamma_label(gamma)
    prefix = run_dir / label / "smoke_n6"
    args = [
        str(binary),
        "10",
        "2",
        "6",
        "1",
        "1",
        "1",
        "0.01",
        "20260630",
        "1",
        str(prefix),
    ]
    result = run_command(args, root)
    summary_path = Path(f"{prefix}_summary.csv")
    samples_path = Path(f"{prefix}_samples.csv")
    parsed: dict = {}
    if result["returncode"] == 0 and summary_path.is_file() and samples_path.is_file():
        with summary_path.open(newline="") as f:
            rows = list(csv.DictReader(f))
        if len(rows) == 1:
            row = rows[0]
            parsed = {
                "summary_path": str(summary_path.relative_to(root)),
                "samples_path": str(samples_path.relative_to(root)),
                "model_version": row.get("model_version"),
                "gamma": float(row["gamma"]),
                "n": int(row["n"]),
                "n_trajectories": int(row["n_trajectories"]),
                "mean_action_current": float(row["mean_action_current"]),
                "standard_error": float(row["standard_error"]),
                "projection_count": int(row["projection_count"]),
                "elapsed_seconds": float(row["elapsed_seconds"]),
            }
            parsed["basic_numeric_check"] = all(
                math.isfinite(parsed[key])
                for key in ("gamma", "mean_action_current", "standard_error")
            )
            parsed["gamma_matches_requested"] = math.isclose(
                parsed["gamma"], gamma, rel_tol=0.0, abs_tol=1e-15
            )
    result["parsed_summary"] = parsed
    return result


def write_markdown(path: Path, payload: dict) -> None:
    lines = [
        "# Gamma robustness smoke report",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "## Scope",
        "",
        "This is a code-path smoke test for thermostat-coupling robustness runs.",
        "It is not manuscript evidence and is not used in the current action-current scaling claim.",
        "",
        "The frozen production source is not edited.  The script verifies its SHA-256,",
        "generates gamma-specific temporary sources under `tmp/`, compiles them, and",
        "runs tiny `n=6` simulations only to check the build/run/output path.",
        "",
        "## Frozen-source check",
        "",
        f"- Source: `{payload['canonical_source']}`",
        f"- Expected SHA-256: `{payload['expected_canonical_sha256']}`",
        f"- Observed SHA-256: `{payload['observed_canonical_sha256']}`",
        f"- Match: `{payload['canonical_sha_matches']}`",
        "",
        "## Smoke results",
        "",
        "| gamma | compile | run | summary gamma | mean current | SE | notes |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    for rec in payload["runs"]:
        parsed = rec.get("run", {}).get("parsed_summary", {})
        notes = []
        if not parsed.get("basic_numeric_check"):
            notes.append("numeric check missing/failed")
        if not parsed.get("gamma_matches_requested"):
            notes.append("gamma mismatch")
        if not notes:
            notes.append("smoke ok")
        lines.append(
            f"| {rec['gamma']} | {rec['compile']['returncode']} | "
            f"{rec['run']['returncode']} | {parsed.get('gamma', '---')} | "
            f"{parsed.get('mean_action_current', '---')} | "
            f"{parsed.get('standard_error', '---')} | {'; '.join(notes)} |"
        )
    lines.extend(
        [
            "",
            "## Production-resolution use",
            "",
            "This smoke report remains a build/run readiness check only.  Manuscript",
            "evidence for gamma robustness must come from production-resolution chains",
            "using the same primary lengths `n=10,20,30,40`, the validated timestep,",
            "the existing burn-in schedule, and trajectory-level bootstrap analysis.",
            "Do not fold smoke results into the primary exponent.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gammas", nargs="+", type=float, default=[0.05, 0.2])
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR_REL)
    args = parser.parse_args()

    root = repo_root()
    run_dir = root / args.run_dir
    run_dir.mkdir(parents=True, exist_ok=True)

    canonical = root / CANONICAL_SOURCE_REL
    observed_sha = sha256_file(canonical)
    canonical_match = observed_sha == EXPECTED_CANONICAL_SHA

    runs = []
    status = "PASS"
    if not canonical_match:
        status = "FAIL"

    for gamma in args.gammas:
        generated = generate_source(root, gamma, run_dir)
        binary = run_dir / gamma_label(gamma) / "flux_canonical_gamma"
        compile_result = compile_source(root, generated, binary)
        run_result = {"returncode": None, "parsed_summary": {}}
        if compile_result["returncode"] == 0:
            run_result = run_smoke(root, binary, gamma, run_dir)
        else:
            status = "FAIL"
        parsed = run_result.get("parsed_summary", {})
        if (
            run_result.get("returncode") != 0
            or not parsed.get("basic_numeric_check")
            or not parsed.get("gamma_matches_requested")
        ):
            status = "FAIL"
        runs.append(
            {
                "gamma": gamma,
                "generated_source": str(generated.relative_to(root)),
                "binary": str(binary.relative_to(root)),
                "generated_source_sha256": sha256_file(generated),
                "compile": compile_result,
                "run": run_result,
            }
        )

    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "status": status,
        "scope": "SMOKE_ONLY_NOT_MANUSCRIPT_EVIDENCE",
        "canonical_source": str(CANONICAL_SOURCE_REL),
        "expected_canonical_sha256": EXPECTED_CANONICAL_SHA,
        "observed_canonical_sha256": observed_sha,
        "canonical_sha_matches": canonical_match,
        "run_dir": str(args.run_dir),
        "runs": runs,
    }
    (root / REPORT_JSON_REL).write_text(json.dumps(payload, indent=2) + "\n")
    write_markdown(root / REPORT_MD_REL, payload)
    print(json.dumps({"status": status, "report": str(REPORT_JSON_REL)}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
