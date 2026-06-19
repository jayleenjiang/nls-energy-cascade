#!/usr/bin/env python3
"""Run the local submission-readiness checks for the revised manuscript.

This runner standardizes the order of the local gates that can be checked
without author-only information or external services:

1. optional LaTeX compile + log scan;
2. data/code and figure path audit;
3. manuscript numerical-claim audit;
4. submission bundle manifest;
5. minimal raw-data archive manifest;
6. final submission bundle manifest refresh.

The script intentionally does not run professional plagiarism checking, upload
raw data, select a journal template, or fill author/funding declarations.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path


REVISION_REL = Path("Paper/revision_2026-06-19")
BUILD_REL = Path("tmp/paper_build/revision")


@dataclass
class CommandResult:
    name: str
    command: list[str]
    returncode: int
    duration_seconds: float
    stdout_tail: str
    stderr_tail: str


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def run_command(root: Path, name: str, command: list[str], *, cwd: Path | None = None) -> CommandResult:
    t0 = time.time()
    proc = subprocess.run(
        command,
        cwd=cwd or root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return CommandResult(
        name=name,
        command=command,
        returncode=proc.returncode,
        duration_seconds=round(time.time() - t0, 3),
        stdout_tail=proc.stdout[-4000:],
        stderr_tail=proc.stderr[-4000:],
    )


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def latex_log_check(log_path: Path) -> dict:
    if not log_path.exists():
        return {"status": "NOT_RUN", "log": str(log_path), "issues": ["LaTeX log not found."]}
    issue_markers = [
        "Overfull \\hbox",
        "Underfull \\hbox",
        "LaTeX Warning:",
        "Package natbib Warning:",
        "Citation",
        "undefined references",
        "Undefined control sequence",
    ]
    issues: list[str] = []
    for line in log_path.read_text(errors="replace").splitlines():
        if any(marker in line for marker in issue_markers):
            # Avoid false positive from package names or explanatory text.
            if "rerunfilecheck" in line:
                continue
            issues.append(line.strip())
    return {"status": "PASS" if not issues else "FAIL", "log": str(log_path), "issues": issues}


def compile_latex(root: Path) -> tuple[CommandResult | None, dict]:
    latexmk = shutil.which("latexmk")
    log_path = root / BUILD_REL / "draft.log"
    if latexmk is None:
        return None, {
            "status": "NOT_RUN",
            "log": str(log_path.relative_to(root)),
            "issues": ["latexmk not found on PATH."],
        }
    build_dir = root / BUILD_REL
    revision_dir = root / REVISION_REL
    build_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        latexmk,
        "-norc",
        "-pdf",
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-synctex=1",
        f"-outdir={build_dir}",
        "draft.tex",
    ]
    result = run_command(root, "latex_compile", cmd, cwd=revision_dir)
    check = latex_log_check(log_path)
    if result.returncode != 0:
        check["status"] = "FAIL"
        check.setdefault("issues", []).append(f"latexmk exited with {result.returncode}")
    return result, check


def write_markdown(path: Path, payload: dict) -> None:
    lines = [
        "# Submission checks summary",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Overall status: **{payload['overall_status']}**",
        "",
        "## Gate summary",
        "",
        "| Gate | Status | Key numbers |",
        "|---|---|---|",
    ]
    for gate in payload["gates"]:
        numbers = ", ".join(f"{k}={v}" for k, v in gate.get("numbers", {}).items())
        lines.append(f"| {gate['name']} | {gate['status']} | {numbers} |")
    lines.extend(["", "## Command results", "", "| Command | Return code | Duration (s) |", "|---|---:|---:|"])
    for result in payload["commands"]:
        lines.append(f"| `{result['name']}` | {result['returncode']} | {result['duration_seconds']} |")
    if payload.get("latex_log"):
        lines.extend(["", "## LaTeX log check", "", f"- Status: {payload['latex_log']['status']}"])
        if payload["latex_log"].get("issues"):
            lines.append("- Issues:")
            for issue in payload["latex_log"]["issues"][:20]:
                lines.append(f"  - `{issue}`")
    lines.extend(
        [
            "",
            "## Scope limitations",
            "",
            "- This runner checks local reproducibility gates only.",
            "- It does not confirm author/funding/competing-interest declarations.",
            "- It does not run professional plagiarism/self-plagiarism screening.",
            "- It does not upload raw data or create a DOI-backed archive.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compile-latex", action="store_true", help="Run latexmk before scanning the LaTeX log.")
    parser.add_argument(
        "--require-latex",
        action="store_true",
        help="Fail if LaTeX was not compiled or the existing LaTeX log is missing.",
    )
    args = parser.parse_args()

    root = repo_root()
    revision = root / REVISION_REL
    py = "python3"
    commands: list[CommandResult] = []

    latex_result: CommandResult | None = None
    if args.compile_latex:
        latex_result, latex_log = compile_latex(root)
        if latex_result is not None:
            commands.append(latex_result)
    else:
        latex_log = latex_log_check(root / BUILD_REL / "draft.log")
        if latex_log["status"] == "PASS":
            latex_log["status"] = "PASS_EXISTING_LOG"
        elif not args.require_latex:
            latex_log["status"] = "NOT_REQUIRED"

    for name, script in [
        ("availability_path_audit", "audit_availability_paths.py"),
        ("manuscript_claim_audit", "audit_manuscript_claims.py"),
        ("submission_bundle_manifest_initial", "build_submission_bundle_manifest.py"),
        ("raw_data_archive_manifest", "build_raw_data_archive_manifest.py"),
        ("submission_bundle_manifest_final", "build_submission_bundle_manifest.py"),
    ]:
        result = run_command(root, name, [py, str(revision / "scripts" / script)])
        commands.append(result)

    availability = load_json(revision / "availability_path_audit.json")
    claims = load_json(revision / "manuscript_claim_audit.json")
    raw = load_json(revision / "raw_data_archive_manifest.json")
    bundle = load_json(revision / "submission_bundle_manifest.json")

    gates = [
        {
            "name": "latex_log",
            "status": latex_log["status"],
            "numbers": {"issues": len(latex_log.get("issues", []))},
        },
        {
            "name": "availability_path_audit",
            "status": "PASS" if availability["missing_paths"] == 0 and availability["untracked_files"] == 0 else "FAIL",
            "numbers": {
                "total_paths": availability["total_paths"],
                "missing_paths": availability["missing_paths"],
                "untracked_files": availability["untracked_files"],
            },
        },
        {
            "name": "manuscript_claim_audit",
            "status": "PASS" if claims["failed"] == 0 else "FAIL",
            "numbers": {"total_claims": claims["total_claims"], "verified": claims["verified"], "failed": claims["failed"]},
        },
        {
            "name": "raw_data_archive_manifest",
            "status": raw["summary"]["status"],
            "numbers": raw["summary"],
        },
        {
            "name": "submission_bundle_manifest",
            "status": bundle["summary"]["status"],
            "numbers": bundle["summary"],
        },
    ]

    hard_fail = any(cmd.returncode != 0 for cmd in commands)
    hard_fail = hard_fail or any(gate["status"] == "FAIL" for gate in gates)
    if args.require_latex and gates[0]["status"] not in {"PASS", "PASS_EXISTING_LOG"}:
        hard_fail = True
    overall = "PASS" if not hard_fail else "FAIL"
    if overall == "PASS" and bundle["summary"]["status"] == "PASS_WITH_LOCAL_RAW_DATA_LIMITATION":
        overall = "PASS_WITH_LOCAL_RAW_DATA_LIMITATION"
    if overall == "PASS" and gates[0]["status"] in {"NOT_REQUIRED", "NOT_RUN"}:
        overall = "PASS_WITHOUT_LATEX"

    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision / "scripts/run_submission_checks.py").relative_to(root)),
        "overall_status": overall,
        "commands": [asdict(c) for c in commands],
        "gates": gates,
        "latex_log": latex_log,
    }
    json_path = revision / "submission_checks_summary.json"
    md_path = revision / "submission_checks_summary.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps({"overall_status": overall, "gates": gates}, indent=2))
    return 0 if overall.startswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
