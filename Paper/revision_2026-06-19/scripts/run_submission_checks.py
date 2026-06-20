#!/usr/bin/env python3
"""Run the local submission-readiness checks for the revised manuscript.

This runner standardizes the order of the local gates that can be checked
without author-only information or external services:

1. optional LaTeX compile + log scan for the generic and SIADS review sources;
2. compiled PDF artifact audit;
3. data/code and figure path audit;
4. manuscript numerical-claim audit;
5. citation/reference integrity audit;
6. author/journal submission-field audit;
7. submission bundle manifest;
8. minimal raw-data archive manifest;
9. final submission bundle manifest refresh;
10. source-only submission bundle packaging dry run.

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
SIADS_BUILD_REL = Path("tmp/paper_build/siads_review")


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


def compile_latex_source(
    root: Path,
    *,
    source_name: str,
    build_rel: Path,
    log_name: str,
    command_name: str,
) -> tuple[CommandResult | None, dict]:
    latexmk = shutil.which("latexmk")
    log_path = root / build_rel / log_name
    if latexmk is None:
        return None, {
            "status": "NOT_RUN",
            "log": str(log_path.relative_to(root)),
            "issues": ["latexmk not found on PATH."],
        }
    build_dir = root / build_rel
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
        source_name,
    ]
    result = run_command(root, command_name, cmd, cwd=revision_dir)
    check = latex_log_check(log_path)
    if result.returncode != 0:
        check["status"] = "FAIL"
        check.setdefault("issues", []).append(f"latexmk exited with {result.returncode}")
    return result, check


def compile_latex(root: Path) -> tuple[CommandResult | None, dict]:
    return compile_latex_source(
        root,
        source_name="draft.tex",
        build_rel=BUILD_REL,
        log_name="draft.log",
        command_name="latex_compile",
    )


def compile_siads_latex(root: Path) -> tuple[CommandResult | None, dict]:
    return compile_latex_source(
        root,
        source_name="draft_siads_review.tex",
        build_rel=SIADS_BUILD_REL,
        log_name="draft_siads_review.log",
        command_name="siads_latex_compile",
    )


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
    latex_checks = [
        ("Generic manuscript", payload.get("latex_log")),
        ("SIADS review source", payload.get("siads_latex_log")),
    ]
    if any(check for _, check in latex_checks):
        lines.extend(["", "## LaTeX log checks", ""])
        for label, check in latex_checks:
            if not check:
                continue
            lines.append(f"- {label}: {check['status']}")
            if check.get("issues"):
                lines.append(f"  - Issues in {label}:")
                for issue in check["issues"][:20]:
                    lines.append(f"    - `{issue}`")
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
    siads_latex_result: CommandResult | None = None
    if args.compile_latex:
        latex_result, latex_log = compile_latex(root)
        if latex_result is not None:
            commands.append(latex_result)
        siads_latex_result, siads_latex_log = compile_siads_latex(root)
        if siads_latex_result is not None:
            commands.append(siads_latex_result)
    else:
        latex_log = latex_log_check(root / BUILD_REL / "draft.log")
        if latex_log["status"] == "PASS":
            latex_log["status"] = "PASS_EXISTING_LOG"
        elif not args.require_latex:
            latex_log["status"] = "NOT_REQUIRED"
        siads_latex_log = latex_log_check(root / SIADS_BUILD_REL / "draft_siads_review.log")
        if siads_latex_log["status"] == "PASS":
            siads_latex_log["status"] = "PASS_EXISTING_LOG"
        elif not args.require_latex:
            siads_latex_log["status"] = "NOT_REQUIRED"

    for name, script in [
        ("compiled_pdf_artifact_audit", "audit_compiled_pdfs.py"),
        ("availability_path_audit", "audit_availability_paths.py"),
        ("manuscript_claim_audit", "audit_manuscript_claims.py"),
        ("reference_integrity_audit", "audit_references.py"),
        ("author_submission_fields_audit", "audit_author_submission_fields.py"),
        ("submission_bundle_manifest_initial", "build_submission_bundle_manifest.py"),
        ("raw_data_archive_manifest", "build_raw_data_archive_manifest.py"),
        ("submission_bundle_manifest_final", "build_submission_bundle_manifest.py"),
        ("submission_source_bundle", "build_submission_source_bundle.py"),
    ]:
        result = run_command(root, name, [py, str(revision / "scripts" / script)])
        commands.append(result)

    pdf_artifacts = load_json(revision / "compiled_pdf_artifact_audit.json")
    availability = load_json(revision / "availability_path_audit.json")
    claims = load_json(revision / "manuscript_claim_audit.json")
    references = load_json(revision / "reference_integrity_audit.json")
    author_fields = load_json(revision / "author_submission_fields_audit.json")
    raw = load_json(revision / "raw_data_archive_manifest.json")
    bundle = load_json(revision / "submission_bundle_manifest.json")
    source_bundle = load_json(revision / "submission_source_bundle_report.json")

    gates = [
        {
            "name": "latex_log",
            "status": latex_log["status"],
            "numbers": {"issues": len(latex_log.get("issues", []))},
        },
        {
            "name": "siads_latex_log",
            "status": siads_latex_log["status"],
            "numbers": {"issues": len(siads_latex_log.get("issues", []))},
        },
        {
            "name": "compiled_pdf_artifact_audit",
            "status": pdf_artifacts["status"],
            "numbers": {
                "pdf_count": len(pdf_artifacts["pdfs"]),
                "failed": sum(1 for rec in pdf_artifacts["pdfs"] if rec["status"] != "PASS"),
            },
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
            "name": "reference_integrity_audit",
            "status": references["status"],
            "numbers": references["summary"],
        },
        {
            "name": "author_submission_fields_audit",
            "status": author_fields["status"],
            "numbers": author_fields["summary"],
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
        {
            "name": "submission_source_bundle",
            "status": source_bundle["status"],
            "numbers": source_bundle["summary"],
        },
    ]

    hard_fail = any(cmd.returncode != 0 for cmd in commands)
    hard_fail = hard_fail or any(gate["status"] == "FAIL" for gate in gates)
    latex_gate_statuses = {
        "latex_log": latex_log["status"],
        "siads_latex_log": siads_latex_log["status"],
    }
    if args.require_latex and any(status not in {"PASS", "PASS_EXISTING_LOG"} for status in latex_gate_statuses.values()):
        hard_fail = True
    overall = "PASS" if not hard_fail else "FAIL"
    if overall == "PASS" and bundle["summary"]["status"] == "PASS_WITH_LOCAL_RAW_DATA_LIMITATION":
        overall = "PASS_WITH_LOCAL_RAW_DATA_LIMITATION"
    if author_fields["status"] == "AUTHOR_CONFIRMATION_PENDING" and overall.startswith("PASS"):
        if overall == "PASS_WITH_LOCAL_RAW_DATA_LIMITATION":
            overall = "PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION"
        else:
            overall = "PASS_WITH_AUTHOR_CONFIRMATION_PENDING"
    if overall == "PASS" and all(status in {"NOT_REQUIRED", "NOT_RUN"} for status in latex_gate_statuses.values()):
        overall = "PASS_WITHOUT_LATEX"

    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision / "scripts/run_submission_checks.py").relative_to(root)),
        "overall_status": overall,
        "commands": [asdict(c) for c in commands],
        "gates": gates,
        "latex_log": latex_log,
        "siads_latex_log": siads_latex_log,
    }
    json_path = revision / "submission_checks_summary.json"
    md_path = revision / "submission_checks_summary.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps({"overall_status": overall, "gates": gates}, indent=2))
    return 0 if overall.startswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
