#!/usr/bin/env python3
"""Audit consistency of submission-facing metadata across handoff documents.

This catches stale values that are easy to miss during finalization: compiled
PDF SHA-256 values, page counts, release-bundle file counts, and source-bundle
included-file counts.  The authoritative values come from the generated JSON
artifacts used by the one-command submission gate; this script only checks
whether the human-facing handoff documents quote those values consistently.
"""

from __future__ import annotations

import datetime as _dt
import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
DEFAULT_JSON = REVISION / "submission_metadata_consistency_audit.json"
DEFAULT_MD = REVISION / "submission_metadata_consistency_audit.md"

SELF_REFERENTIAL_SOURCE_BUNDLE_OUTPUTS = {
    "Paper/revision_2026-06-19/submission_bundle_manifest.json",
    "Paper/revision_2026-06-19/submission_bundle_manifest.md",
    "Paper/revision_2026-06-19/submission_checks_summary.json",
    "Paper/revision_2026-06-19/submission_checks_summary.md",
    "Paper/revision_2026-06-19/submission_source_bundle_report.json",
    "Paper/revision_2026-06-19/submission_source_bundle_report.md",
}


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def git_files_under(rel_dir: str) -> list[str]:
    prefix = rel_dir.rstrip("/") + "/"
    proc = subprocess.run(
        ["git", "ls-files", prefix],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return [line for line in proc.stdout.splitlines() if line]


def predicted_source_bundle_included_count(manifest: dict[str, Any]) -> int:
    """Mirror the count logic in build_submission_source_bundle.py.

    The source-bundle report itself is generated after this audit in the
    one-command gate, so this function predicts the included-file count from
    the current release manifest rather than relying on a potentially stale
    prior source-bundle report.
    """

    seen: set[str] = set()
    copied_count = 0
    for rec in manifest["release_files"]:
        path = rec["path"]
        if path in SELF_REFERENTIAL_SOURCE_BUNDLE_OUTPUTS:
            continue
        if path in seen:
            continue
        seen.add(path)
        copied_count += 1

    for directory in manifest["release_directories"]:
        for path in git_files_under(directory["path"]):
            if path in SELF_REFERENTIAL_SOURCE_BUNDLE_OUTPUTS:
                continue
            if path in seen:
                continue
            seen.add(path)
            copied_count += 1

    # build_submission_source_bundle.py adds SOURCE_BUNDLE_CONTENTS.sha256.
    return copied_count + 1


def formatted_int(value: int) -> str:
    return f"{value:,}"


def check_contains(checks: list[dict[str, Any]], path: Path, expected: str, label: str) -> None:
    rel_path = rel(path)
    if not path.exists():
        checks.append(
            {
                "path": rel_path,
                "label": label,
                "expected": expected,
                "status": "FAIL",
                "issue": "file missing",
            }
        )
        return
    text = path.read_text(errors="replace")
    status = "PASS" if expected in text else "FAIL"
    checks.append(
        {
            "path": rel_path,
            "label": label,
            "expected": expected,
            "status": status,
            "issue": "" if status == "PASS" else "expected text not found",
        }
    )


def build_checks(values: dict[str, Any]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    generic = values["generic_pdf"]
    siads = values["siads_pdf"]
    release_count = values["release_file_count"]
    source_count = values["predicted_source_bundle_included_count"]

    final_sheet = REVISION / "final_submission_decision_sheet_2026-06-20.md"
    upload_index = REVISION / "journal_upload_file_index_2026-06-20.md"
    siads_packet = REVISION / "siads_first_submission_packet_2026-06-20.md"
    progress = REVISION / "progress_report.md"
    shortlist = REVISION / "target_journal_shortlist_2026-06-19.md"
    final_audit = REVISION / "final_pre_submission_audit_2026-06-20.md"
    integrity = REVISION / "integrity_audit_2026-06-19.md"
    layout = REVISION / "pdf_layout_qa_2026-06-19.md"
    pdf_audit = REVISION / "compiled_pdf_artifact_audit.md"

    check_contains(checks, final_sheet, f"| Release-bundle files | {release_count} |", "final decision release count")
    check_contains(checks, final_sheet, f"| Source-only bundle included files | {source_count} |", "final decision source count")
    check_contains(checks, final_sheet, siads["sha256"], "final decision SIADS PDF SHA")

    check_contains(checks, upload_index, f"| Release-bundle files | {release_count} |", "upload index release count")
    check_contains(checks, upload_index, f"| Source-bundle included files | {source_count} |", "upload index source count")
    check_contains(checks, upload_index, generic["sha256"], "upload index generic PDF SHA")
    check_contains(checks, upload_index, siads["sha256"], "upload index SIADS PDF SHA")

    check_contains(checks, siads_packet, siads["sha256"], "SIADS packet PDF SHA")

    check_contains(checks, progress, f"current build artifact to {generic['pages_from_log']} pages.", "progress report PDF pages")
    check_contains(checks, progress, f"packaging run includes {source_count} regular files", "progress report source count")
    check_contains(checks, progress, siads["sha256"], "progress report SIADS PDF SHA")

    check_contains(checks, shortlist, f"current {generic['pages_from_log']}-page generic article", "shortlist generic PDF pages")
    check_contains(checks, integrity, f"the {generic['pages_from_log']}-page A4 PDF", "integrity audit layout pages")

    check_contains(
        checks,
        final_audit,
        f"`PASS_WITH_LOCAL_RAW_DATA_LIMITATION`, `{release_count}` release files, `0` missing, `0` untracked release files",
        "final audit release count",
    )
    check_contains(checks, final_audit, f"PASS, `{source_count}` included files", "final audit source count")

    check_contains(checks, layout, generic["sha256"], "layout QA generic PDF SHA")
    check_contains(checks, layout, siads["sha256"], "layout QA SIADS PDF SHA")
    check_contains(checks, pdf_audit, generic["sha256"], "compiled PDF audit generic SHA")
    check_contains(checks, pdf_audit, siads["sha256"], "compiled PDF audit SIADS SHA")
    check_contains(checks, pdf_audit, f"| Generic revised manuscript PDF | `{generic['source']}` | {generic['pages_from_log']} |", "compiled PDF audit generic pages")
    check_contains(checks, pdf_audit, f"| SIADS review-preparation PDF | `{siads['source']}` | {siads['pages_from_log']} |", "compiled PDF audit SIADS pages")

    return checks


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    values = payload["values"]
    lines = [
        "# Submission metadata consistency audit",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "## Authoritative values",
        "",
        "| Value | Current value |",
        "|---|---:|",
        f"| Generic PDF pages | {values['generic_pdf']['pages_from_log']} |",
        f"| Generic PDF bytes | {formatted_int(values['generic_pdf']['size_bytes'])} |",
        f"| Generic PDF SHA-256 | `{values['generic_pdf']['sha256']}` |",
        f"| SIADS PDF pages | {values['siads_pdf']['pages_from_log']} |",
        f"| SIADS PDF bytes | {formatted_int(values['siads_pdf']['size_bytes'])} |",
        f"| SIADS PDF SHA-256 | `{values['siads_pdf']['sha256']}` |",
        f"| Release-bundle files | {values['release_file_count']} |",
        f"| Predicted source-bundle included files | {values['predicted_source_bundle_included_count']} |",
        "",
        "## Document checks",
        "",
        "| Document | Check | Status | Expected text |",
        "|---|---|---|---|",
    ]
    for check in payload["checks"]:
        expected = check["expected"].replace("|", "\\|")
        lines.append(f"| `{check['path']}` | {check['label']} | {check['status']} | `{expected}` |")
    failing = [check for check in payload["checks"] if check["status"] != "PASS"]
    if failing:
        lines.extend(["", "## Issues", ""])
        for check in failing:
            lines.append(f"- `{check['path']}` ({check['label']}): {check['issue']}")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    pdf_artifacts = load_json(REVISION / "compiled_pdf_artifact_audit.json")
    manifest = load_json(REVISION / "submission_bundle_manifest.json")

    pdfs = {rec["id"]: rec for rec in pdf_artifacts["pdfs"]}
    values = {
        "generic_pdf": pdfs["generic_manuscript_pdf"],
        "siads_pdf": pdfs["siads_review_pdf"],
        "release_file_count": manifest["summary"]["release_file_count"],
        "predicted_source_bundle_included_count": predicted_source_bundle_included_count(manifest),
    }
    checks = build_checks(values)
    failed = [check for check in checks if check["status"] != "PASS"]
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "status": "PASS" if not failed else "FAIL",
        "summary": {
            "checked_documents": len({check["path"] for check in checks}),
            "total_checks": len(checks),
            "failed_checks": len(failed),
            "release_file_count": values["release_file_count"],
            "predicted_source_bundle_included_count": values["predicted_source_bundle_included_count"],
        },
        "values": values,
        "checks": checks,
    }
    DEFAULT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(DEFAULT_MD, payload)
    print(json.dumps({"status": payload["status"], "summary": payload["summary"]}, indent=2, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
