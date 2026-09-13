#!/usr/bin/env python3
"""Audit local compiled manuscript PDF artifacts.

The manuscript PDFs are build artifacts under ``tmp/`` and are intentionally
not committed.  This script records enough metadata to make the current local
build auditable: path, source, log, page count reported by TeX, file size, and
SHA-256.  It also scans the corresponding LaTeX logs for the same warning
markers used by the submission gate.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
DEFAULT_JSON = REVISION / "compiled_pdf_artifact_audit.json"
DEFAULT_MD = REVISION / "compiled_pdf_artifact_audit.md"

PDF_TARGETS = [
    {
        "id": "generic_manuscript_pdf",
        "source": REVISION / "draft.tex",
        "pdf": ROOT / "tmp/paper_build/revision/draft.pdf",
        "log": ROOT / "tmp/paper_build/revision/draft.log",
        "role": "Generic revised manuscript PDF",
    },
    {
        "id": "siads_review_pdf",
        "source": REVISION / "draft_siads_review.tex",
        "pdf": ROOT / "tmp/paper_build/siads_review/draft_siads_review.pdf",
        "log": ROOT / "tmp/paper_build/siads_review/draft_siads_review.log",
        "role": "SIADS review-preparation PDF",
    },
]

ISSUE_MARKERS = [
    "Overfull \\hbox",
    "Underfull \\hbox",
    "LaTeX Warning:",
    "Package natbib Warning:",
    "Citation",
    "undefined references",
    "Undefined control sequence",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scan_log(log_path: Path) -> tuple[list[str], int | None, int | None]:
    if not log_path.exists():
        return ["LaTeX log not found."], None, None

    log_text = log_path.read_text(errors="replace")
    issues: list[str] = []
    for line in log_text.splitlines():
        if any(marker in line for marker in ISSUE_MARKERS):
            if "rerunfilecheck" in line:
                continue
            issues.append(line.strip())

    # TeX may wrap the output path across display lines, so parse the stable
    # parenthesized summary rather than the path itself.
    matches = re.findall(r"\((\d+)\s+pages,\s+(\d+)\s+bytes\)", log_text)
    pages: int | None = None
    bytes_from_log: int | None = None
    if matches:
        pages = int(matches[-1][0])
        bytes_from_log = int(matches[-1][1])
    else:
        issues.append("Could not find TeX output page/byte summary in log.")
    return issues, pages, bytes_from_log


def audit_target(target: dict[str, Path | str]) -> dict[str, Any]:
    source = target["source"]
    pdf = target["pdf"]
    log = target["log"]
    assert isinstance(source, Path)
    assert isinstance(pdf, Path)
    assert isinstance(log, Path)

    issues, pages, bytes_from_log = scan_log(log)
    pdf_exists = pdf.exists()
    source_exists = source.exists()
    log_exists = log.exists()
    size_bytes = pdf.stat().st_size if pdf_exists else None
    sha256 = sha256_file(pdf) if pdf_exists else None
    if not source_exists:
        issues.append("TeX source not found.")
    if not pdf_exists:
        issues.append("Compiled PDF not found.")
    if pdf_exists and bytes_from_log is not None and size_bytes != bytes_from_log:
        issues.append(f"PDF size {size_bytes} does not match TeX log byte count {bytes_from_log}.")

    return {
        "id": target["id"],
        "role": target["role"],
        "source": rel(source),
        "pdf": rel(pdf),
        "log": rel(log),
        "source_exists": source_exists,
        "pdf_exists": pdf_exists,
        "log_exists": log_exists,
        "pages_from_log": pages,
        "size_bytes": size_bytes,
        "bytes_from_log": bytes_from_log,
        "bytes_match_log": pdf_exists and bytes_from_log is not None and size_bytes == bytes_from_log,
        "sha256": sha256,
        "issues": issues,
        "status": "PASS" if not issues else "FAIL",
    }


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# Compiled PDF artifact audit",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "The PDFs listed here are local build artifacts under `tmp/` and are not committed to Git.",
        "Regenerate them after final author declarations, journal-template changes, or release metadata edits.",
        "",
        "| PDF | Source | Pages | Size bytes | SHA-256 | Status |",
        "|---|---|---:|---:|---|---|",
    ]
    for rec in payload["pdfs"]:
        sha = rec["sha256"] or ""
        lines.append(
            f"| {rec['role']} | `{rec['source']}` | {rec['pages_from_log'] or ''} | "
            f"{rec['size_bytes'] or ''} | `{sha}` | {rec['status']} |"
        )
    failing = [rec for rec in payload["pdfs"] if rec["issues"]]
    if failing:
        lines.extend(["", "## Issues", ""])
        for rec in failing:
            lines.append(f"### {rec['id']}")
            for issue in rec["issues"]:
                lines.append(f"- `{issue}`")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    records = [audit_target(target) for target in PDF_TARGETS]
    status = "PASS" if all(rec["status"] == "PASS" for rec in records) else "FAIL"
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "status": status,
        "pdfs": records,
    }
    DEFAULT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(DEFAULT_MD, payload)
    print(json.dumps({"status": status, "pdfs": len(records), "failed": [r["id"] for r in records if r["status"] != "PASS"]}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
