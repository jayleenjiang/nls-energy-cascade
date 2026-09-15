#!/usr/bin/env python3
"""Compile and audit the SIADS cover-letter template.

The output is a local template PDF under ``tmp/``.  It contains bracketed
placeholders and an explicit template warning, so it is not a final
author-approved cover letter.  The purpose is to make the SIADS submission
handoff concrete without inventing author declarations.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
SOURCE = REVISION / "siads_cover_letter_template.tex"
BUILD_DIR = ROOT / "tmp/siads_cover_letter_template"
PDF = BUILD_DIR / "siads_cover_letter_template.pdf"
LOG = BUILD_DIR / "siads_cover_letter_template.log"
REPORT_JSON = REVISION / "siads_cover_letter_template_build.json"
REPORT_MD = REVISION / "siads_cover_letter_template_build.md"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scan_log(path: Path) -> list[str]:
    if not path.exists():
        return ["LaTeX log not found."]
    markers = ["LaTeX Warning:", "Overfull \\hbox", "Underfull \\hbox", "Undefined control sequence"]
    issues: list[str] = []
    for line in path.read_text(errors="replace").splitlines():
        if any(marker in line for marker in markers):
            if "rerunfilecheck" in line:
                continue
            issues.append(line.strip())
    return issues


def compile_template() -> dict[str, Any]:
    latexmk = shutil.which("latexmk")
    if latexmk is None:
        return {
            "status": "NOT_RUN",
            "issues": ["latexmk not found on PATH."],
            "source": rel(SOURCE),
            "pdf": rel(PDF),
            "log": rel(LOG),
        }

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        latexmk,
        "-norc",
        "-pdf",
        "-interaction=nonstopmode",
        "-halt-on-error",
        f"-outdir={BUILD_DIR}",
        SOURCE.name,
    ]
    proc = subprocess.run(cmd, cwd=REVISION, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    issues = scan_log(LOG)
    if proc.returncode != 0:
        issues.append(f"latexmk exited with {proc.returncode}")

    pdf_exists = PDF.exists()
    if not pdf_exists:
        issues.append("Compiled cover-letter template PDF not found.")

    return {
        "status": "PASS" if not issues else "FAIL",
        "source": rel(SOURCE),
        "pdf": rel(PDF),
        "log": rel(LOG),
        "pdf_exists": pdf_exists,
        "size_bytes": PDF.stat().st_size if pdf_exists else None,
        "sha256": sha256_file(PDF) if pdf_exists else None,
        "issues": issues,
        "command": cmd,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
        "template_not_final": True,
    }


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# SIADS cover-letter template build",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "This is a template build only.  The compiled PDF contains placeholder",
        "fields and an explicit template warning; it must not be submitted until",
        "authors replace bracketed fields and approve the final cover letter.",
        "",
        "## Artifact",
        "",
        "| Item | Value |",
        "|---|---|",
        f"| Source | `{payload['source']}` |",
        f"| PDF | `{payload['pdf']}` |",
        f"| Size bytes | `{payload.get('size_bytes', '')}` |",
        f"| SHA-256 | `{payload.get('sha256', '')}` |",
        f"| Template-not-final flag | `{payload.get('template_not_final')}` |",
        "",
        "## Issues",
        "",
    ]
    if payload.get("issues"):
        for issue in payload["issues"]:
            lines.append(f"- `{issue}`")
    else:
        lines.append("- None.")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    result = compile_template()
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        **result,
    }
    REPORT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(REPORT_MD, payload)
    print(json.dumps({"status": payload["status"], "pdf": payload["pdf"], "sha256": payload.get("sha256")}, indent=2, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
