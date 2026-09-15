#!/usr/bin/env python3
"""Audit author/journal-only submission fields.

This is deliberately not a hard scientific reproducibility gate.  It detects
provisional manuscript wording and external submission tasks that cannot be
completed from local code/data alone.  The script exits with status 0 so the
local numerical/LaTeX gate can still pass, while the generated audit makes the
remaining author-confirmation blockers explicit.
"""

from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict, dataclass
from pathlib import Path


REVISION_REL = Path("Paper/revision_2026-06-19")


@dataclass
class SubmissionFieldCheck:
    id: str
    category: str
    status: str
    evidence: str
    required_action: str


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def phrase_present(texts: dict[str, str], phrase: str) -> list[str]:
    return [name for name, text in texts.items() if phrase in text]


def make_check(id: str, category: str, pending: bool, evidence: str, required_action: str) -> SubmissionFieldCheck:
    return SubmissionFieldCheck(
        id=id,
        category=category,
        status="PENDING_AUTHOR_OR_EXTERNAL_ACTION" if pending else "PASS",
        evidence=evidence,
        required_action=required_action,
    )


def write_markdown(path: Path, payload: dict) -> None:
    lines = [
        "# Author/submission field audit",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "## Summary",
        "",
        "| Metric | Count |",
        "|---|---:|",
        f"| Checks | {payload['summary']['total_checks']} |",
        f"| Pending author/external items | {payload['summary']['pending_count']} |",
        f"| Passed items | {payload['summary']['pass_count']} |",
        "",
        "## Checks",
        "",
        "| ID | Category | Status | Evidence | Required action |",
        "|---|---|---|---|---|",
    ]
    for rec in payload["checks"]:
        lines.append(
            f"| `{rec['id']}` | {rec['category']} | {rec['status']} | "
            f"{rec['evidence']} | {rec['required_action']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This audit records author, journal, and external-service blockers.  A",
            "`PENDING_AUTHOR_OR_EXTERNAL_ACTION` item is not a local numerical or",
            "LaTeX failure, but the manuscript should not be formally submitted",
            "until the item is resolved and the full local gate is rerun.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    root = repo_root()
    revision = root / REVISION_REL
    source_paths = {
        "draft.tex": revision / "draft.tex",
        "draft_siads_review.tex": revision / "draft_siads_review.tex",
    }
    texts = {name: path.read_text() for name, path in source_paths.items()}
    draft = texts["draft.tex"]

    checks: list[SubmissionFieldCheck] = []

    author_phrase = "Both authors should review and approve the final submitted version."
    hits = phrase_present(texts, author_phrase)
    checks.append(
        make_check(
            "author_approval_placeholder",
            "manuscript-declarations",
            bool(hits),
            f"Provisional approval sentence present in {', '.join(hits)}." if hits else "No provisional author-approval sentence found.",
            "After both authors approve, replace with final target-journal wording.",
        )
    )

    competing_phrase = "No competing interests are declared in the materials supplied for this draft."
    hits = phrase_present(texts, competing_phrase)
    checks.append(
        make_check(
            "competing_interests_placeholder",
            "manuscript-declarations",
            bool(hits),
            f"Provisional competing-interest sentence present in {', '.join(hits)}." if hits else "No provisional competing-interest sentence found.",
            "Confirm with authors and replace with final competing-interest declaration.",
        )
    )

    funding_phrase = "Funding information was not supplied in the current manuscript materials"
    hits = phrase_present(texts, funding_phrase)
    checks.append(
        make_check(
            "funding_placeholder",
            "manuscript-declarations",
            bool(hits),
            f"Provisional funding sentence present in {', '.join(hits)}." if hits else "No provisional funding sentence found.",
            "Confirm funding/no-funding statement with authors and insert final wording.",
        )
    )

    author_block_pending = r"\author{Jayleen Jiang \and Yao Li\thanks" in draft
    checks.append(
        make_check(
            "author_metadata_incomplete",
            "front-matter",
            author_block_pending,
            "Current author block gives Yao Li affiliation/email only." if author_block_pending else "Author block no longer matches the known incomplete pattern.",
            "Confirm final author order, Jayleen affiliation/email, corresponding author, and ORCID choices.",
        )
    )

    repository_only = "github.com/jayleenjiang/nls-energy-cascade" in draft and "zenodo" not in draft.lower() and "doi.org" not in draft.lower()
    checks.append(
        make_check(
            "data_release_route_unfinalized",
            "data-availability",
            repository_only,
            "Data availability currently cites the GitHub repository but no immutable release tag or DOI." if repository_only else "Data availability appears to include a DOI or non-GitHub archive language.",
            "Choose GitHub release-only or DOI-backed raw-data archive route; update data availability and rerun the gate.",
        )
    )

    for item_id, category, evidence, action in [
        (
            "target_journal_confirmation",
            "journal-system",
            "No target-journal choice can be proven from local manuscript text.",
            "Authors must confirm target journal, article type, and whether SIADS review formatting is final.",
        ),
        (
            "professional_similarity_screening",
            "external-service",
            "No iThenticate/Turnitin/journal similarity report is available in the repository.",
            "Run professional similarity/self-plagiarism screening after final edits.",
        ),
        (
            "raw_data_doi_upload",
            "external-service",
            "A local raw-data .tar.gz build exists, but no DOI/upload can be proven locally.",
            "If the DOI route is chosen, upload the archive to Zenodo/OSF/journal storage and insert the DOI.",
        ),
        (
            "final_post_edit_pdf_review",
            "journal-system",
            "Current PDFs are local build artifacts before final author/journal edits.",
            "After final declarations and release metadata are inserted, rerun the gate and visually inspect the final PDF.",
        ),
    ]:
        checks.append(make_check(item_id, category, True, evidence, action))

    pending = [rec for rec in checks if rec.status != "PASS"]
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision / "scripts/audit_author_submission_fields.py").relative_to(root)),
        "status": "PASS" if not pending else "AUTHOR_CONFIRMATION_PENDING",
        "summary": {
            "total_checks": len(checks),
            "pending_count": len(pending),
            "pass_count": len(checks) - len(pending),
        },
        "checks": [asdict(rec) for rec in checks],
    }

    json_path = revision / "author_submission_fields_audit.json"
    md_path = revision / "author_submission_fields_audit.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps({"status": payload["status"], **payload["summary"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
