#!/usr/bin/env python3
"""Apply confirmed author/submission fields to manuscript TeX sources.

Default behavior is a dry run.  The script writes to ``draft.tex`` and
``draft_siads_review.tex`` only when ``--apply`` is supplied and all required
fields are complete.  Before writing, it backs up both TeX sources under
``Paper/revision_2026-06-19/backups/author_submission_fields_<timestamp>/``.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import shutil
from pathlib import Path


REVISION_REL = Path("Paper/revision_2026-06-19")
DEFAULT_INPUT_REL = REVISION_REL / "author_submission_fields.json"
SOURCE_NAMES = ("draft.tex", "draft_siads_review.tex")

CURRENT_AUTHOR_LINE = r"\author{Jayleen Jiang \and Yao Li\thanks{Department of Mathematics and Statistics, University of Massachusetts Amherst, Amherst, MA 01003, USA. \texttt{liyao@umass.edu}}}"
CURRENT_CONTRIBUTIONS = (
    "Jayleen Jiang performed the numerical experiments, assembled the computational\n"
    "artifacts, and drafted the manuscript. Yao Li supervised the project and\n"
    "contributed to the model formulation, theoretical framing, and interpretation.\n"
    "Both authors should review and approve the final submitted version."
)
CURRENT_COMPETING = "No competing interests are declared in the materials supplied for this draft."
CURRENT_FUNDING = (
    "Funding information was not supplied in the current manuscript materials and\n"
    "should be completed by the authors before submission if required by the target\n"
    "journal."
)
CURRENT_DATA_AVAILABILITY = (
    "A public repository snapshot is maintained at\n"
    "\\href{https://github.com/jayleenjiang/nls-energy-cascade}{the project\n"
    "repository on GitHub}."
)

REQUIRED_TEXT_FIELDS = [
    "target_journal",
    "article_type",
    "author_latex",
    "author_contributions_tex",
    "competing_interests_tex",
    "funding_tex",
    "data_release_route",
    "data_availability_tex",
    "repository_release_url",
]
REQUIRED_BOOLEAN_FIELDS = [
    "target_journal_confirmed",
    "similarity_screening_completed",
    "final_pdf_review_completed",
    "author_approval_confirmed",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_fields(path: Path) -> dict:
    return json.loads(path.read_text())


def incomplete_text(value: object) -> bool:
    if not isinstance(value, str):
        return True
    stripped = value.strip()
    if not stripped:
        return True
    lowered = stripped.lower()
    return any(marker in stripped for marker in ("[", "]", "TODO", "TBD")) or "placeholder" in lowered


def validate(fields: dict) -> list[str]:
    problems: list[str] = []
    for key in REQUIRED_TEXT_FIELDS:
        if incomplete_text(fields.get(key)):
            problems.append(f"{key} is missing or still contains placeholder text.")
    route = str(fields.get("data_release_route", "")).strip()
    if route not in {"github_release", "doi_archive"}:
        problems.append("data_release_route must be either 'github_release' or 'doi_archive'.")
    if route == "doi_archive":
        if incomplete_text(fields.get("raw_data_doi_or_accession")):
            problems.append("raw_data_doi_or_accession is required when data_release_route='doi_archive'.")
        if not fields.get("raw_data_upload_completed_if_applicable"):
            problems.append("raw_data_upload_completed_if_applicable must be true when data_release_route='doi_archive'.")
    for key in REQUIRED_BOOLEAN_FIELDS:
        if fields.get(key) is not True:
            problems.append(f"{key} must be true before applying final submission fields.")
    if "should review and approve" in str(fields.get("author_contributions_tex", "")):
        problems.append("author_contributions_tex still contains provisional approval wording.")
    if "not supplied" in str(fields.get("funding_tex", "")).lower():
        problems.append("funding_tex still contains provisional not-supplied wording.")
    if "materials supplied for this draft" in str(fields.get("competing_interests_tex", "")):
        problems.append("competing_interests_tex still contains provisional draft wording.")
    return problems


def replace_exact(text: str, old: str, new: str, label: str, problems: list[str]) -> str:
    if old not in text:
        problems.append(f"Could not find expected text block for {label}.")
        return text
    return text.replace(old, new, 1)


def render_updates(text: str, fields: dict, source_name: str, problems: list[str]) -> str:
    text = replace_exact(text, CURRENT_AUTHOR_LINE, fields["author_latex"], f"{source_name}: author line", problems)
    text = replace_exact(text, CURRENT_CONTRIBUTIONS, fields["author_contributions_tex"], f"{source_name}: author contributions", problems)
    text = replace_exact(text, CURRENT_COMPETING, fields["competing_interests_tex"], f"{source_name}: competing interests", problems)
    text = replace_exact(text, CURRENT_FUNDING, fields["funding_tex"], f"{source_name}: funding", problems)
    text = replace_exact(text, CURRENT_DATA_AVAILABILITY, fields["data_availability_tex"], f"{source_name}: data availability", problems)
    return text


def backup_sources(revision: Path, sources: list[Path]) -> Path:
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup_dir = revision / "backups" / f"author_submission_fields_{stamp}"
    backup_dir.mkdir(parents=True, exist_ok=False)
    for src in sources:
        shutil.copy2(src, backup_dir / src.name)
    return backup_dir


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(DEFAULT_INPUT_REL), help="Path to filled author_submission_fields.json.")
    parser.add_argument("--apply", action="store_true", help="Actually edit draft.tex and draft_siads_review.tex after validation.")
    args = parser.parse_args()

    root = repo_root()
    revision = root / REVISION_REL
    input_path = root / args.input if not Path(args.input).is_absolute() else Path(args.input)
    fields = load_fields(input_path)

    problems = validate(fields)
    sources = [revision / name for name in SOURCE_NAMES]
    updated_texts: dict[Path, str] = {}
    replacement_problems: list[str] = []
    for src in sources:
        updated_texts[src] = render_updates(src.read_text(), fields, src.name, replacement_problems)

    ready = not problems and not replacement_problems
    payload = {
        "input": str(input_path.relative_to(root) if input_path.is_relative_to(root) else input_path),
        "apply_requested": args.apply,
        "ready_to_apply": ready,
        "validation_problems": problems,
        "replacement_problems": replacement_problems,
        "sources": [str(src.relative_to(root)) for src in sources],
    }
    print(json.dumps(payload, indent=2))

    if not args.apply:
        return 0
    if not ready:
        return 1

    backup_dir = backup_sources(revision, sources)
    for src, text in updated_texts.items():
        src.write_text(text)
    print(json.dumps({"applied": True, "backup_dir": str(backup_dir.relative_to(root))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
