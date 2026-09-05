#!/usr/bin/env python3
"""Audit manuscript citations and the BibTeX reference list.

This audit is deliberately offline and deterministic.  It checks the local
manuscript/source pair for dangling citations, orphan references, missing
required BibTeX fields, and missing DOI/URL identifiers.  It also records the
external publisher/authority URLs that were manually checked for this revision,
without making the local gate depend on network availability.
"""

from __future__ import annotations

import datetime as _dt
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path


REVISION_REL = Path("Paper/revision_2026-06-19")

REFERENCE_VERIFICATION_SOURCES = {
    "CKSTT": {
        "source_url": "https://link.springer.com/article/10.1007/s00222-010-0242-2",
        "source_type": "Springer Nature article page",
        "checked_on": "2026-06-19",
    },
    "HLNS": {
        "source_url": "https://arxiv.org/abs/2505.16018",
        "source_type": "arXiv abstract page",
        "checked_on": "2026-06-19",
    },
    "ZhaiDobsonLi": {
        "source_url": "https://proceedings.mlr.press/v145/zhai22a.html",
        "source_type": "PMLR proceedings page",
        "checked_on": "2026-06-19",
    },
    "Li2019": {
        "source_url": "https://intlpress.com/JDetail/1806262739393794050",
        "source_type": "International Press journal metadata page",
        "checked_on": "2026-06-19",
    },
    "DobsonLiZhai": {
        "source_url": "https://link.intlpress.com/JDetail/1806261569648545793",
        "source_type": "International Press journal metadata page",
        "checked_on": "2026-06-19",
    },
    "GallavottiCohen": {
        "source_url": "https://link.aps.org/doi/10.1103/PhysRevLett.74.2694",
        "source_type": "APS DOI article page",
        "checked_on": "2026-06-19",
    },
    "LepriLiviPoliti": {
        "source_url": "https://doi.org/10.1016/S0370-1573(02)00558-6",
        "source_type": "DOI resolver for Elsevier Physics Reports article",
        "checked_on": "2026-06-19",
    },
    "Dhar2008": {
        "source_url": "https://arxiv.org/abs/0808.3256",
        "source_type": "arXiv abstract page with Advances in Physics journal reference and related DOI",
        "checked_on": "2026-06-21",
    },
    "Spohn2014": {
        "source_url": "https://arxiv.org/abs/1305.6412",
        "source_type": "arXiv abstract page with Journal of Statistical Physics reference and related DOI",
        "checked_on": "2026-06-21",
    },
    "LebowitzSpohn1999": {
        "source_url": "https://arxiv.org/abs/cond-mat/9811220",
        "source_type": "arXiv abstract page with related Journal of Statistical Physics DOI",
        "checked_on": "2026-06-21",
    },
    "Nazarenko": {
        "source_url": "https://link.springer.com/book/10.1007/978-3-642-15942-8",
        "source_type": "Springer Nature book page",
        "checked_on": "2026-06-19",
    },
}

REQUIRED_FIELDS = {
    "article": {"author", "title", "journal", "year"},
    "book": {"author", "title", "publisher", "year"},
    "inproceedings": {"author", "title", "booktitle", "year"},
    "misc": {"author", "title", "year"},
}


@dataclass
class BibEntry:
    key: str
    entry_type: str
    fields: dict[str, str]


@dataclass
class CitationUse:
    key: str
    line: int
    command: str


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def strip_comments(tex: str) -> str:
    out_lines: list[str] = []
    for line in tex.splitlines():
        cut = len(line)
        for idx, char in enumerate(line):
            if char != "%":
                continue
            backslashes = 0
            j = idx - 1
            while j >= 0 and line[j] == "\\":
                backslashes += 1
                j -= 1
            if backslashes % 2 == 0:
                cut = idx
                break
        out_lines.append(line[:cut])
    return "\n".join(out_lines)


def split_bib_entries(text: str) -> list[str]:
    entries: list[str] = []
    i = 0
    while i < len(text):
        start = text.find("@", i)
        if start == -1:
            break
        brace = text.find("{", start)
        if brace == -1:
            break
        depth = 0
        end = brace
        while end < len(text):
            if text[end] == "{":
                depth += 1
            elif text[end] == "}":
                depth -= 1
                if depth == 0:
                    entries.append(text[start : end + 1])
                    break
            end += 1
        i = end + 1
    return entries


def clean_value(value: str) -> str:
    value = value.strip().rstrip(",").strip()
    if value.startswith("{") and value.endswith("}"):
        value = value[1:-1]
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return " ".join(value.replace("\n", " ").split())


def parse_fields(rest: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    i = 0
    while i < len(rest):
        while i < len(rest) and rest[i] in " \t\r\n,":
            i += 1
        name_start = i
        while i < len(rest) and (rest[i].isalnum() or rest[i] in "_-"):
            i += 1
        if i == name_start:
            break
        name = rest[name_start:i].lower()
        while i < len(rest) and rest[i].isspace():
            i += 1
        if i >= len(rest) or rest[i] != "=":
            break
        i += 1
        while i < len(rest) and rest[i].isspace():
            i += 1
        value_start = i
        depth = 0
        quote = False
        while i < len(rest):
            char = rest[i]
            if char == '"' and (i == 0 or rest[i - 1] != "\\"):
                quote = not quote
            elif not quote:
                if char == "{":
                    depth += 1
                elif char == "}":
                    if depth == 0:
                        break
                    depth -= 1
                elif char == "," and depth == 0:
                    break
            i += 1
        fields[name] = clean_value(rest[value_start:i])
        if i < len(rest) and rest[i] == ",":
            i += 1
    return fields


def parse_bib(path: Path) -> list[BibEntry]:
    entries: list[BibEntry] = []
    for raw in split_bib_entries(path.read_text()):
        match = re.match(r"@([A-Za-z]+)\s*\{\s*([^,]+)\s*,(.*)\}\s*$", raw, re.S)
        if not match:
            continue
        entry_type, key, rest = match.groups()
        entries.append(BibEntry(key=key.strip(), entry_type=entry_type.lower(), fields=parse_fields(rest)))
    return entries


def extract_citations(tex_path: Path) -> list[CitationUse]:
    uses: list[CitationUse] = []
    cite_re = re.compile(r"\\(cite[A-Za-z*]*)\s*(?:\[[^\]]*\]\s*)*\{([^}]+)\}")
    for line_no, raw_line in enumerate(strip_comments(tex_path.read_text()).splitlines(), start=1):
        for match in cite_re.finditer(raw_line):
            command, keys = match.groups()
            for key in keys.split(","):
                clean = key.strip()
                if clean:
                    uses.append(CitationUse(key=clean, line=line_no, command=command))
    return uses


def write_markdown(path: Path, payload: dict) -> None:
    lines = [
        "# Reference integrity audit",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "## Summary",
        "",
        "| Metric | Count |",
        "|---|---:|",
        f"| BibTeX entries | {payload['summary']['bib_entries']} |",
        f"| Citation commands | {payload['summary']['citation_commands']} |",
        f"| Unique cited keys | {payload['summary']['unique_cited_keys']} |",
        f"| Citation uses | {payload['summary']['citation_uses']} |",
        f"| Dangling cited keys | {payload['summary']['dangling_citation_count']} |",
        f"| Orphan BibTeX entries | {payload['summary']['orphan_reference_count']} |",
        f"| Entries missing required fields | {payload['summary']['missing_required_field_count']} |",
        f"| Entries without DOI or URL | {payload['summary']['missing_identifier_count']} |",
        f"| Entries without recorded external source | {payload['summary']['missing_verification_source_count']} |",
        "",
        "## Reference records",
        "",
        "| Key | Type | Cited uses | Identifier | Recorded external source | Status |",
        "|---|---|---:|---|---|---|",
    ]
    for rec in payload["references"]:
        identifier = rec.get("doi") or rec.get("url") or ""
        source = rec["verification_source"]["source_url"] if rec.get("verification_source") else ""
        status = "PASS" if rec["passes"] else "CHECK"
        lines.append(
            f"| `{rec['key']}` | `{rec['entry_type']}` | {rec['citation_count']} | "
            f"`{identifier}` | {source} | {status} |"
        )
    lines.extend(["", "## Citation use by key", "", "| Key | Lines |", "|---|---|"])
    for key, lines_used in payload["citation_lines_by_key"].items():
        lines.append(f"| `{key}` | {', '.join(str(line) for line in lines_used)} |")
    lines.extend(["", "## Issues"])
    issue_sections = [
        ("Dangling cited keys", payload["issues"]["dangling_citations"]),
        ("Orphan BibTeX entries", payload["issues"]["orphan_references"]),
        ("Missing required fields", payload["issues"]["missing_required_fields"]),
        ("Missing DOI/URL identifier", payload["issues"]["missing_identifiers"]),
        ("Missing recorded external source", payload["issues"]["missing_verification_sources"]),
        ("Duplicate BibTeX keys", payload["issues"]["duplicate_bib_keys"]),
    ]
    for title, items in issue_sections:
        lines.extend(["", f"### {title}", ""])
        if not items:
            lines.append("- None.")
        else:
            for item in items:
                lines.append(f"- `{item}`")
    lines.extend(
        [
            "",
            "## Scope limitations",
            "",
            "- This is an offline structural audit of the local TeX/BibTeX pair.",
            "- The recorded external source URLs were checked manually on 2026-06-19; the script does not fetch them during the local gate.",
            "- This audit does not replace a professional plagiarism/self-plagiarism screen or target-journal reference-style conversion.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    root = repo_root()
    revision = root / REVISION_REL
    bib_path = revision / "references.bib"
    tex_path = revision / "draft.tex"

    entries = parse_bib(bib_path)
    uses = extract_citations(tex_path)
    by_key: dict[str, list[CitationUse]] = {}
    for use in uses:
        by_key.setdefault(use.key, []).append(use)

    entry_keys = [entry.key for entry in entries]
    entry_key_set = set(entry_keys)
    cited_key_set = set(by_key)

    duplicate_keys = sorted({key for key in entry_keys if entry_keys.count(key) > 1})
    dangling = sorted(cited_key_set - entry_key_set)
    orphan = sorted(entry_key_set - cited_key_set)

    missing_required: list[str] = []
    missing_identifiers: list[str] = []
    missing_sources: list[str] = []
    reference_records: list[dict] = []

    for entry in entries:
        required = REQUIRED_FIELDS.get(entry.entry_type, {"author", "title", "year"})
        missing = sorted(required - set(entry.fields))
        if missing:
            missing_required.append(f"{entry.key}: {', '.join(missing)}")
        has_identifier = bool(entry.fields.get("doi") or entry.fields.get("url"))
        if not has_identifier:
            missing_identifiers.append(entry.key)
        source = REFERENCE_VERIFICATION_SOURCES.get(entry.key)
        if not source:
            missing_sources.append(entry.key)
        passes = not missing and has_identifier and bool(source) and entry.key not in orphan
        reference_records.append(
            {
                "key": entry.key,
                "entry_type": entry.entry_type,
                "title": entry.fields.get("title", ""),
                "year": entry.fields.get("year", ""),
                "doi": entry.fields.get("doi", ""),
                "url": entry.fields.get("url", ""),
                "citation_count": len(by_key.get(entry.key, [])),
                "citation_lines": [use.line for use in by_key.get(entry.key, [])],
                "missing_required_fields": missing,
                "has_identifier": has_identifier,
                "verification_source": source,
                "passes": passes,
            }
        )

    status = "PASS"
    if any([duplicate_keys, dangling, orphan, missing_required, missing_identifiers, missing_sources]):
        status = "FAIL"

    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision / "scripts/audit_references.py").relative_to(root)),
        "status": status,
        "summary": {
            "bib_entries": len(entries),
            "citation_commands": len(re.findall(r"\\cite[A-Za-z*]*", strip_comments(tex_path.read_text()))),
            "unique_cited_keys": len(cited_key_set),
            "citation_uses": len(uses),
            "dangling_citation_count": len(dangling),
            "orphan_reference_count": len(orphan),
            "missing_required_field_count": len(missing_required),
            "missing_identifier_count": len(missing_identifiers),
            "missing_verification_source_count": len(missing_sources),
        },
        "references": reference_records,
        "citation_lines_by_key": {
            key: [use.line for use in sorted(uses_for_key, key=lambda item: item.line)]
            for key, uses_for_key in sorted(by_key.items())
        },
        "issues": {
            "duplicate_bib_keys": duplicate_keys,
            "dangling_citations": dangling,
            "orphan_references": orphan,
            "missing_required_fields": missing_required,
            "missing_identifiers": missing_identifiers,
            "missing_verification_sources": missing_sources,
        },
    }

    json_path = revision / "reference_integrity_audit.json"
    md_path = revision / "reference_integrity_audit.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps({"status": status, "summary": payload["summary"]}, indent=2, sort_keys=True))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
