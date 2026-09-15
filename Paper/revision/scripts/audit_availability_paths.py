#!/usr/bin/env python3
r"""Audit manuscript data/code-availability and figure paths.

The manuscript mixes two natural path conventions:

* ``\path{...}`` entries in the data/code availability statement should be
  repository-root relative when they begin with a known top-level directory
  such as ``Paper/`` or ``flux/``; otherwise they are interpreted relative to
  the revision directory.
* ``\includegraphics{...}`` entries are LaTeX-source relative, hence relative
  to ``Paper/revision_2026-06-19`` for this draft.

The script writes a machine-readable JSON audit and a compact Markdown summary
next to the manuscript.  It exits non-zero if any referenced path is missing.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path


REVISION_REL = Path("Paper/revision_2026-06-19")
TOP_LEVEL_REPO_DIRS = {
    "Energy Cascade",
    "KDE",
    "Paper",
    "cpp",
    "flux",
    "gibbs_mcmc",
    "lte",
    "matlab",
    "python",
}
SELF_GENERATED_AUDIT_FILES = {
    str(REVISION_REL / "availability_path_audit.json"),
    str(REVISION_REL / "availability_path_audit.md"),
}


@dataclass
class PathRecord:
    kind: str
    manuscript_path: str
    resolved_path: str
    resolution_rule: str
    exists: bool
    path_type: str
    size_bytes: int | None
    sha256: str | None
    git_tracked: bool


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_tracked(root: Path, path: Path) -> bool:
    try:
        rel = path.relative_to(root)
    except ValueError:
        return False
    proc = subprocess.run(
        ["git", "ls-files", "--error-unmatch", str(rel)],
        cwd=root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return proc.returncode == 0


def resolve_manuscript_path(root: Path, revision_dir: Path, raw: str, kind: str) -> tuple[Path, str]:
    raw_path = Path(raw)
    if raw_path.is_absolute():
        return raw_path, "absolute"
    if kind == "figure":
        return revision_dir / raw_path, "latex-source-relative"
    first = raw.split("/", 1)[0]
    if first in TOP_LEVEL_REPO_DIRS:
        return root / raw_path, "repo-root-relative"
    return revision_dir / raw_path, "revision-dir-relative"


def make_record(root: Path, revision_dir: Path, raw: str, kind: str) -> PathRecord:
    resolved, rule = resolve_manuscript_path(root, revision_dir, raw, kind)
    exists = resolved.exists()
    if exists and resolved.is_file():
        path_type = "file"
        size = resolved.stat().st_size
        try:
            rel_to_root = str(resolved.relative_to(root))
        except ValueError:
            rel_to_root = ""
        digest = None if rel_to_root in SELF_GENERATED_AUDIT_FILES else sha256_file(resolved)
    elif exists and resolved.is_dir():
        path_type = "directory"
        size = None
        digest = None
    else:
        path_type = "missing"
        size = None
        digest = None
    try:
        resolved_display = str(resolved.relative_to(root))
    except ValueError:
        resolved_display = str(resolved)
    return PathRecord(
        kind=kind,
        manuscript_path=raw,
        resolved_path=resolved_display,
        resolution_rule=rule,
        exists=exists,
        path_type=path_type,
        size_bytes=size,
        sha256=digest,
        git_tracked=git_tracked(root, resolved) if exists and resolved.is_file() else False,
    )


def extract_paths(tex: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for match in re.finditer(r"\\path\{([^}]+)\}", tex):
        records.append(("availability", match.group(1)))
    for match in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex):
        records.append(("figure", match.group(1)))
    return records


def write_markdown(path: Path, payload: dict) -> None:
    records = payload["records"]
    lines = [
        "# Availability path audit",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Draft: `{payload['draft']}`",
        "",
        "## Summary",
        "",
        f"- Total paths checked: {payload['total_paths']}",
        f"- Missing paths: {payload['missing_paths']}",
        f"- Untracked files among existing file paths: {payload['untracked_files']}",
        "",
        "## Path records",
        "",
        "| Kind | Manuscript path | Resolved path | Rule | Status | Git-tracked | SHA-256 |",
        "|---|---|---|---|---|---:|---|",
    ]
    for rec in records:
        status = "PASS" if rec["exists"] else "MISSING"
        digest = rec["sha256"][:12] if rec["sha256"] else ""
        tracked = "yes" if rec["git_tracked"] else ("n/a" if rec["path_type"] == "directory" else "no")
        lines.append(
            "| {kind} | `{manuscript_path}` | `{resolved_path}` | {resolution_rule} | "
            "{status} ({path_type}) | {tracked} | {digest} |".format(
                **rec,
                status=status,
                tracked=tracked,
                digest=digest,
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This audit checks local path availability only.  It does not prove that",
            "the GitHub branch, journal supplement, or archival DOI contains exactly",
            "the same files unless rerun against that release artifact.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    root = repo_root()
    revision_dir = root / REVISION_REL
    draft = revision_dir / "draft.tex"
    tex = draft.read_text()
    records = [
        asdict(make_record(root, revision_dir, raw, kind))
        for kind, raw in extract_paths(tex)
    ]
    missing = [rec for rec in records if not rec["exists"]]
    untracked_files = [
        rec
        for rec in records
        if rec["exists"] and rec["path_type"] == "file" and not rec["git_tracked"]
    ]
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision_dir / "scripts/audit_availability_paths.py").relative_to(root)),
        "draft": str(draft.relative_to(root)),
        "total_paths": len(records),
        "missing_paths": len(missing),
        "untracked_files": len(untracked_files),
        "records": records,
    }
    json_path = revision_dir / "availability_path_audit.json"
    md_path = revision_dir / "availability_path_audit.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps({k: payload[k] for k in ("total_paths", "missing_paths", "untracked_files")}, indent=2))
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
