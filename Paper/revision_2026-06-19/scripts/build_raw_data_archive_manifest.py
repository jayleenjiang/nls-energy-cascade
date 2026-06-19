#!/usr/bin/env python3
"""Build a minimal raw-data archive manifest for the revised paper.

This script reads ``submission_bundle_manifest.json`` and extracts only the
local raw-data files that are actually referenced by the source-trace artifacts
used in the manuscript.  It does not copy, compress, upload, or delete data.

The resulting manifest is intended to support a future Zenodo/OSF/journal
supplement decision: it identifies a compact raw-data subset (currently much
smaller than the full local raw-data roots) and records SHA-256 checksums for
each file.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


REVISION_REL = Path("Paper/revision_2026-06-19")


@dataclass
class RawArchiveRecord:
    archive_path: str
    local_path: str
    source_artifacts: list[str]
    exists: bool
    size_bytes: int | None
    sha256: str | None


@dataclass
class RawRootSummary:
    root: str
    referenced_file_count: int
    referenced_total_bytes: int
    full_local_root_file_count: int | None
    full_local_root_total_bytes: int | None


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_git(root: Path, args: list[str]) -> str:
    proc = subprocess.run(["git", *args], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout.strip()


def archive_path_for(local_path: str) -> str:
    return f"raw_data/{local_path}"


def group_records(records: Iterable[dict]) -> list[RawArchiveRecord]:
    grouped: dict[str, set[str]] = {}
    for rec in records:
        grouped.setdefault(rec["path"], set()).add(rec["source_artifact"])
    root = repo_root()
    out: list[RawArchiveRecord] = []
    for local_path, sources in sorted(grouped.items()):
        path = root / local_path
        exists = path.is_file()
        out.append(
            RawArchiveRecord(
                archive_path=archive_path_for(local_path),
                local_path=local_path,
                source_artifacts=sorted(sources),
                exists=exists,
                size_bytes=path.stat().st_size if exists else None,
                sha256=sha256_file(path) if exists else None,
            )
        )
    return out


def summarize_roots(records: list[RawArchiveRecord], bundle_roots: list[dict]) -> list[RawRootSummary]:
    full = {item["root"]: item for item in bundle_roots}
    prefixes = sorted({rec.local_path.split("/", 1)[0] for rec in records})
    summaries: list[RawRootSummary] = []
    for prefix in prefixes:
        subset = [r for r in records if r.local_path == prefix or r.local_path.startswith(prefix + "/")]
        full_item = full.get(prefix, {})
        summaries.append(
            RawRootSummary(
                root=prefix,
                referenced_file_count=len(subset),
                referenced_total_bytes=sum(r.size_bytes or 0 for r in subset),
                full_local_root_file_count=full_item.get("file_count"),
                full_local_root_total_bytes=full_item.get("total_bytes"),
            )
        )
    return summaries


def write_markdown(path: Path, payload: dict) -> None:
    lines = [
        "# Raw-data archive manifest",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Git branch: `{payload['git_branch']}`",
        f"Git HEAD: `{payload['git_head']}`",
        "",
        "## Summary",
        "",
        f"- Referenced raw files: {payload['summary']['referenced_file_count']}",
        f"- Missing referenced raw files: {payload['summary']['missing_file_count']}",
        f"- Referenced raw-data bytes: {payload['summary']['referenced_total_bytes']}",
        f"- Suggested archive root: `raw_data/`",
        "",
        "This manifest is a preparation aid only.  It does not create or upload",
        "an archive.  If the target journal requires raw data, archive the files",
        "listed below with their relative paths preserved under `raw_data/`, then",
        "replace or supplement the GitHub-only availability statement with the",
        "archive DOI.",
        "",
        "## Referenced subset versus full local roots",
        "",
        "| Root | Referenced files | Referenced bytes | Full local files | Full local bytes |",
        "|---|---:|---:|---:|---:|",
    ]
    for rec in payload["root_summaries"]:
        lines.append(
            f"| `{rec['root']}` | {rec['referenced_file_count']} | {rec['referenced_total_bytes']} | "
            f"{rec['full_local_root_file_count'] if rec['full_local_root_file_count'] is not None else ''} | "
            f"{rec['full_local_root_total_bytes'] if rec['full_local_root_total_bytes'] is not None else ''} |"
        )
    lines.extend(
        [
            "",
            "## Raw file records",
            "",
            "| Archive path | Local path | Size | SHA-256 | Source artifacts |",
            "|---|---|---:|---|---|",
        ]
    )
    for rec in payload["raw_files"]:
        lines.append(
            f"| `{rec['archive_path']}` | `{rec['local_path']}` | "
            f"{rec['size_bytes'] if rec['size_bytes'] is not None else ''} | "
            f"`{(rec['sha256'] or '')[:16]}` | {', '.join(f'`{s}`' for s in rec['source_artifacts'])} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    root = repo_root()
    revision = root / REVISION_REL
    bundle_path = revision / "submission_bundle_manifest.json"
    bundle = json.loads(bundle_path.read_text())
    records = group_records(bundle["local_raw_dependencies"])
    summaries = summarize_roots(records, bundle.get("local_raw_roots", []))
    missing = [rec for rec in records if not rec.exists]
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision / "scripts/build_raw_data_archive_manifest.py").relative_to(root)),
        "source_manifest": str(bundle_path.relative_to(root)),
        "git_branch": run_git(root, ["branch", "--show-current"]),
        "git_head": run_git(root, ["rev-parse", "HEAD"]),
        "summary": {
            "referenced_file_count": len(records),
            "missing_file_count": len(missing),
            "referenced_total_bytes": sum(rec.size_bytes or 0 for rec in records),
            "status": "PASS" if not missing else "FAIL",
        },
        "root_summaries": [asdict(s) for s in summaries],
        "raw_files": [asdict(r) for r in records],
    }
    json_path = revision / "raw_data_archive_manifest.json"
    md_path = revision / "raw_data_archive_manifest.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps(payload["summary"], indent=2))
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
