#!/usr/bin/env python3
"""Build a source/release manifest for the revised NLS manuscript.

The manifest has two deliberately separate layers:

1. ``release_files``: files that are expected to be present in the GitHub
   source bundle for compiling/auditing the manuscript.  These must exist and
   be git-tracked.
2. ``local_raw_dependencies``: larger raw or historical data files referenced
   by source-trace JSON artifacts.  These are checked for local existence and
   git-tracking status, but are not forced into the GitHub source bundle because
   the largest local roots are hundreds of MB to multiple GB.  They should be
   handled by a deliberate archive/Zenodo/OSF policy if full raw-data release is
   required.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable


REVISION_REL = Path("Paper/revision_2026-06-19")
RAW_ROOT_PREFIXES = ("Energy Cascade/", "KDE/", "lte/")
SELF_GENERATED_MANIFEST_FILES = {
    str(REVISION_REL / "submission_bundle_manifest.json"),
    str(REVISION_REL / "submission_bundle_manifest.md"),
}
VOLATILE_GENERATED_FILES = SELF_GENERATED_MANIFEST_FILES | {
    str(REVISION_REL / "submission_checks_summary.json"),
    str(REVISION_REL / "submission_checks_summary.md"),
    str(REVISION_REL / "submission_source_bundle_report.json"),
    str(REVISION_REL / "submission_source_bundle_report.md"),
}


@dataclass
class FileRecord:
    path: str
    roles: list[str]
    exists: bool
    git_tracked: bool
    size_bytes: int | None
    sha256: str | None


@dataclass
class DirectoryRecord:
    path: str
    roles: list[str]
    exists: bool
    git_tracked_file_count: int
    git_tracked_total_bytes: int
    tree_sha256: str | None


@dataclass
class RawDependencyRecord:
    path: str
    source_artifact: str
    exists: bool
    git_tracked: bool
    size_bytes: int | None


@dataclass
class RawRootRecord:
    root: str
    exists: bool
    file_count: int
    total_bytes: int
    git_tracked_file_count: int


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def run_git(root: Path, args: list[str]) -> str:
    proc = subprocess.run(["git", *args], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout


def git_tracked(root: Path, relpath: str) -> bool:
    proc = subprocess.run(
        ["git", "ls-files", "--error-unmatch", relpath],
        cwd=root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return proc.returncode == 0


def git_files_under(root: Path, rel_dir: str) -> list[str]:
    prefix = rel_dir.rstrip("/") + "/"
    out = run_git(root, ["ls-files", prefix])
    return [line for line in out.splitlines() if line]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_digest(root: Path, files: Iterable[str]) -> str:
    h = hashlib.sha256()
    for rel in sorted(files):
        p = root / rel
        h.update(rel.encode())
        h.update(b"\0")
        h.update(sha256_file(p).encode())
        h.update(b"\0")
    return h.hexdigest()


def add_role(mapping: dict[str, set[str]], path: str, role: str) -> None:
    mapping.setdefault(path, set()).add(role)


def extract_tex_paths(root: Path, revision_dir: Path, release: dict[str, set[str]]) -> None:
    draft = revision_dir / "draft.tex"
    tex = draft.read_text()
    add_role(release, str(draft.relative_to(root)), "manuscript-source")
    siads_review = revision_dir / "draft_siads_review.tex"
    if siads_review.exists():
        add_role(release, str(siads_review.relative_to(root)), "siads-review-source")
    for match in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex):
        rel = str((revision_dir / match.group(1)).relative_to(root))
        add_role(release, rel, "manuscript-figure")


def collect_availability_paths(root: Path, revision_dir: Path, release: dict[str, set[str]], directories: dict[str, set[str]]) -> None:
    audit = json.loads((revision_dir / "availability_path_audit.json").read_text())
    add_role(release, str((revision_dir / "availability_path_audit.json").relative_to(root)), "availability-audit")
    add_role(release, str((revision_dir / "availability_path_audit.md").relative_to(root)), "availability-audit")
    for rec in audit["records"]:
        resolved = rec["resolved_path"]
        if rec["path_type"] == "file":
            add_role(release, resolved, "availability-declared")
        elif rec["path_type"] == "directory":
            directories.setdefault(resolved, set()).add("availability-declared-directory")


def collect_claim_evidence(root: Path, revision_dir: Path, release: dict[str, set[str]]) -> None:
    claim = json.loads((revision_dir / "manuscript_claim_audit.json").read_text())
    add_role(release, str((revision_dir / "manuscript_claim_audit.json").relative_to(root)), "claim-audit")
    add_role(release, str((revision_dir / "manuscript_claim_audit.md").relative_to(root)), "claim-audit")
    for item in claim["claims"]:
        for evidence in item.get("evidence", []):
            add_role(release, evidence, f"claim-evidence:{item['id']}")


def collect_handoff_docs(root: Path, revision_dir: Path, release: dict[str, set[str]]) -> None:
    for name in [
        "references.bib",
        "progress_report.md",
        "integrity_audit_2026-06-19.md",
        "submission_readiness_checklist_2026-06-19.md",
        "final_pre_submission_audit_2026-06-20.md",
        "final_submission_decision_sheet_2026-06-20.md",
        "originality_spotcheck_2026-06-19.md",
        "material_inventory.md",
        "revision_roadmap.md",
        "audit_report.md",
        "submission_bundle_manifest.json",
        "submission_bundle_manifest.md",
        "raw_data_archive_manifest.json",
        "raw_data_archive_manifest.md",
        "raw_data_archive_build_report.json",
        "raw_data_archive_build_report.md",
        "author_submission_fields_audit.json",
        "author_submission_fields_audit.md",
        "reference_integrity_audit.json",
        "reference_integrity_audit.md",
        "submission_checks_summary.json",
        "submission_checks_summary.md",
        "submission_source_bundle_report.json",
        "submission_source_bundle_report.md",
        "submission_metadata_consistency_audit.json",
        "submission_metadata_consistency_audit.md",
        "siads_cover_letter_template.tex",
        "siads_cover_letter_template_build.json",
        "siads_cover_letter_template_build.md",
        "gamma_robustness_smoke_report.json",
        "gamma_robustness_smoke_report.md",
        "compiled_pdf_artifact_audit.json",
        "compiled_pdf_artifact_audit.md",
        "pdf_layout_qa_2026-06-19.md",
        "pre_submission_reviewer_audit_2026-06-21.md",
        "author_submission_action_packet_2026-06-19.md",
        "final_author_submission_fields_request_2026-06-20.md",
        "target_journal_shortlist_2026-06-19.md",
        "siads_first_submission_packet_2026-06-20.md",
        "journal_upload_file_index_2026-06-20.md",
        "author_submission_fields_template.json",
        "submission_reproducibility_readme_2026-06-19.md",
    ]:
        add_role(release, str((revision_dir / name).relative_to(root)), "handoff-document")
    add_role(
        release,
        str((revision_dir / "scripts/build_submission_bundle_manifest.py").relative_to(root)),
        "manifest-builder",
    )
    add_role(
        release,
        str((revision_dir / "scripts/build_raw_data_archive_manifest.py").relative_to(root)),
        "manifest-builder",
    )
    add_role(
        release,
        str((revision_dir / "scripts/build_raw_data_archive.py").relative_to(root)),
        "raw-data-archive-builder",
    )
    add_role(
        release,
        str((revision_dir / "scripts/audit_author_submission_fields.py").relative_to(root)),
        "author-submission-audit",
    )
    add_role(
        release,
        str((revision_dir / "scripts/apply_author_submission_fields.py").relative_to(root)),
        "author-field-applier",
    )
    add_role(
        release,
        str((revision_dir / "scripts/run_submission_checks.py").relative_to(root)),
        "submission-check-runner",
    )
    add_role(
        release,
        str((revision_dir / "scripts/build_submission_source_bundle.py").relative_to(root)),
        "source-bundle-builder",
    )
    add_role(
        release,
        str((revision_dir / "scripts/build_journal_upload_package.py").relative_to(root)),
        "journal-upload-package-builder",
    )
    add_role(
        release,
        str((revision_dir / "scripts/build_siads_cover_letter_template.py").relative_to(root)),
        "cover-letter-template-builder",
    )
    add_role(
        release,
        str((revision_dir / "scripts/run_gamma_robustness_smoke.py").relative_to(root)),
        "gamma-robustness-smoke-runner",
    )
    add_role(
        release,
        str((revision_dir / "scripts/run_gamma_robustness_production.py").relative_to(root)),
        "gamma-robustness-production-runner",
    )
    add_role(
        release,
        str((revision_dir / "scripts/audit_submission_metadata_consistency.py").relative_to(root)),
        "submission-metadata-audit",
    )
    add_role(
        release,
        str((revision_dir / "scripts/audit_references.py").relative_to(root)),
        "reference-audit",
    )
    add_role(
        release,
        str((revision_dir / "scripts/audit_compiled_pdfs.py").relative_to(root)),
        "compiled-pdf-audit",
    )


def find_paths_in_json(obj: Any, *, path_context: bool = False) -> Iterable[str]:
    if isinstance(obj, dict):
        for key, value in obj.items():
            lower = key.lower()
            if isinstance(value, str) and (lower.endswith("path") or "source" in lower or "file" in lower):
                yield value
            child_context = path_context or lower.endswith("paths") or "source" in lower or "file" in lower
            yield from find_paths_in_json(value, path_context=child_context)
    elif isinstance(obj, list):
        for value in obj:
            if path_context and isinstance(value, str):
                yield value
            else:
                yield from find_paths_in_json(value, path_context=path_context)


def collect_raw_dependencies(root: Path, revision_dir: Path) -> list[RawDependencyRecord]:
    records: dict[tuple[str, str], RawDependencyRecord] = {}
    for artifact_name in [
        "manuscript_figure_metrics.json",
        "source_trace_metrics.json",
        "report_assets/compare_residual_mesh_metrics.json",
    ]:
        artifact = revision_dir / artifact_name
        data = json.loads(artifact.read_text())
        for value in find_paths_in_json(data):
            if not any(value.startswith(prefix) for prefix in RAW_ROOT_PREFIXES):
                continue
            rel = value
            path = root / rel
            key = (rel, artifact_name)
            records[key] = RawDependencyRecord(
                path=rel,
                source_artifact=str(artifact.relative_to(root)),
                exists=path.exists(),
                git_tracked=git_tracked(root, rel) if path.is_file() else False,
                size_bytes=path.stat().st_size if path.is_file() else None,
            )
    return sorted(records.values(), key=lambda r: (r.source_artifact, r.path))


def collect_raw_roots(root: Path) -> list[RawRootRecord]:
    out: list[RawRootRecord] = []
    for prefix in RAW_ROOT_PREFIXES:
        rel_root = prefix.rstrip("/")
        path = root / rel_root
        files = [p for p in path.rglob("*") if p.is_file()] if path.is_dir() else []
        tracked = git_files_under(root, rel_root) if path.is_dir() else []
        out.append(
            RawRootRecord(
                root=rel_root,
                exists=path.is_dir(),
                file_count=len(files),
                total_bytes=sum(p.stat().st_size for p in files),
                git_tracked_file_count=len(tracked),
            )
        )
    return out


def make_file_records(root: Path, release: dict[str, set[str]]) -> list[FileRecord]:
    records: list[FileRecord] = []
    for rel, roles in sorted(release.items()):
        path = root / rel
        if path.is_dir():
            continue
        exists = path.exists()
        records.append(
            FileRecord(
                path=rel,
                roles=sorted(roles),
                exists=exists,
                git_tracked=git_tracked(root, rel) if exists else False,
                size_bytes=None if rel in VOLATILE_GENERATED_FILES else (path.stat().st_size if exists and path.is_file() else None),
                sha256=None if rel in VOLATILE_GENERATED_FILES else (sha256_file(path) if exists and path.is_file() else None),
            )
        )
    return records


def make_directory_records(root: Path, directories: dict[str, set[str]]) -> list[DirectoryRecord]:
    out: list[DirectoryRecord] = []
    for rel, roles in sorted(directories.items()):
        path = root / rel
        files = git_files_under(root, rel) if path.is_dir() else []
        total = sum((root / f).stat().st_size for f in files if (root / f).is_file())
        out.append(
            DirectoryRecord(
                path=rel,
                roles=sorted(roles),
                exists=path.is_dir(),
                git_tracked_file_count=len(files),
                git_tracked_total_bytes=total,
                tree_sha256=tree_digest(root, files) if files else None,
            )
        )
    return out


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    file_records = payload["release_files"]
    dir_records = payload["release_directories"]
    raw_records = payload["local_raw_dependencies"]
    raw_roots = payload["local_raw_roots"]
    lines = [
        "# Submission bundle manifest",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Git branch: `{payload['git_branch']}`",
        f"Git HEAD at manifest generation: `{payload['git_head']}`",
        "",
        "## Summary",
        "",
        f"- Release files: {payload['summary']['release_file_count']}",
        f"- Release directories: {payload['summary']['release_directory_count']}",
        f"- Missing release files: {payload['summary']['missing_release_files']}",
        f"- Untracked release files: {payload['summary']['untracked_release_files']}",
        f"- Local raw dependency records: {payload['summary']['local_raw_dependency_count']}",
        f"- Untracked local raw dependency records: {payload['summary']['untracked_local_raw_dependencies']}",
        f"- Manifest status: **{payload['summary']['status']}**",
        "",
        "## Release directories",
        "",
        "| Path | Roles | Tracked files | Total bytes | Tree SHA-256 |",
        "|---|---|---:|---:|---|",
    ]
    for rec in dir_records:
        lines.append(
            f"| `{rec['path']}` | {', '.join(rec['roles'])} | {rec['git_tracked_file_count']} | "
            f"{rec['git_tracked_total_bytes']} | `{(rec['tree_sha256'] or '')[:12]}` |"
        )
    lines.extend(
        [
            "",
            "## Release file records",
            "",
            "| Path | Roles | Size | SHA-256 |",
            "|---|---|---:|---|",
        ]
    )
    for rec in file_records:
        lines.append(
            f"| `{rec['path']}` | {', '.join(rec['roles'][:3])}{' ...' if len(rec['roles']) > 3 else ''} | "
            f"{rec['size_bytes'] if rec['size_bytes'] is not None else ''} | `{(rec['sha256'] or '')[:12]}` |"
        )
    lines.extend(
        [
            "",
            "## Local raw-data dependency limitation",
            "",
            "The source-trace JSON files also point to local raw-data roots such as",
            "`Energy Cascade/`, `KDE/`, and `lte/`.  These local files are checked",
            "for existence here but are not all git-tracked; the largest local root",
            "is multi-GB.  If the target journal requires complete raw-data release,",
            "create an archival DOI-backed supplement and rerun this manifest against",
            "that archive.",
            "",
            "### Local raw roots",
            "",
            "| Root | Exists | Files | Total bytes | Git-tracked files |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for rec in raw_roots:
        lines.append(
            f"| `{rec['root']}` | {str(rec['exists']).lower()} | {rec['file_count']} | "
            f"{rec['total_bytes']} | {rec['git_tracked_file_count']} |"
        )
    lines.extend(
        [
            "",
            "### Referenced local raw dependency records",
            "",
            "| Source artifact | Raw path | Exists | Git-tracked | Size |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for rec in raw_records:
        lines.append(
            f"| `{rec['source_artifact']}` | `{rec['path']}` | {str(rec['exists']).lower()} | "
            f"{str(rec['git_tracked']).lower()} | {rec['size_bytes'] if rec['size_bytes'] is not None else ''} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    root = repo_root()
    revision_dir = root / REVISION_REL
    release: dict[str, set[str]] = {}
    directories: dict[str, set[str]] = {}
    collect_handoff_docs(root, revision_dir, release)
    extract_tex_paths(root, revision_dir, release)
    collect_availability_paths(root, revision_dir, release, directories)
    collect_claim_evidence(root, revision_dir, release)
    file_records = make_file_records(root, release)
    dir_records = make_directory_records(root, directories)
    raw_records = collect_raw_dependencies(root, revision_dir)
    raw_roots = collect_raw_roots(root)

    missing = [r for r in file_records if not r.exists]
    untracked = [r for r in file_records if r.exists and not r.git_tracked]
    raw_untracked = [r for r in raw_records if r.exists and not r.git_tracked]
    status = "PASS" if not missing and not untracked else "FAIL"
    if status == "PASS" and raw_untracked:
        status = "PASS_WITH_LOCAL_RAW_DATA_LIMITATION"

    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": str((revision_dir / "scripts/build_submission_bundle_manifest.py").relative_to(root)),
        "git_branch": run_git(root, ["branch", "--show-current"]).strip(),
        "git_head": run_git(root, ["rev-parse", "HEAD"]).strip(),
        "summary": {
            "release_file_count": len(file_records),
            "release_directory_count": len(dir_records),
            "missing_release_files": len(missing),
            "untracked_release_files": len(untracked),
            "local_raw_dependency_count": len(raw_records),
            "untracked_local_raw_dependencies": len(raw_untracked),
            "status": status,
        },
        "release_files": [asdict(r) for r in file_records],
        "release_directories": [asdict(r) for r in dir_records],
        "local_raw_dependencies": [asdict(r) for r in raw_records],
        "local_raw_roots": [asdict(r) for r in raw_roots],
    }
    json_path = revision_dir / "submission_bundle_manifest.json"
    md_path = revision_dir / "submission_bundle_manifest.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps(payload["summary"], indent=2))
    return 1 if missing or untracked else 0


if __name__ == "__main__":
    raise SystemExit(main())
