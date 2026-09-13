#!/usr/bin/env python3
"""Build a source-only submission bundle from the release manifest.

The source bundle is a practical handoff artifact for coauthors, editors, or a
journal submission portal.  It copies the tracked release files and tracked
files in declared release directories into a fresh staging directory under
``tmp/`` and writes a deterministic ``.tar.gz`` archive plus checksum report.

Large local raw-data roots are intentionally not copied.  They remain governed
by ``raw_data_archive_manifest.json`` and any future Zenodo/OSF/journal archive
decision.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import json
import shutil
import subprocess
import tarfile
from dataclasses import asdict, dataclass
from pathlib import Path


REVISION_REL = Path("Paper/revision_2026-06-19")
DEFAULT_RUNS_REL = Path("tmp/submission_source_bundle/runs")
ARCHIVE_ROOT = "NLS_numerical_study_source"
MANIFEST_REL = REVISION_REL / "submission_bundle_manifest.json"
REPORT_JSON_REL = REVISION_REL / "submission_source_bundle_report.json"
REPORT_MD_REL = REVISION_REL / "submission_source_bundle_report.md"
CONTENTS_SHA_REL = Path("SOURCE_BUNDLE_CONTENTS.sha256")

# These files are generated during or after the local gate.  Excluding them from
# the archive avoids self-referential bundles whose checksums depend on the
# report about the bundle itself.  They remain tracked handoff files in Git.
SELF_REFERENTIAL_OUTPUTS = {
    str(REVISION_REL / "submission_bundle_manifest.json"),
    str(REVISION_REL / "submission_bundle_manifest.md"),
    str(REVISION_REL / "submission_checks_summary.json"),
    str(REVISION_REL / "submission_checks_summary.md"),
    str(REPORT_JSON_REL),
    str(REPORT_MD_REL),
}


@dataclass
class BundleFile:
    path: str
    source: str
    size_bytes: int
    sha256: str


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def run_git(root: Path, args: list[str]) -> str:
    proc = subprocess.run(["git", *args], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout


def git_files_under(root: Path, rel_dir: str) -> list[str]:
    prefix = rel_dir.rstrip("/") + "/"
    return [line for line in run_git(root, ["ls-files", prefix]).splitlines() if line]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_file(root: Path, rel: str, staging_root: Path, source: str) -> BundleFile:
    src = root / rel
    dst = staging_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return BundleFile(path=rel, source=source, size_bytes=dst.stat().st_size, sha256=sha256_file(dst))


def deterministic_tar_gz(staging_root: Path, archive_path: Path) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w") as tar:
                for path in sorted(staging_root.rglob("*")):
                    rel = path.relative_to(staging_root)
                    info = tar.gettarinfo(str(path), arcname=str(Path(ARCHIVE_ROOT) / rel))
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    info.mtime = 0
                    if path.is_file():
                        with path.open("rb") as f:
                            tar.addfile(info, f)
                    else:
                        tar.addfile(info)


def write_contents_checksum(staging_root: Path, files: list[BundleFile]) -> BundleFile:
    contents = staging_root / CONTENTS_SHA_REL
    lines = [f"{rec.sha256}  {rec.path}" for rec in sorted(files, key=lambda item: item.path)]
    contents.write_text("\n".join(lines) + "\n")
    return BundleFile(
        path=str(CONTENTS_SHA_REL),
        source="generated-checksum-index",
        size_bytes=contents.stat().st_size,
        sha256=sha256_file(contents),
    )


def write_markdown(path: Path, payload: dict) -> None:
    summary = payload["summary"]
    lines = [
        "# Submission source bundle report",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "## Archive",
        "",
        f"- Archive path: `{payload['archive']['path']}`",
        f"- Archive bytes: {payload['archive']['size_bytes']}",
        f"- Archive SHA-256: `{payload['archive']['sha256']}`",
        f"- Staging directory: `{payload['staging_directory']}`",
        "",
        "## Summary",
        "",
        "| Metric | Count |",
        "|---|---:|",
        f"| Manifest release-file records | {summary['manifest_release_file_records']} |",
        f"| Manifest release directories | {summary['manifest_release_directories']} |",
        f"| Included regular files | {summary['included_file_count']} |",
        f"| Included bytes | {summary['included_total_bytes']} |",
        f"| Excluded volatile/self-referential files | {summary['excluded_self_referential_count']} |",
        f"| Missing files | {summary['missing_file_count']} |",
        f"| Directory-tracked files copied | {summary['directory_tracked_file_count']} |",
        f"| Local raw dependency records excluded | {summary['local_raw_dependency_count']} |",
        "",
        "## Excluded self-referential files",
        "",
    ]
    if payload["excluded_self_referential_files"]:
        for rel in payload["excluded_self_referential_files"]:
            lines.append(f"- `{rel}`")
    else:
        lines.append("- None.")
    lines.extend(["", "## Missing files", ""])
    if payload["missing_files"]:
        for rel in payload["missing_files"]:
            lines.append(f"- `{rel}`")
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- The archive is source-only and intentionally excludes the large local raw-data roots.",
            "- If a target journal requires raw data, use `raw_data_archive_manifest.md` to prepare a DOI-backed raw-data supplement.",
            "- The archive is written under `tmp/` and is not intended to be committed to Git.",
            "- The file `SOURCE_BUNDLE_CONTENTS.sha256` inside the archive records SHA-256 checksums for copied files.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-runs-dir", default=str(DEFAULT_RUNS_REL), help="Directory under the repo root for timestamped bundle runs.")
    args = parser.parse_args()

    root = repo_root()
    manifest_path = root / MANIFEST_REL
    manifest = json.loads(manifest_path.read_text())

    now = _dt.datetime.now(_dt.timezone.utc)
    run_id = now.strftime("%Y%m%dT%H%M%SZ")
    run_dir = root / args.output_runs_dir / run_id
    staging_root = run_dir / ARCHIVE_ROOT
    archive_path = run_dir / f"{ARCHIVE_ROOT}.tar.gz"
    staging_root.mkdir(parents=True, exist_ok=False)

    missing: list[str] = []
    excluded: list[str] = []
    copied: list[BundleFile] = []
    seen: set[str] = set()

    for rec in manifest["release_files"]:
        rel = rec["path"]
        if rel in SELF_REFERENTIAL_OUTPUTS:
            excluded.append(rel)
            continue
        if rel in seen:
            continue
        seen.add(rel)
        if not (root / rel).is_file():
            missing.append(rel)
            continue
        copied.append(copy_file(root, rel, staging_root, "release-file"))

    directory_file_count = 0
    for directory in manifest["release_directories"]:
        rel_dir = directory["path"]
        for rel in git_files_under(root, rel_dir):
            if rel in SELF_REFERENTIAL_OUTPUTS:
                excluded.append(rel)
                continue
            if rel in seen:
                continue
            seen.add(rel)
            if not (root / rel).is_file():
                missing.append(rel)
                continue
            copied.append(copy_file(root, rel, staging_root, f"release-directory:{rel_dir}"))
            directory_file_count += 1

    contents_record = write_contents_checksum(staging_root, copied)
    copied_with_index = [*copied, contents_record]
    deterministic_tar_gz(staging_root, archive_path)

    status = "PASS" if not missing else "FAIL"
    payload = {
        "generated_at_utc": now.isoformat(),
        "script": str((root / REVISION_REL / "scripts/build_submission_source_bundle.py").relative_to(root)),
        "status": status,
        "manifest": str(MANIFEST_REL),
        "staging_directory": str(staging_root.relative_to(root)),
        "archive": {
            "path": str(archive_path.relative_to(root)),
            "size_bytes": archive_path.stat().st_size,
            "sha256": sha256_file(archive_path),
        },
        "summary": {
            "manifest_release_file_records": len(manifest["release_files"]),
            "manifest_release_directories": len(manifest["release_directories"]),
            "included_file_count": len(copied_with_index),
            "included_total_bytes": sum(rec.size_bytes for rec in copied_with_index),
            "excluded_self_referential_count": len(sorted(set(excluded))),
            "missing_file_count": len(missing),
            "directory_tracked_file_count": directory_file_count,
            "local_raw_dependency_count": manifest["summary"]["local_raw_dependency_count"],
        },
        "excluded_self_referential_files": sorted(set(excluded)),
        "missing_files": missing,
        "files": [asdict(rec) for rec in copied_with_index],
    }

    json_path = root / REPORT_JSON_REL
    md_path = root / REPORT_MD_REL
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(md_path, payload)
    print(json.dumps({"status": status, "summary": payload["summary"], "archive": payload["archive"]}, indent=2, sort_keys=True))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
