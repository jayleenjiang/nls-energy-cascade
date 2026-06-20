#!/usr/bin/env python3
"""Build a local journal-upload package from the current verified artifacts.

The output is local-only under ``tmp/``.  It does not submit the manuscript,
upload data, create a DOI, or modify the manuscript source.  Its purpose is to
put the current PDF, source archive, and handoff documents in one timestamped
directory with checksums so the authors/advisor can see exactly which local
files are ready for a journal upload route.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import json
import shutil
import tarfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
REVISION_REL = Path("Paper/revision_2026-06-19")
RUNS_REL = Path("tmp/journal_upload_package/runs")
ARCHIVE_ROOT = "NLS_numerical_study_journal_upload_package"
CONTENTS_SHA_REL = Path("UPLOAD_PACKAGE_CONTENTS.sha256")


@dataclass
class PackageFile:
    package_path: str
    source_path: str
    role: str
    size_bytes: int
    sha256: str


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(rel: str | Path) -> dict[str, Any]:
    return json.loads((ROOT / rel).read_text())


def copy_into_package(source: Path, staging_root: Path, package_rel: str, role: str) -> PackageFile:
    if not source.is_file():
        raise FileNotFoundError(f"Required upload-package source not found: {source.relative_to(ROOT)}")
    dst = staging_root / package_rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dst)
    return PackageFile(
        package_path=package_rel,
        source_path=str(source.relative_to(ROOT)),
        role=role,
        size_bytes=dst.stat().st_size,
        sha256=sha256_file(dst),
    )


def write_contents_checksum(staging_root: Path, files: list[PackageFile]) -> PackageFile:
    contents = staging_root / CONTENTS_SHA_REL
    lines = [f"{rec.sha256}  {rec.package_path}" for rec in sorted(files, key=lambda item: item.package_path)]
    contents.write_text("\n".join(lines) + "\n")
    return PackageFile(
        package_path=str(CONTENTS_SHA_REL),
        source_path="generated checksum index",
        role="checksum-index",
        size_bytes=contents.stat().st_size,
        sha256=sha256_file(contents),
    )


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


def write_package_manifest(staging_root: Path, payload: dict[str, Any]) -> list[PackageFile]:
    json_path = staging_root / "UPLOAD_PACKAGE_MANIFEST.json"
    md_path = staging_root / "UPLOAD_PACKAGE_MANIFEST.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Journal upload package manifest",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        f"Route: `{payload['route']}`",
        "",
        f"Status: **{payload['status']}**",
        "",
        "This package is local-only.  It does not certify author approval, target-journal submission,",
        "professional similarity screening, DOI creation, or raw-data upload.",
        "",
        "## Files",
        "",
        "| Package path | Role | Source | Size bytes | SHA-256 |",
        "|---|---|---|---:|---|",
    ]
    for rec in payload["files"]:
        lines.append(
            f"| `{rec['package_path']}` | {rec['role']} | `{rec['source_path']}` | "
            f"{rec['size_bytes']} | `{rec['sha256']}` |"
        )
    lines.extend(
        [
            "",
            "## Remaining external blockers",
            "",
            "- Final author metadata, declarations, and author approval.",
            "- Target-journal confirmation and any journal-template conversion.",
            "- Professional similarity/self-plagiarism screening.",
            "- Repository release URL and/or DOI-backed raw-data archive, depending on route.",
            "- Final page-by-page PDF review after all final edits.",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n")
    return [
        PackageFile(
            package_path="UPLOAD_PACKAGE_MANIFEST.json",
            source_path="generated package manifest",
            role="package-manifest",
            size_bytes=json_path.stat().st_size,
            sha256=sha256_file(json_path),
        ),
        PackageFile(
            package_path="UPLOAD_PACKAGE_MANIFEST.md",
            source_path="generated package manifest",
            role="package-manifest",
            size_bytes=md_path.stat().st_size,
            sha256=sha256_file(md_path),
        ),
    ]


def pdf_record_for_route(pdf_artifacts: dict[str, Any], route: str) -> dict[str, Any]:
    pdfs = {rec["id"]: rec for rec in pdf_artifacts["pdfs"]}
    if route.startswith("siads"):
        return pdfs["siads_review_pdf"]
    return pdfs["generic_manuscript_pdf"]


def build_package(route: str, include_raw_data: bool, output_runs_dir: Path) -> dict[str, Any]:
    pdf_artifacts = load_json(REVISION_REL / "compiled_pdf_artifact_audit.json")
    source_bundle = load_json(REVISION_REL / "submission_source_bundle_report.json")
    raw_build = load_json(REVISION_REL / "raw_data_archive_build_report.json")
    checks = load_json(REVISION_REL / "submission_checks_summary.json")

    now = _dt.datetime.now(_dt.timezone.utc)
    run_id = now.strftime("%Y%m%dT%H%M%SZ")
    run_dir = ROOT / output_runs_dir / run_id
    staging_root = run_dir / ARCHIVE_ROOT
    archive_path = run_dir / f"{ARCHIVE_ROOT}.tar.gz"
    staging_root.mkdir(parents=True, exist_ok=False)

    selected_pdf = pdf_record_for_route(pdf_artifacts, route)
    files: list[PackageFile] = []
    pdf_name = "draft_siads_review.pdf" if route.startswith("siads") else "draft.pdf"
    files.append(copy_into_package(ROOT / selected_pdf["pdf"], staging_root, f"manuscript/{pdf_name}", "manuscript-pdf"))

    source_archive = ROOT / source_bundle["archive"]["path"]
    files.append(copy_into_package(source_archive, staging_root, f"source/{source_archive.name}", "source-bundle"))

    handoff_docs = [
        "journal_upload_file_index_2026-06-20.md",
        "final_author_submission_fields_request_2026-06-20.md",
        "siads_first_submission_packet_2026-06-20.md" if route.startswith("siads") else "target_journal_shortlist_2026-06-19.md",
        "submission_checks_summary.md",
        "submission_metadata_consistency_audit.md",
        "compiled_pdf_artifact_audit.md",
        "submission_source_bundle_report.md",
        "raw_data_archive_build_report.md",
    ]
    for name in handoff_docs:
        files.append(copy_into_package(ROOT / REVISION_REL / name, staging_root, f"handoff/{name}", "handoff-document"))

    raw_archive_record: dict[str, Any] | None = None
    if include_raw_data:
        raw_archive = ROOT / raw_build["archive"]["path"]
        files.append(copy_into_package(raw_archive, staging_root, f"raw_data/{raw_archive.name}", "raw-data-archive"))
        raw_archive_record = raw_build["archive"]

    checksum_record = write_contents_checksum(staging_root, files)
    payload_files = [*files, checksum_record]
    payload = {
        "generated_at_utc": now.isoformat(),
        "script": str((ROOT / REVISION_REL / "scripts/build_journal_upload_package.py").relative_to(ROOT)),
        "route": route,
        "include_raw_data": include_raw_data,
        "status": "PASS",
        "selected_pdf": selected_pdf,
        "source_bundle_archive": source_bundle["archive"],
        "raw_data_archive": raw_archive_record,
        "submission_checks_overall_status": checks["overall_status"],
        "files": [asdict(rec) for rec in payload_files],
        "summary": {
            "package_file_count": len(payload_files),
            "package_total_bytes": sum(rec.size_bytes for rec in payload_files),
            "contains_raw_data_archive": include_raw_data,
        },
        "notes": [
            "Local package only; no journal submission or external upload was performed.",
            "Regenerate after any final author/journal/declaration/data-release edits.",
            "UPLOAD_PACKAGE_CONTENTS.sha256 records checksums for package inputs; the manifest itself is not self-checksummed.",
        ],
    }
    write_package_manifest(staging_root, payload)
    deterministic_tar_gz(staging_root, archive_path)
    payload["archive"] = {
        "path": str(archive_path.relative_to(ROOT)),
        "size_bytes": archive_path.stat().st_size,
        "sha256": sha256_file(archive_path),
    }
    payload["staging_directory"] = str(staging_root.relative_to(ROOT))
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--route",
        choices=["siads-repository", "siads-doi-raw", "generic-repository"],
        default="siads-repository",
        help="Upload route to package locally.",
    )
    parser.add_argument("--include-raw-data", action="store_true", help="Copy the local minimal raw-data archive into the package.")
    parser.add_argument("--output-runs-dir", default=str(RUNS_REL), help="Directory under the repo root for timestamped upload packages.")
    args = parser.parse_args()

    include_raw = args.include_raw_data or args.route.endswith("doi-raw")
    payload = build_package(args.route, include_raw, Path(args.output_runs_dir))
    print(json.dumps({"status": payload["status"], "archive": payload["archive"], "summary": payload["summary"]}, indent=2, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
