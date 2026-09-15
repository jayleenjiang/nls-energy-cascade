#!/usr/bin/env python3
"""Audit generated artifacts and build a self-contained Section 4.2 package."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import zipfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"
REPORT = HERE / "report"
FIGURES = HERE / "figures"
PACKAGE_ROOT = HERE / "package"
PACKAGE_NAME = "section4_2_stabilization_rebuild_2026-09-15"

CSV_EXPECTED_ROWS = {
    "support_summary.csv": 12,
    "asymmetry_pointwise.csv": 10,
    "asymmetry_angular_tv.csv": 12,
    "phase_locking.csv": 360,
    "phase_histograms.csv": 2592,
    "current_balance.csv": 36,
    "density_trajectory_comparisons.csv": 62,
    "current_balance_full_trajectory_audit.csv": 6,
    "current_balance_masked_comparison_audit.csv": 2,
    "report_asymmetry_summary.csv": 6,
    "report_phase_summary.csv": 45,
    "report_phase_failures.csv": 3,
    "report_current_summary.csv": 24,
    "report_gate_summary.csv": 5,
}

REQUIRED = [
    HERE / "PROTOCOL.md",
    HERE / "POST_HOC_CURRENT_SCOPE_AUDIT.md",
    HERE / "FINAL_VERDICT.md",
    HERE / "VALIDATION_REPORT.md",
    HERE / "COMMANDS.txt",
    HERE / "VISUAL_QA.md",
    HERE / "PACKAGE_CONTENTS.md",
    HERE / "analyze_stabilization.py",
    HERE / "audit_current_scope.py",
    HERE / "build_report_summary.py",
    HERE / "audit_and_package.py",
    HERE / "run_analysis.sh",
    HERE / "run_current_scope_audit.sh",
    HERE / "run_report_summary.sh",
    REPORT / "stabilization_density_rebuild_report.tex",
    REPORT / "section4_2_stabilization.tex",
    REPORT / "section4_2_smoke.tex",
    REPORT / "stabilization_density_rebuild_report.pdf",
    REPORT / "section4_2_smoke.pdf",
    FIGURES / "exchange_asymmetry_rebuilt.pdf",
    FIGURES / "exchange_asymmetry_rebuilt.png",
    FIGURES / "phase_locking_rebuilt.pdf",
    FIGURES / "phase_locking_rebuilt.png",
    FIGURES / "current_balance_rebuilt.pdf",
    FIGURES / "current_balance_rebuilt.png",
    ANALYSIS / "input_hashes.csv",
    ANALYSIS / "provenance.json",
    ANALYSIS / "verdict.json",
    ANALYSIS / "verdict_frozen_v1.json",
    ANALYSIS / "current_scope_audit_verdict.json",
    ANALYSIS / "report_scalar_summary.json",
] + [ANALYSIS / name for name in CSV_EXPECTED_ROWS]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_csv(path: Path, expected: int) -> dict:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    nonfinite = []
    for row_index, row in enumerate(rows, start=2):
        for key, value in row.items():
            if value in ("", "N/A", "True", "False"):
                continue
            try:
                number = float(value)
            except ValueError:
                continue
            if not math.isfinite(number):
                nonfinite.append({"row": row_index, "column": key, "value": value})
    return {
        "path": str(path.relative_to(HERE)),
        "rows": len(rows),
        "expected_rows": expected,
        "row_count_pass": len(rows) == expected,
        "nonfinite_numeric_values": nonfinite,
        "finite_pass": not nonfinite,
    }


def main() -> None:
    shutil.copy2(
        REPORT / "build/stabilization_density_rebuild_report.pdf",
        REPORT / "stabilization_density_rebuild_report.pdf",
    )
    shutil.copy2(
        REPORT / "build/section4_2_smoke.pdf",
        REPORT / "section4_2_smoke.pdf",
    )
    missing = [str(path) for path in REQUIRED if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing required artifacts: {missing}")

    csv_audits = [audit_csv(ANALYSIS / name, count)
                  for name, count in CSV_EXPECTED_ROWS.items()]
    integrity = {
        "analysis_only": True,
        "new_simulation": False,
        "required_file_count": len(REQUIRED),
        "missing": missing,
        "csv_audits": csv_audits,
        "complete_pass": all(
            row["row_count_pass"] and row["finite_pass"] for row in csv_audits
        ),
    }
    (ANALYSIS / "output_integrity.json").write_text(json.dumps(integrity, indent=2) + "\n")
    if not integrity["complete_pass"]:
        raise RuntimeError("output integrity audit failed")

    hash_targets = sorted(set(REQUIRED + [ANALYSIS / "output_integrity.json"]))
    hash_rows = [{
        "path": str(path.relative_to(HERE)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    } for path in hash_targets]
    with (ANALYSIS / "artifact_hashes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("path", "bytes", "sha256"), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(hash_rows)

    package_dir = PACKAGE_ROOT / PACKAGE_NAME
    if package_dir.exists():
        shutil.rmtree(package_dir)
    package_dir.mkdir(parents=True)
    package_targets = hash_targets + [ANALYSIS / "artifact_hashes.csv"]
    for source in package_targets:
        destination = package_dir / source.relative_to(HERE)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    archive = PACKAGE_ROOT / f"{PACKAGE_NAME}.zip"
    if archive.exists():
        archive.unlink()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(package_dir.rglob("*")):
            if path.is_file():
                handle.write(path, Path(PACKAGE_NAME) / path.relative_to(package_dir))
    archive_hash = sha256(archive)
    (PACKAGE_ROOT / f"{PACKAGE_NAME}.zip.sha256").write_text(
        f"{archive_hash}  {archive.name}\n"
    )
    print(json.dumps({
        "integrity_pass": integrity["complete_pass"],
        "package": str(archive),
        "package_bytes": archive.stat().st_size,
        "package_sha256": archive_hash,
    }, indent=2))


if __name__ == "__main__":
    main()
