#!/usr/bin/env python3
"""Validate paper-module manifests and local navigation links."""

from __future__ import annotations

import csv
from pathlib import Path


MODULES = Path(__file__).resolve().parent
ROOT = MODULES.parents[1]


def resolve_registered_path(raw: str) -> Path:
    path = Path(raw)
    return path if path.is_absolute() else ROOT / path


def main() -> int:
    manifests = sorted(MODULES.glob("[0-9][0-9]_*/manifest.tsv"))
    missing: list[tuple[Path, str, str]] = []
    records = 0

    for manifest in manifests:
        with manifest.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["role", "path", "status", "note"]:
                print(f"BAD HEADER: {manifest.relative_to(ROOT)}")
                return 2
            for row in reader:
                records += 1
                target = resolve_registered_path(row["path"])
                if not target.exists():
                    missing.append((manifest, row["path"], row["status"]))

    broken_links = [
        path
        for path in MODULES.glob("[0-9][0-9]_*/local/**/*")
        if path.is_symlink() and not path.exists()
    ]

    print(f"modules={len(manifests)} registered_paths={records}")
    print(f"missing_registered_paths={len(missing)} broken_local_links={len(broken_links)}")

    for manifest, raw, status in missing:
        print(f"MISSING [{status}] {manifest.parent.name}: {raw}")
    for path in broken_links:
        print(f"BROKEN LINK: {path.relative_to(ROOT)} -> {path.readlink()}")

    return 1 if missing or broken_links else 0


if __name__ == "__main__":
    raise SystemExit(main())
