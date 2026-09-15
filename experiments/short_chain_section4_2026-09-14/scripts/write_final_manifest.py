#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for the Section 4 package."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "provenance" / "FINAL_SHA256.tsv"

EXCLUDED_SUFFIXES = {".aux", ".fdb_latexmk", ".fls", ".log", ".synctex.gz"}


def excluded(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    if rel == Path("provenance/FINAL_SHA256.tsv"):
        return True
    if "__pycache__" in rel.parts:
        return True
    if any(str(rel).endswith(suffix) for suffix in EXCLUDED_SUFFIXES):
        return True
    if path.name.endswith(".pid"):
        return True
    return False


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def main() -> None:
    paths = sorted(path for path in ROOT.rglob("*") if path.is_file() and not excluded(path))
    lines = ["sha256\tbytes\tpath"]
    for path in paths:
        lines.append(f"{digest(path)}\t{path.stat().st_size}\t{path.relative_to(ROOT)}")
    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT} with {len(paths)} files")


if __name__ == "__main__":
    main()

