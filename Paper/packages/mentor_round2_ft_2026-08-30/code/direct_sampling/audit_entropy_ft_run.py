#!/usr/bin/env python3
"""Strict streaming audit for NLS entropy/action FT sampler outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path


EXPECTED_HEADER = [
    "stream_id",
    "block_id",
    "q_left",
    "q_right",
    "delta_energy",
    "entropy_medium",
    "entropy_rate",
    "action_current",
    "energy_balance_error",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--expected-n", default="10,20,30,40")
    parser.add_argument("--expected-blocks", type=int)
    parser.add_argument("--require-completed", action="store_true")
    parser.add_argument("--require-hash-match", action="store_true")
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument("--balance-rms-limit", type=float, default=5.0e-5)
    return parser.parse_args()


@dataclass
class RunningMoments:
    count: int = 0
    total: float = 0.0
    square_total: float = 0.0

    def add(self, value: float) -> None:
        self.count += 1
        self.total += value
        self.square_total += value * value

    @property
    def mean(self) -> float:
        return self.total / self.count

    @property
    def rms(self) -> float:
        return math.sqrt(self.square_total / self.count)


def close(a: float, b: float, atol: float = 5.0e-11, rtol: float = 5.0e-10) -> bool:
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def one_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one summary row, found {len(rows)}")
    return rows[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_manifest(
    run_dir: Path, require_completed: bool, require_hash_match: bool
) -> dict[str, object]:
    candidates = sorted(run_dir.glob("*manifest.txt"))
    errors: list[str] = []
    if not candidates:
        return {"status": "FAIL", "errors": ["no manifest file found"]}
    manifest = candidates[0]
    text = manifest.read_text()
    if require_completed and "completed_utc=" not in text:
        errors.append("manifest has no completed_utc marker")
    hashes = []
    for line in text.splitlines():
        match = re.match(r"^([0-9a-f]{64})\s+(.+)$", line)
        if not match:
            continue
        recorded, raw_path = match.groups()
        path = Path(raw_path)
        current = sha256(path) if path.exists() else None
        hashes.append(
            {
                "path": str(path),
                "recorded_sha256": recorded,
                "current_sha256": current,
                "matches_current": current == recorded,
            }
        )
        if require_hash_match and current != recorded:
            errors.append(f"source hash mismatch or missing file: {path}")
    return {
        "status": "PASS" if not errors else "FAIL",
        "path": str(manifest),
        "errors": errors,
        "source_hashes": hashes,
    }


def audit_run(
    blocks_path: Path,
    expected_blocks: int | None,
    balance_rms_limit: float,
) -> dict[str, object]:
    summary_path = blocks_path.with_name(blocks_path.name.replace("_blocks.csv", "_summary.csv"))
    if not summary_path.exists():
        return {"status": "FAIL", "blocks": str(blocks_path), "errors": ["missing summary"]}
    summary = one_row(summary_path)
    errors: list[str] = []
    n = int(summary["n"])
    temperature_left = float(summary["T1"])
    temperature_right = float(summary["Tn"])
    block_time = float(summary["block_time"])
    expected_from_summary = int(summary["n_blocks"])
    stream_count = int(summary["n_streams"])
    blocks_per_stream = int(summary["blocks_per_stream"])

    if expected_blocks is not None and expected_from_summary != expected_blocks:
        errors.append(
            f"summary n_blocks={expected_from_summary}, expected {expected_blocks}"
        )
    if expected_from_summary != stream_count * blocks_per_stream:
        errors.append("summary n_blocks != n_streams * blocks_per_stream")
    if summary.get("coordinates") != "cartesian":
        errors.append(f"coordinates={summary.get('coordinates')!r}, expected cartesian")
    if summary.get("energy_convention") != "E=H/2":
        errors.append(
            f"energy_convention={summary.get('energy_convention')!r}, expected E=H/2"
        )
    if int(summary["midpoint_failure_count"]) != 0:
        errors.append(f"midpoint failures={summary['midpoint_failure_count']}")

    q_left_rate = RunningMoments()
    q_right_rate = RunningMoments()
    drift_rate = RunningMoments()
    entropy_rate = RunningMoments()
    action_current = RunningMoments()
    balance_rate = RunningMoments()
    formula_max = {"entropy": 0.0, "entropy_rate": 0.0, "balance": 0.0}
    nonfinite = 0
    ordering_errors = 0
    previous_stream = -1
    previous_block = -1

    with blocks_path.open(newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, None)
        if header != EXPECTED_HEADER:
            errors.append(f"unexpected header: {header}")
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(EXPECTED_HEADER):
                errors.append(f"row {row_number}: expected 9 columns, found {len(row)}")
                continue
            stream_id = int(row[0])
            block_id = int(row[1])
            values = [float(item) for item in row[2:]]
            if not all(math.isfinite(value) for value in values):
                nonfinite += 1
                continue
            q_left, q_right, delta_energy, entropy, reported_entropy_rate, action, balance = values
            expected_entropy = -q_left / temperature_left - q_right / temperature_right
            expected_entropy_rate = entropy / block_time
            expected_balance = q_left + q_right - delta_energy
            formula_max["entropy"] = max(
                formula_max["entropy"], abs(entropy - expected_entropy)
            )
            formula_max["entropy_rate"] = max(
                formula_max["entropy_rate"],
                abs(reported_entropy_rate - expected_entropy_rate),
            )
            formula_max["balance"] = max(
                formula_max["balance"], abs(balance - expected_balance)
            )

            if stream_id == previous_stream:
                if block_id != previous_block + 1:
                    ordering_errors += 1
            elif stream_id == previous_stream + 1:
                if block_id != 0:
                    ordering_errors += 1
            else:
                ordering_errors += 1
            previous_stream = stream_id
            previous_block = block_id

            q_left_rate.add(q_left / block_time)
            q_right_rate.add(q_right / block_time)
            drift_rate.add(delta_energy / block_time)
            entropy_rate.add(reported_entropy_rate)
            action_current.add(action)
            balance_rate.add(balance / block_time)

    observed_blocks = q_left_rate.count
    if observed_blocks != expected_from_summary:
        errors.append(
            f"observed {observed_blocks} finite rows, summary expects {expected_from_summary}"
        )
    if nonfinite:
        errors.append(f"found {nonfinite} rows with NaN/Inf")
    if ordering_errors:
        errors.append(f"found {ordering_errors} stream/block ordering errors")
    if previous_stream != stream_count - 1 or previous_block != blocks_per_stream - 1:
        errors.append(
            "last stream/block does not match summary: "
            f"({previous_stream},{previous_block}) vs "
            f"({stream_count - 1},{blocks_per_stream - 1})"
        )
    for name, maximum in formula_max.items():
        if maximum > 5.0e-10:
            errors.append(f"{name} formula max error {maximum:.6g} exceeds tolerance")

    recomputed = {
        "mean_q_left_rate": q_left_rate.mean,
        "mean_q_right_rate": q_right_rate.mean,
        "mean_energy_drift_rate": drift_rate.mean,
        "mean_entropy_rate": entropy_rate.mean,
        "mean_action_current": action_current.mean,
        "mean_energy_balance_error_rate": balance_rate.mean,
        "rms_energy_balance_error_rate": balance_rate.rms,
    }
    for name, value in recomputed.items():
        reported = float(summary[name])
        if not close(value, reported):
            errors.append(
                f"{name}: recomputed {value:.17g}, summary {reported:.17g}"
            )
    if balance_rate.rms > balance_rms_limit:
        errors.append(
            f"balance RMS rate {balance_rate.rms:.6g} exceeds {balance_rms_limit:.6g}"
        )

    return {
        "status": "PASS" if not errors else "FAIL",
        "n": n,
        "blocks": str(blocks_path),
        "summary": str(summary_path),
        "observed_blocks": observed_blocks,
        "expected_blocks": expected_from_summary,
        "n_streams": stream_count,
        "blocks_per_stream": blocks_per_stream,
        "nonfinite_rows": nonfinite,
        "ordering_errors": ordering_errors,
        "formula_max_abs_error": formula_max,
        "recomputed": recomputed,
        "errors": errors,
    }


def markdown(report: dict[str, object]) -> str:
    lines = ["# Entropy/action FT run audit", ""]
    lines.append(f"Overall status: **{report['status']}**")
    lines.append("")
    manifest = report["manifest"]
    lines.extend(
        [
            "## Manifest",
            "",
            f"Status: **{manifest['status']}**",
            f"Path: `{manifest.get('path', 'missing')}`",
            "",
        ]
    )
    for item in report["runs"]:
        lines.extend(
            [
                f"## n={item.get('n', '?')}",
                "",
                f"Status: **{item['status']}**",
                f"Rows: {item.get('observed_blocks', 'unknown')}",
            ]
        )
        if "recomputed" in item:
            values = item["recomputed"]
            lines.extend(
                [
                    f"Mean action current: {values['mean_action_current']:.10g}",
                    f"Mean entropy rate: {values['mean_entropy_rate']:.10g}",
                    f"Energy-balance RMS rate: {values['rms_energy_balance_error_rate']:.6g}",
                ]
            )
        if item.get("errors"):
            lines.append("")
            lines.extend(f"- {error}" for error in item["errors"])
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    expected_n = [int(value) for value in args.expected_n.split(",") if value]
    runs = []
    for n in expected_n:
        path = args.run_dir / f"n{n}_blocks.csv"
        if not path.exists():
            runs.append(
                {"status": "FAIL", "n": n, "blocks": str(path), "errors": ["missing blocks file"]}
            )
        else:
            runs.append(audit_run(path, args.expected_blocks, args.balance_rms_limit))
    manifest = audit_manifest(
        args.run_dir, args.require_completed, args.require_hash_match
    )
    status = "PASS" if manifest["status"] == "PASS" and all(
        run["status"] == "PASS" for run in runs
    ) else "FAIL"
    report = {
        "status": status,
        "run_dir": str(args.run_dir),
        "manifest": manifest,
        "runs": runs,
    }
    prefix = args.output_prefix or (args.run_dir / "audit")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    prefix.with_suffix(".md").write_text(markdown(report) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if status == "PASS" else 1)


if __name__ == "__main__":
    main()
