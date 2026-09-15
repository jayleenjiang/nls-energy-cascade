#!/usr/bin/env python3
"""Independent consistency audit for direct entropy-SCGF output."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def close(a: float, b: float, tolerance: float = 2.0e-11) -> bool:
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    scgf = rows(args.output_dir / "scgf_points.csv")
    symmetry = rows(args.output_dir / "gc_symmetry_residuals.csv")
    summary = rows(args.output_dir / "scgf_reliability_summary.csv")
    errors: list[str] = []
    checks = 0

    lookup: dict[tuple[int, float, float], dict[str, str]] = {}
    for row in scgf:
        key = (int(row["n"]), float(row["tau"]), float(row["k"]))
        lookup[key] = row
        psi = float(row["psi"])
        ci_low = float(row["psi_ci_low"])
        ci_high = float(row["psi_ci_high"])
        sample_count = int(row["sample_count"])
        sample_ess = float(row["sample_weight_ess"])
        stream_ess = float(row["stream_weight_ess"])
        checks += 6
        if not (math.isfinite(psi) and ci_low <= ci_high):
            errors.append(f"nonfinite or reversed SCGF interval: {key}")
        if not (1.0 <= sample_ess <= sample_count * (1.0 + 1e-10)):
            errors.append(f"sample ESS out of range: {key}")
        if not (1.0 <= stream_ess <= int(row["n_streams"]) * (1.0 + 1e-10)):
            errors.append(f"stream ESS out of range: {key}")
        if not (0.0 < float(row["max_sample_weight_fraction"]) <= 1.0):
            errors.append(f"sample weight fraction out of range: {key}")
        if not (0.0 < float(row["max_stream_weight_fraction"]) <= 1.0):
            errors.append(f"stream weight fraction out of range: {key}")
        if close(float(row["k"]), 0.0) and not close(psi, 0.0, 1e-13):
            errors.append(f"psi(0) is not zero: {key}")

    for row in symmetry:
        n = int(row["n"])
        tau = float(row["tau"])
        k = float(row["k"])
        partner = float(row["one_minus_k"])
        expected = float(lookup[(n, tau, k)]["psi"]) - float(
            lookup[(n, tau, partner)]["psi"]
        )
        observed = float(row["symmetry_residual"])
        reverse = next(
            candidate
            for candidate in symmetry
            if int(candidate["n"]) == n
            and close(float(candidate["tau"]), tau)
            and close(float(candidate["k"]), partner)
        )
        checks += 3
        if not close(observed, expected):
            errors.append(f"residual mismatch: {(n, tau, k)}")
        if not close(observed, -float(reverse["symmetry_residual"])):
            errors.append(f"residual antisymmetry failed: {(n, tau, k)}")
        if close(k, 0.5) and not close(observed, 0.0, 1e-13):
            errors.append(f"midpoint residual is not zero: {(n, tau)}")

    expected_summary_rows = len({(int(r["n"]), float(r["tau"])) for r in scgf})
    checks += 1
    if len(summary) != expected_summary_rows:
        errors.append("summary row count mismatch")

    # Recompute one raw point without importing the analysis implementation.
    target = lookup[(10, 20.0, 0.5)]
    raw_path = Path(target["source"])
    entropy = np.loadtxt(raw_path, delimiter=",", skiprows=1, usecols=(5,))
    z = -0.5 * entropy
    maximum = float(np.max(z))
    direct = (
        maximum + math.log(float(np.exp(z - maximum).sum()))
        - math.log(entropy.size)
    ) / 20.0
    checks += 1
    if not close(direct, float(target["psi"]), 5e-13):
        errors.append("independent raw n=10,t=20,k=0.5 recomputation failed")

    status = "PASS" if not errors else "FAIL"
    report = [
        "# Direct SCGF audit",
        "",
        f"Status: **{status}**",
        f"Checks: {checks}",
        f"Errors: {len(errors)}",
    ]
    if errors:
        report.extend(["", "## Errors", *[f"- {error}" for error in errors]])
    (args.output_dir / "audit.md").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
