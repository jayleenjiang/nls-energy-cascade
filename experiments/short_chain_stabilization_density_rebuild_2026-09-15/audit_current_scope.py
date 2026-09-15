#!/usr/bin/env python3
"""Audit the scope of the frozen current-balance gate without density extrapolation."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"
DENSITY_SCRIPTS = (
    ROOT / "experiments/short_chain_density_ratio_constrained_2026-09-14/scripts"
)
sys.path.insert(0, str(DENSITY_SCRIPTS))

from v6_common import evaluation_roles, load_rows  # noqa: E402


BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260915042


def derived_seed(label: str) -> int:
    digest = hashlib.sha256(f"{BOOTSTRAP_SEED}|{label}".encode()).hexdigest()
    return int(digest[:16], 16)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def direct_full_rows(case: str) -> list[dict]:
    target_rows, _ = evaluation_roles(case, "fine", "test")
    _, streams = load_rows(target_rows)
    sums = np.zeros((len(streams), 3), dtype=np.float64)
    counts = np.zeros(len(streams), dtype=np.float64)
    for stream_index, state in enumerate(streams):
        j12 = state[:, 0] * state[:, 1] * np.sin(state[:, 3])
        j23 = state[:, 1] * state[:, 2] * np.sin(state[:, 4])
        sums[stream_index] = (j12.sum(), j23.sum(), (j12 + j23).sum())
        counts[stream_index] = len(state)

    point = sums.sum(axis=0) / counts.sum()
    rng = np.random.default_rng(derived_seed(f"full_current|{case}"))
    draws = rng.multinomial(
        len(streams), np.full(len(streams), 1.0 / len(streams)),
        size=BOOTSTRAP_REPLICATES,
    )
    boot = (draws @ sums) / (draws @ counts)[:, None]
    rows = []
    for index, observable in enumerate(("J12", "J23", "J12_plus_J23")):
        low, high = np.quantile(boot[:, index], (0.025, 0.975))
        se = float(np.std(boot[:, index], ddof=1))
        rows.append({
            "case": case,
            "scope": "complete_held_out_trajectory",
            "source": "direct_trajectory",
            "observable": observable,
            "estimate": float(point[index]),
            "ci_low": float(low),
            "ci_high": float(high),
            "bootstrap_se": se,
            "z_from_zero": float(point[index] / se),
            "streams": len(streams),
            "rows": int(counts.sum()),
            "reference_zero_applicable": observable == "J12_plus_J23",
        })
    return rows


def main() -> None:
    full_rows = direct_full_rows("driven") + direct_full_rows("equilibrium")
    write_csv(ANALYSIS / "current_balance_full_trajectory_audit.csv", full_rows)

    with (ANALYSIS / "current_balance.csv").open() as handle:
        masked = list(csv.DictReader(handle))
    with (ANALYSIS / "density_trajectory_comparisons.csv").open() as handle:
        comparisons = list(csv.DictReader(handle))

    masked_comparisons = []
    for case in ("driven", "equilibrium"):
        comparison = next(
            row for row in comparisons
            if row["case"] == case and row["diagnostic"] == "current_balance"
        )
        density = next(
            row for row in masked
            if row["case"] == case and row["source"] == "new_ensemble"
            and row["observable"] == "J12_plus_J23"
        )
        trajectory = next(
            row for row in masked
            if row["case"] == case and row["source"] == "direct_trajectory"
            and row["observable"] == "J12_plus_J23"
        )
        z = float(comparison["difference_z_approx"])
        masked_comparisons.append({
            "case": case,
            "scope": "section4_1_support_mask",
            "density_estimate": density["estimate"],
            "density_ci_low": density["ci_low"],
            "density_ci_high": density["ci_high"],
            "trajectory_estimate": trajectory["estimate"],
            "trajectory_ci_low": trajectory["ci_low"],
            "trajectory_ci_high": trajectory["ci_high"],
            "density_minus_trajectory": comparison["difference"],
            "difference_z": z,
            "density_trajectory_gate_pass": abs(z) <= 1.96,
            "reference_zero_applicable_after_conditioning": False,
        })
    write_csv(ANALYSIS / "current_balance_masked_comparison_audit.csv", masked_comparisons)

    outcome = {
        "frozen_masked_zero_gate_retained": False,
        "scope_correction": (
            "E[J12+J23]=0 is tested only on the unconditioned held-out trajectories; "
            "the support-conditioned density is tested only against the identically "
            "conditioned held-out trajectory estimate."
        ),
        "global_direct_zero_gate": {
            row["case"]: bool(float(row["ci_low"]) <= 0.0 <= float(row["ci_high"]))
            for row in full_rows if row["observable"] == "J12_plus_J23"
        },
        "masked_density_trajectory_gate": {
            row["case"]: bool(row["density_trajectory_gate_pass"])
            for row in masked_comparisons
        },
        "density_evaluated_outside_support": False,
    }
    (ANALYSIS / "current_scope_audit_verdict.json").write_text(
        json.dumps(outcome, indent=2) + "\n"
    )
    print(json.dumps(outcome, indent=2))


if __name__ == "__main__":
    main()
