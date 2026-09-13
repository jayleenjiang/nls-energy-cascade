#!/usr/bin/env python3
"""Independent integrity and acceptance audit for entropy-cloning pilots."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}")
    return rows[0]


def cloning_rows(root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(root.rglob("*_summary.csv")):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "clone_count" not in reader.fieldnames:
                continue
            data = list(reader)
        if len(data) != 1:
            raise ValueError(f"expected one row in {path}")
        row = data[0]
        rows.append(
            {
                "path": str(path),
                "mode": row.get("mode", "naive"),
                "n": int(row["n"]),
                "clone_count": int(row["clone_count"]),
                "observation_time": float(row["observation_time"]),
                "selection_time": float(row["selection_time"]),
                "dt": float(row["dt"]),
                "k": float(row["k"]),
                "scgf": float(row["scgf"]),
                "entropy_rate": float(row["mean_population_entropy_rate"]),
                "action_current": float(row["mean_population_action_current"]),
                "minimum_weight_ess": float(row["minimum_weight_ess"]),
                "minimum_root_ess": float(row["minimum_root_count_ess"]),
                "minimum_unique_roots": int(row["minimum_unique_roots"]),
                "midpoint_failures": int(row["midpoint_failures"]),
            }
        )
    return rows


def mean_se(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    mean = float(array.mean())
    se = float(array.std(ddof=1) / math.sqrt(array.size)) if array.size > 1 else math.nan
    return mean, se


def direct_lookup(path: Path) -> dict[tuple[int, float, float], dict[str, str]]:
    with path.open(newline="") as handle:
        return {
            (int(row["n"]), float(row["tau"]), round(float(row["k"]), 12)): row
            for row in csv.DictReader(handle)
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment_dir", type=Path)
    parser.add_argument("production_summary", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    root = args.experiment_dir
    rows = cloning_rows(root)
    if not rows:
        raise SystemExit("no cloning summary rows found")
    production = read_row(args.production_summary)
    direct = direct_lookup(root / "direct_pilot" / "scgf_points.csv")

    checks: list[dict[str, object]] = []

    def add(name: str, status: str, detail: str) -> None:
        checks.append({"check": name, "status": status, "detail": detail})

    midpoint_total = sum(int(row["midpoint_failures"]) for row in rows)
    add(
        "successful_summary_midpoint_integrity",
        "PASS" if midpoint_total == 0 else "FAIL",
        f"{midpoint_total} failures over {len(rows)} successful summaries",
    )

    numerical_rejections: list[str] = []
    for path in sorted(root.rglob("*.log")):
        contents = path.read_text(errors="replace")
        if "ERROR:" in contents and "numerical gate" in contents:
            numerical_rejections.append(str(path.relative_to(root)))
    add(
        "numerical_gate_rejections",
        "PASS" if not numerical_rejections else "BLOCKED",
        "none" if not numerical_rejections else (
            f"{len(numerical_rejections)} attempted run(s) rejected before summary: "
            + ", ".join(numerical_rejections)
        ),
    )

    direct_audit = (root / "direct_pilot" / "audit.md").read_text()
    add(
        "direct_scgf_audit",
        "PASS" if "Status: **PASS**" in direct_audit else "FAIL",
        "independent direct-SCGF audit status",
    )

    k0 = [
        row for row in rows
        if "/k0_recovery/" in str(row["path"]) and float(row["k"]) == 0.0
    ]
    entropy_mean, entropy_se = mean_se([float(row["entropy_rate"]) for row in k0])
    action_mean, action_se = mean_se([float(row["action_current"]) for row in k0])
    entropy_reference = float(production["mean_entropy_rate"])
    action_reference = float(production["mean_action_current"])
    entropy_z = abs(entropy_mean - entropy_reference) / entropy_se
    action_z = abs(action_mean - action_reference) / action_se
    add(
        "k0_unbiased_recovery",
        "PASS" if len(k0) >= 4 and entropy_z <= 3.0 and action_z <= 3.0 else "FAIL",
        f"runs={len(k0)}, entropy z={entropy_z:.3f}, action-current z={action_z:.3f}",
    )

    identity = [row for row in rows if "/generator_identity/" in str(row["path"])]
    identity_groups: defaultdict[float, dict[str, float]] = defaultdict(dict)
    for row in identity:
        identity_groups[float(row["k"])][str(row["mode"])] = float(row["scgf"])
    identity_differences = {
        k: abs(values["naive"] - values["guided"])
        for k, values in identity_groups.items()
        if {"naive", "guided"}.issubset(values)
    }
    identity_ok = bool(identity_differences) and all(
        difference <= 0.02 for difference in identity_differences.values()
    )
    add(
        "one_interval_generator_identity",
        "PASS" if identity_ok else "FAIL",
        ", ".join(f"k={k:g}: |difference|={value:.5f}"
                  for k, value in sorted(identity_differences.items())),
    )

    groups: defaultdict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (
            row["mode"], row["n"], row["clone_count"],
            row["observation_time"], row["selection_time"], row["dt"], row["k"],
        )
        groups[key].append(row)

    accepted_groups: list[dict[str, object]] = []
    for key, group in sorted(groups.items()):
        mode, n, clone_count, observation, selection, dt, k = key
        mean, se = mean_se([float(row["scgf"]) for row in group])
        support = (
            len(group) >= 4
            and min(float(row["minimum_weight_ess"]) for row in group)
            >= 0.1 * int(clone_count)
            and min(int(row["minimum_unique_roots"]) for row in group) >= 32
            and min(float(row["minimum_root_ess"]) for row in group) >= 16.0
            and sum(int(row["midpoint_failures"]) for row in group) == 0
        )
        direct_row = direct.get((int(n), float(observation), round(float(k), 12)))
        direct_consistent: bool | None = None
        direct_z = math.nan
        if direct_row and int(direct_row["direct_reliability_gate"]) == 1 and math.isfinite(se) and se > 0:
            direct_z = abs(mean - float(direct_row["psi"])) / se
            direct_consistent = direct_z <= 3.0
        accepted_groups.append(
            {
                "mode": mode,
                "n": n,
                "clone_count": clone_count,
                "observation_time": observation,
                "selection_time": selection,
                "dt": dt,
                "k": k,
                "runs": len(group),
                "mean_scgf": mean,
                "run_scgf_se": se,
                "support_gate": support,
                "direct_consistent": direct_consistent,
                "direct_z": direct_z,
            }
        )

    low_tilt = [
        row for row in accepted_groups
        if row["mode"] == "guided" and row["n"] == 10
        and row["observation_time"] == 20.0 and row["clone_count"] >= 512
        and row["selection_time"] == 0.5 and row["dt"] == 0.0005
        and float(row["k"]) <= 0.3 and row["direct_consistent"] is not None
    ]
    add(
        "low_tilt_direct_crosscheck",
        "PASS" if low_tilt and all(bool(row["direct_consistent"]) for row in low_tilt) else "FAIL",
        "; ".join(
            f"N={row['clone_count']}, k={row['k']}: z={float(row['direct_z']):.2f}"
            for row in low_tilt
        ),
    )

    population_series = sorted(
        (
            row for row in accepted_groups
            if row["mode"] == "guided" and row["n"] == 10
            and row["observation_time"] == 20.0
            and row["selection_time"] == 0.5 and row["dt"] == 0.0005
            and abs(float(row["k"]) - 0.3) < 1.0e-12
            and int(row["runs"]) >= 4
        ),
        key=lambda row: int(row["clone_count"]),
    )
    population_comparisons: list[str] = []
    population_converged = False
    if len(population_series) >= 2:
        lower, upper = population_series[-2:]
        difference = abs(float(upper["mean_scgf"]) - float(lower["mean_scgf"]))
        combined_se = math.hypot(
            float(lower["run_scgf_se"]), float(upper["run_scgf_se"])
        )
        population_converged = (
            difference <= 2.0 * combined_se and difference <= 0.02
            and bool(upper["support_gate"])
        )
        population_comparisons.append(
            f"N={lower['clone_count']}->{upper['clone_count']}: "
            f"difference={difference:.5f}, 2SE={2.0 * combined_se:.5f}, "
            f"upper support={upper['support_gate']}"
        )
    add(
        "population_convergence_k0.3",
        "PASS" if population_converged else "BLOCKED",
        "; ".join(population_comparisons) or "fewer than two populations",
    )

    selection_rows = {
        float(row["selection_time"]): row
        for row in accepted_groups
        if row["mode"] == "guided" and row["n"] == 10
        and row["clone_count"] == 1024 and row["observation_time"] == 20.0
        and row["dt"] == 0.0005 and abs(float(row["k"]) - 0.3) < 1.0e-12
        and int(row["runs"]) >= 4
    }
    selection_details: list[str] = []
    selection_pass = all(value in selection_rows for value in (0.25, 0.5, 1.0))
    if selection_pass:
        reference = selection_rows[0.5]
        for value in (0.25, 1.0):
            candidate = selection_rows[value]
            difference = abs(
                float(candidate["mean_scgf"]) - float(reference["mean_scgf"])
            )
            combined_se = math.hypot(
                float(candidate["run_scgf_se"]),
                float(reference["run_scgf_se"]),
            )
            comparison_pass = (
                difference <= 2.0 * combined_se and difference <= 0.01
                and bool(candidate["support_gate"])
            )
            selection_pass = selection_pass and comparison_pass
            selection_details.append(
                f"Delta={value:g} vs 0.5: difference={difference:.5f}, "
                f"2SE={2.0 * combined_se:.5f}, support={candidate['support_gate']}"
            )
    add(
        "selection_interval_convergence_k0.3",
        "PASS" if selection_pass else "BLOCKED",
        "; ".join(selection_details) or "requires Delta=0.25,0.5,1.0",
    )

    timestep_rows = {
        float(row["dt"]): row
        for row in accepted_groups
        if row["mode"] == "guided" and row["n"] == 10
        and row["clone_count"] == 1024 and row["observation_time"] == 20.0
        and row["selection_time"] == 0.5
        and abs(float(row["k"]) - 0.3) < 1.0e-12
        and int(row["runs"]) >= 4
    }
    timestep_pass = all(value in timestep_rows for value in (0.00025, 0.0005))
    timestep_detail = "requires dt=0.00025 and 0.0005"
    if timestep_pass:
        fine = timestep_rows[0.00025]
        coarse = timestep_rows[0.0005]
        difference = abs(float(fine["mean_scgf"]) - float(coarse["mean_scgf"]))
        combined_se = math.hypot(
            float(fine["run_scgf_se"]), float(coarse["run_scgf_se"])
        )
        timestep_pass = (
            difference <= 2.0 * combined_se and difference <= 0.01
            and bool(fine["support_gate"])
        )
        timestep_detail = (
            f"difference={difference:.5f}, 2SE={2.0 * combined_se:.5f}, "
            f"fine support={fine['support_gate']}"
        )
    add(
        "timestep_convergence_k0.3",
        "PASS" if timestep_pass else "BLOCKED",
        timestep_detail,
    )

    observation_rows = {
        float(row["observation_time"]): row
        for row in accepted_groups
        if row["mode"] == "guided" and row["n"] == 10
        and row["clone_count"] == 1024 and row["selection_time"] == 0.5
        and row["dt"] == 0.0005 and abs(float(row["k"]) - 0.3) < 1.0e-12
        and int(row["runs"]) >= 4
    }
    observation_pass = 40.0 in observation_rows and bool(
        observation_rows[40.0]["support_gate"]
    )
    if 40.0 in observation_rows:
        extended = observation_rows[40.0]
        observation_detail = (
            f"t=40 support={extended['support_gate']}, "
            f"mean={float(extended['mean_scgf']):.5f} +/- "
            f"{float(extended['run_scgf_se']):.5f}"
        )
    else:
        observation_detail = "requires four t=40 runs"
    add(
        "observation_time_extension_k0.3",
        "PASS" if observation_pass else "BLOCKED",
        observation_detail,
    )

    high_tilt_supported = [
        row for row in accepted_groups
        if row["mode"] == "guided" and row["n"] == 10
        and row["observation_time"] >= 20.0 and float(row["k"]) >= 0.5
        and bool(row["support_gate"])
    ]
    add(
        "high_tilt_sampling_support",
        "PASS" if high_tilt_supported else "BLOCKED",
        "no k>=0.5 group passes support" if not high_tilt_supported
        else f"{len(high_tilt_supported)} supported group(s)",
    )

    standard_groups = {
        round(float(row["k"]), 12): row
        for row in accepted_groups
        if row["mode"] == "guided" and row["n"] == 10
        and row["clone_count"] == 1024 and row["observation_time"] == 20.0
        and row["selection_time"] == 0.5 and row["dt"] == 0.0005
        and int(row["runs"]) >= 4 and bool(row["support_gate"])
    }
    reliable_pairs: list[tuple[float, float]] = []
    for k in sorted(standard_groups):
        partner = round(1.0 - k, 12)
        if k < partner and partner in standard_groups:
            reliable_pairs.append((k, partner))
    pair_details: list[str] = []
    pair_consistent = bool(reliable_pairs)
    for k, partner in reliable_pairs:
        left = standard_groups[k]
        right = standard_groups[partner]
        residual = float(left["mean_scgf"]) - float(right["mean_scgf"])
        combined_se = math.hypot(
            float(left["run_scgf_se"]), float(right["run_scgf_se"])
        )
        consistent = abs(residual) <= 2.0 * combined_se
        pair_consistent = pair_consistent and consistent
        pair_details.append(
            f"k={k:g}<->{partner:g}: residual={residual:.5f}, "
            f"2SE={2.0 * combined_se:.5f}, consistent={consistent}"
        )
    add(
        "gc_pair_ready",
        "PASS" if pair_consistent else "BLOCKED",
        "; ".join(pair_details) or "0 standard-setting supported pairs",
    )

    integrity_failures = [
        item for item in checks
        if item["status"] == "FAIL" and item["check"] not in {"low_tilt_direct_crosscheck"}
    ]
    if integrity_failures:
        overall = "FAIL"
    elif pair_consistent:
        overall = "PASS_WITH_SUPPORTED_GC_PAIR"
    else:
        overall = "PASS_WITH_FT_BLOCKED"
    payload = {
        "overall": overall,
        "successful_cloning_summaries": len(rows),
        "checks": checks,
        "accepted_group_diagnostics": accepted_groups,
        "reliable_gc_pairs": reliable_pairs,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = [
        "# Entropy cloning audit",
        "",
        f"Overall: **{overall}**",
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    lines.extend(
        f"| {item['check']} | {item['status']} | {item['detail']} |"
        for item in checks
    )
    lines.extend(["", (
        "`PASS_WITH_FT_BLOCKED` means that the implemented diagnostics pass "
        "their integrity checks, but no fluctuation-theorem claim is accepted "
        "because a nontrivial `k <-> 1-k` pair lacks reliable support."
    ), ""])
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output} and {args.output.with_suffix('.json')}: {overall}")


if __name__ == "__main__":
    main()
