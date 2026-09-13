#!/usr/bin/env python3
"""Freeze the Phase-III production grid using support diagnostics only."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PILOT = ROOT / "raw" / "pilot"
EXPECTED = {
    10: [(0.25, 0.75), (0.35, 0.65), (0.45, 0.55)],
    20: [(0.25, 0.75), (0.35, 0.65), (0.45, 0.55)],
    30: [(0.25, 0.75), (0.35, 0.65), (0.45, 0.55)],
    40: [(0.25, 0.75), (0.35, 0.65), (0.45, 0.55)],
}
SETTINGS = {
    10: (2048, 60, 0.3, 0.7),
    20: (1024, 60, 0.4, 0.6),
    30: (4096, 60, 0.4, 0.6),
    40: (1024, 120, 0.4, 0.6),
}


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"expected one summary row in {path}")
    return rows[0]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def member_support(row: dict[str, str], n: int, k: float) -> tuple[bool, str]:
    checks = {
        "mode": row.get("mode", "").startswith("controlled_exact"),
        "n": int(row["n"]) == n,
        "clone_count": int(row["clone_count"]) == 512,
        "horizon": close(float(row["observation_time"]), 20.0),
        "selection": close(float(row["selection_time"]), 2.0),
        "dt": close(float(row["dt"]), 0.0005),
        "k": close(float(row["k"]), k),
        "gauge": close(float(row["gauge_shift"]), 0.1),
        "control": close(float(row["control_scale"]), 0.5),
        "midpoint": int(row["midpoint_failures"]) == 0,
        "weight_ess": float(row["minimum_weight_ess"]) >= 51.2,
        "unique_roots": int(row["minimum_unique_roots"]) >= 32,
        "root_count_ess": float(row["minimum_root_count_ess"]) >= 16.0,
        "root_weight_ess": float(row["minimum_root_weight_ess"]) >= 16.0,
    }
    failed = [name for name, passed in checks.items() if not passed]
    return not failed, "pass" if not failed else ";".join(failed)


def main() -> None:
    if (ROOT / "PRODUCTION_STARTED").exists():
        raise SystemExit("production already started; refusing to rewrite grid")
    summaries: dict[tuple[int, float], tuple[Path, dict[str, str]]] = {}
    for path in sorted(PILOT.rglob("*_summary.csv")):
        row = read_one(path)
        key = (int(row["n"]), round(float(row["k"]), 12))
        if key in summaries:
            raise ValueError(f"duplicate pilot member {key}")
        summaries[key] = (path, row)

    support_rows: list[dict[str, object]] = []
    eligible_by_n: dict[int, list[tuple[float, float]]] = {}
    for n, pairs in EXPECTED.items():
        eligible_by_n[n] = []
        for low, high in pairs:
            member_results: list[tuple[bool, str, Path, dict[str, str]]] = []
            for k in (low, high):
                found = summaries.get((n, round(k, 12)))
                if found is None:
                    raise FileNotFoundError(f"missing pilot summary n={n}, k={k}")
                path, row = found
                timeseries = path.with_name(
                    path.name.replace("_summary.csv", "_timeseries.csv")
                )
                if not timeseries.is_file() or timeseries.stat().st_size == 0:
                    raise FileNotFoundError(f"missing timeseries for {path}")
                passed, reason = member_support(row, n, k)
                member_results.append((passed, reason, path, row))
                support_rows.append(
                    {
                        "n": n,
                        "k": k,
                        "support_pass": int(passed),
                        "reason": reason,
                        "minimum_weight_ess": row["minimum_weight_ess"],
                        "minimum_unique_roots": row["minimum_unique_roots"],
                        "minimum_root_count_ess": row["minimum_root_count_ess"],
                        "minimum_root_weight_ess": row["minimum_root_weight_ess"],
                        "midpoint_failures": row["midpoint_failures"],
                        "summary_path": str(path.relative_to(ROOT)),
                    }
                )
            if all(result[0] for result in member_results):
                eligible_by_n[n].append((low, high))

    write_csv(ROOT / "PILOT_SUPPORT.csv", support_rows)
    blocked = {n: pairs for n, pairs in eligible_by_n.items() if len(pairs) < 2}
    if blocked:
        detail = ", ".join(f"n={n}: {len(pairs)} pairs" for n, pairs in blocked.items())
        (ROOT / "PILOT_VERDICT.md").write_text(
            "# Pilot verdict\n\nStatus: **BLOCKED**\n\n" + detail + "\n"
        )
        raise SystemExit(f"support pilot blocked: {detail}")

    frozen: list[dict[str, object]] = []
    verdict_lines = [
        "# Pilot verdict",
        "",
        "Status: **PASS**",
        "",
        "Selection used only the predeclared support fields; SCGF values were",
        "not read by this program.",
        "",
        "| n | eligible pairs | selected outer | selected inner |",
        "|---:|:---|:---:|:---:|",
    ]
    for n, eligible in eligible_by_n.items():
        ordered = sorted(eligible, key=lambda pair: pair[0])
        selected = [("outer", ordered[0]), ("inner", ordered[-1])]
        clones, horizon, baseline_low, baseline_high = SETTINGS[n]
        verdict_lines.append(
            f"| {n} | {', '.join(f'({a:g},{b:g})' for a,b in ordered)} "
            f"| ({selected[0][1][0]:g},{selected[0][1][1]:g}) "
            f"| ({selected[1][1][0]:g},{selected[1][1][1]:g}) |"
        )
        for role_index, (role, (low, high)) in enumerate(selected):
            frozen.append(
                {
                    "n": n,
                    "role": role,
                    "k_low": low,
                    "k_high": high,
                    "clone_count": clones,
                    "horizon": horizon,
                    "baseline_k_low": baseline_low,
                    "baseline_k_high": baseline_high,
                    "production_seed_base": 960000 + n * 100 + role_index * 20,
                    "population_seed_base": 970000 + n * 100 + role_index * 20,
                }
            )
    write_csv(ROOT / "FROZEN_GRID.csv", frozen)
    digest = hashlib.sha256((ROOT / "FROZEN_GRID.csv").read_bytes()).hexdigest()
    verdict_lines.extend(["", f"`FROZEN_GRID.csv` SHA-256: `{digest}`", ""])
    (ROOT / "PILOT_VERDICT.md").write_text("\n".join(verdict_lines))
    print(f"frozen grid: {len(frozen)} rows; sha256={digest}")


if __name__ == "__main__":
    main()

