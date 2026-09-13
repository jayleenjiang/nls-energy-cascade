#!/usr/bin/env python3
"""Generate the prospectively frozen seed/argument matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
POPULATIONS = (512, 1024, 2048, 4096)
OPTIONAL_POPULATION = 8192
RUNS = range(1, 5)


def seed_n2(k: float, population_index: int, run: int) -> int:
    return 210_000_000 + int(round(10 * k)) * 10_000 + population_index * 100 + run


def seed_long(n: int, k: float, population_index: int, run: int) -> int:
    return 310_000_000 + n * 100_000 + population_index * 1_000 + int(round(10 * k)) * 10 + run


def main() -> None:
    fields = [
        "study", "stage", "n", "clone_count", "burnin", "horizon",
        "selection_time", "dt", "k", "gauge_shift", "control_scale",
        "seed", "threads", "run", "prefix",
    ]
    rows: list[dict[str, object]] = []
    for p_index, population in enumerate((*POPULATIONS, OPTIONAL_POPULATION), 1):
        stage = "primary" if population in POPULATIONS else "conditional_8192"
        for k in (0.3, 0.5, 0.7):
            for run in RUNS:
                rows.append({
                    "study": "n2_known_answer",
                    "stage": stage,
                    "n": 2,
                    "clone_count": population,
                    "burnin": 500,
                    "horizon": 0.1,
                    "selection_time": 0.1,
                    "dt": 0.0005,
                    "k": k,
                    "gauge_shift": 0.1,
                    "control_scale": 0.5,
                    "seed": seed_n2(k, p_index, run),
                    "threads": 5,
                    "run": run,
                    "prefix": f"raw/n2_known_answer/N{population}/n2_k{str(k).replace('.', 'p')}_run{run}",
                })

    specifications = (
        (10, (0.3, 0.7), 60),
        (20, (0.4, 0.6), 60),
        (30, (0.4, 0.6), 60),
        (40, (0.4, 0.6), 120),
    )
    for n, tilts, horizon in specifications:
        for p_index, population in enumerate((*POPULATIONS, OPTIONAL_POPULATION), 1):
            stage = "primary" if population in POPULATIONS else "conditional_8192"
            for k in tilts:
                for run in RUNS:
                    rows.append({
                        "study": "long_chain",
                        "stage": stage,
                        "n": n,
                        "clone_count": population,
                        "burnin": 500,
                        "horizon": horizon,
                        "selection_time": 2,
                        "dt": 0.0005,
                        "k": k,
                        "gauge_shift": 0.1,
                        "control_scale": 0.5,
                        "seed": seed_long(n, k, p_index, run),
                        "threads": 5,
                        "run": run,
                        "prefix": f"raw/long_chain/n{n}/N{population}/n{n}_k{str(k).replace('.', 'p')}_run{run}",
                    })

    seeds = [int(row["seed"]) for row in rows]
    if len(seeds) != len(set(seeds)):
        raise RuntimeError("seed collision")
    with (ROOT / "RUN_MATRIX.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
