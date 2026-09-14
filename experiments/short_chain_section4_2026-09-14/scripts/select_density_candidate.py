#!/usr/bin/env python3
"""Apply the frozen equilibrium-validation model-selection rule."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
ROOT = HERE / "stationary_density" / "equilibrium_candidates"
OUT = HERE / "stationary_density" / "FROZEN_DENSITY_CHOICE.json"
SEEDS = (4101, 4102, 4103)


def main() -> None:
    records = []
    for path in sorted(ROOT.glob("*/validation_metrics.json")):
        if "_failed_serialization_" in path.parent.name:
            continue
        record = json.loads(path.read_text())
        record["metrics_path"] = str(path.relative_to(HERE))
        records.append(record)
    if not records:
        raise RuntimeError(f"no candidate metrics under {ROOT}")

    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for record in records:
        key = (
            int(record["components"]),
            tuple(int(x) for x in record["widths"]),
            float(record["lambda_fp"]),
        )
        grouped[key].append(record)

    expected = {
        (8, (64, 64), 0.01), (8, (64, 64), 0.1),
        (16, (64, 64), 0.01), (16, (64, 64), 0.1),
        (16, (128, 128), 0.01), (16, (128, 128), 0.1),
    }
    if set(grouped) != expected:
        raise RuntimeError(f"candidate matrix mismatch: {sorted(grouped)}")

    summaries = []
    for key in sorted(grouped):
        items = sorted(grouped[key], key=lambda x: int(x["seed"]))
        if tuple(int(x["seed"]) for x in items) != SEEDS:
            raise RuntimeError(f"seed matrix mismatch for {key}")
        nll = np.asarray([float(x["best_validation_nll"]) for x in items])
        admissible = all(bool(x["selection_admissible"]) for x in items)
        summaries.append({
            "components": key[0],
            "widths": list(key[1]),
            "lambda_fp": key[2],
            "all_seed_admissible": admissible,
            "median_validation_nll": float(np.median(nll)),
            "mean_validation_nll": float(np.mean(nll)),
            "across_seed_se_nll": float(np.std(nll, ddof=1)/math.sqrt(len(nll))),
            "seeds": items,
        })

    eligible = [x for x in summaries if x["all_seed_admissible"]]
    result = {
        "protocol_version": "section4-v1",
        "selection_data": "equilibrium_dt2p5e-4 validation streams only",
        "candidate_summaries": summaries,
        "status": "NO_ADMISSIBLE_CANDIDATE",
        "selected": None,
    }
    if eligible:
        numerical_best = min(eligible, key=lambda x: x["median_validation_nll"])
        tied = [
            x for x in eligible
            if x["median_validation_nll"] <= numerical_best["median_validation_nll"]
            + max(x["across_seed_se_nll"], numerical_best["across_seed_se_nll"])
        ]
        selected = min(
            tied,
            key=lambda x: (
                x["components"], tuple(x["widths"]), x["lambda_fp"],
                x["median_validation_nll"],
            ),
        )
        result.update({
            "status": "FROZEN",
            "numerical_best": {
                k: numerical_best[k]
                for k in ("components", "widths", "lambda_fp",
                          "median_validation_nll", "across_seed_se_nll")
            },
            "tie_set": [
                {k: x[k] for k in ("components", "widths", "lambda_fp",
                                    "median_validation_nll", "across_seed_se_nll")}
                for x in tied
            ],
            "selected": {
                k: selected[k]
                for k in ("components", "widths", "lambda_fp",
                          "median_validation_nll", "across_seed_se_nll")
            },
        })

    OUT.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["selected"] is None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
