#!/usr/bin/env python3
"""Apply the frozen BC1 endpoint-reproduction gate without changing data."""

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REF = {
    25: (0.737907, 0.165409),
    50: (0.514213, 0.105599),
    100: (0.361686, 0.0728212),
}


def rows(path):
    with path.open() as handle:
        clean = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(clean))


def summary(path):
    with path.open() as handle:
        return next(csv.DictReader(handle))


def pooled(parts):
    # Recover each replicate's sample variance from SE and trajectory count,
    # then pool independent trajectory samples exactly from sufficient stats.
    total_n = sum(n for n, _, _ in parts)
    mean = sum(n * m for n, m, _ in parts) / total_n
    ss = 0.0
    for n, m, se in parts:
        variance = se * se * n
        ss += (n - 1) * variance + n * (m - mean) ** 2
    variance = ss / (total_n - 1)
    return mean, math.sqrt(variance / total_n), total_n


def evaluate(n):
    profile_rows = []
    nonfinite = 0
    discarded = 0
    seeds = []
    for rep in (0, 1):
        prefix = ROOT / "raw" / f"BC1_T10_T2_n{n}_rep{rep}"
        pfile = Path(str(prefix) + "_profile.csv")
        sfile = Path(str(prefix) + "_summary.csv")
        if not pfile.exists() or not sfile.exists():
            raise FileNotFoundError(f"missing completed gate output for {prefix}")
        p = rows(pfile)
        s = summary(sfile)
        if int(s["valid_trajectories"]) <= 1:
            raise RuntimeError(f"insufficient valid trajectories in {sfile}")
        profile_rows.append((p, int(s["valid_trajectories"])))
        nonfinite += int(s["nonfinite_trajectories"])
        discarded += int(s["discarded_trajectories"])
        seeds.append(int(s["seed"]))

    checks = []
    for side, index, reference in (("left", 0, REF[n][0]),
                                   ("right", n - 1, REF[n][1])):
        parts = []
        for p, count in profile_rows:
            row = p[index]
            parts.append((count, float(row["mean_I"]), float(row["se_mean_I"])))
        mean, se, count = pooled(parts)
        discrepancy = mean - reference
        tolerance = max(3.0 * se, 0.02 * abs(reference), 0.002)
        checks.append({
            "side": side,
            "reference": reference,
            "new_mean": mean,
            "new_se": se,
            "trajectory_count": count,
            "signed_discrepancy": discrepancy,
            "absolute_discrepancy": abs(discrepancy),
            "tolerance": tolerance,
            "pass": abs(discrepancy) <= tolerance,
        })
    passed = all(item["pass"] for item in checks) and nonfinite == 0 and discarded == 0
    result = {
        "n": n,
        "seeds": seeds,
        "nonfinite_trajectories": nonfinite,
        "discarded_trajectories": discarded,
        "endpoint_checks": checks,
        "pass": passed,
    }
    out = ROOT / "analysis" / f"bc1_gate_n{n}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return passed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=sorted(REF), required=True)
    args = parser.parse_args()
    raise SystemExit(0 if evaluate(args.n) else 2)


if __name__ == "__main__":
    main()

