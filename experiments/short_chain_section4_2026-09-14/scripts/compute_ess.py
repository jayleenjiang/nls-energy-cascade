#!/usr/bin/env python3
"""Compute predeclared trajectory-level integrated autocorrelation times."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
SPLITS = HERE / "DIRECT_STATE_SPLITS.csv"
OUT_CSV = HERE / "analysis" / "effective_sample_size.csv"
OUT_JSON = HERE / "analysis" / "effective_sample_size_gate.json"
N_SNAP = 100_001
N_COL = 5

PROBES = (
    "I1", "I2", "I3", "log_I1", "log_I2", "log_I3",
    "cos_theta1", "sin_theta1", "cos_theta3", "sin_theta3",
    "I1_minus_I3", "I1_plus_I3",
)


def load_probes(path: Path) -> np.ndarray:
    raw = np.fromfile(path, dtype="<f8").reshape(N_SNAP, N_COL)
    i1, i2, i3, theta1, theta3 = raw.T
    c1, s1 = np.cos(theta1), np.sin(theta1)
    c3, s3 = np.cos(theta3), np.sin(theta3)
    idiff = i1-i3
    isum = i1+i3
    return np.column_stack((
        i1, i2, i3, np.log(i1), np.log(i2), np.log(i3),
        c1, s1, c3, s3, idiff, isum,
    ))


def ips_tau(values: np.ndarray) -> np.ndarray:
    """Initial-positive-sequence IAT for every column."""
    n = values.shape[0]
    centered = values-values.mean(axis=0, keepdims=True)
    nfft = 1 << (2*n-1).bit_length()
    spectrum = np.fft.rfft(centered, n=nfft, axis=0)
    acov = np.fft.irfft(spectrum.conj()*spectrum, n=nfft, axis=0)[:n]
    acov /= np.arange(n, 0, -1, dtype=np.float64)[:, None]
    variance = acov[0]
    if np.any(variance <= 0) or not np.isfinite(variance).all():
        raise RuntimeError("nonpositive or nonfinite probe variance")
    rho = acov/variance
    tau = np.ones(values.shape[1], dtype=np.float64)
    active = np.ones(values.shape[1], dtype=bool)
    # Geyer's initial-positive sequence of adjacent autocorrelation pairs.
    for lag in range(1, n-1, 2):
        pair = rho[lag]+rho[lag+1]
        keep = active & np.isfinite(pair) & (pair > 0)
        tau[keep] += 2*pair[keep]
        active &= keep
        if not active.any():
            break
    return np.maximum(tau, 1.0)


def main() -> None:
    if not SPLITS.exists():
        raise FileNotFoundError("run audit_inputs.py first")
    with SPLITS.open(newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    rows = [row for row in all_rows if row["branch"] == "density"]
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["case"], row["split"])].append(row)

    output_rows = []
    gate_records = []
    for (case, split), members in sorted(grouped.items()):
        per_stream_tau = []
        for row in sorted(members, key=lambda item: int(item["stream_id"])):
            values = load_probes(REPO/row["relative_path"])
            per_stream_tau.append(ips_tau(values))
        taus = np.asarray(per_stream_tau)
        total_n = len(members)*N_SNAP
        combined_tau = np.mean(taus, axis=0)
        combined_ess = total_n/combined_tau
        for index, probe in enumerate(PROBES):
            output_rows.append({
                "case": case,
                "split": split,
                "probe": probe,
                "streams": len(members),
                "raw_snapshots": total_n,
                "mean_iat_snapshots": combined_tau[index],
                "median_stream_iat": np.median(taus[:, index]),
                "max_stream_iat": np.max(taus[:, index]),
                "effective_samples": combined_ess[index],
            })
        minimum_index = int(np.argmin(combined_ess))
        primary_fine_test = split == "test" and case in (
            "driven_dt2p5e-4", "equilibrium_dt2p5e-4"
        )
        gate_records.append({
            "case": case,
            "split": split,
            "streams": len(members),
            "raw_snapshots": total_n,
            "minimum_effective_samples": float(combined_ess[minimum_index]),
            "limiting_probe": PROBES[minimum_index],
            "primary_fine_test": primary_fine_test,
            "gate_threshold": 5000 if primary_fine_test else None,
            "gate_pass": bool(combined_ess[minimum_index] >= 5000) if primary_fine_test else None,
        })

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    gate = {
        "protocol_version": "section4-v1",
        "method": "per-stream Geyer initial-positive-sequence IAT; stationary-time weighted combination",
        "records": gate_records,
        "primary_gate_pass": all(
            item["gate_pass"] for item in gate_records if item["primary_fine_test"]
        ),
    }
    OUT_JSON.write_text(json.dumps(gate, indent=2, sort_keys=True)+"\n")
    print(json.dumps(gate, indent=2, sort_keys=True))
    if not gate["primary_gate_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
