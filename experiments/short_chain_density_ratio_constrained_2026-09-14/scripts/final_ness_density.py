#!/usr/bin/env python3
"""Evaluate the final absolutely normalized five-dimensional NESS density."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp

from double_calibrated_transport_ratio_model import load_bundle as load_fine_bundle
from triple_calibrated_transport_ratio_model import load_bundle as load_coarse_bundle

ROOT = Path(__file__).resolve().parent.parent


def energy(state):
    state = np.asarray(state, np.float64)
    i1, i2, i3, th1, th3 = state.T
    mass = i1 + i2 + i3
    return (0.5 * mass * mass - 0.25 * (i1 * i1 + i2 * i2 + i3 * i3)
            + i1 * i2 * np.cos(th1) + i2 * i3 * np.cos(th3))


class FinalNESSDensity:
    """Arithmetic ensemble of three independently fitted normalized densities."""

    def __init__(self, timestep="fine"):
        if timestep not in ("fine", "coarse"):
            raise ValueError(timestep)
        tf.keras.backend.set_floatx("float64")
        if timestep == "fine":
            paths = [ROOT / f"recovery_r9/factor_0.625/models/seed_{s}"
                     for s in (8201, 8202, 8203)]
            loader = load_fine_bundle
        else:
            paths = [ROOT / f"recovery_r10_coarse/models/seed_{s}"
                     for s in (8201, 8202, 8203)]
            loader = load_coarse_bundle
        self.models = [loader(path) for path in paths]
        with (ROOT / "final_analysis/normalizer_summary.csv").open(newline="") as handle:
            rows = [row for row in csv.DictReader(handle) if row["timestep"] == timestep]
        rows.sort(key=lambda row: int(row["model"]))
        self.log_ratio_normalizers = np.asarray(
            [float(row["log_ratio_normalizer"]) for row in rows]
        )
        self.log_gibbs_partition = float(rows[0]["gibbs_log_partition"])

    def component_log_density(self, state, batch=8192):
        state = np.asarray(state, np.float64)
        if state.ndim != 2 or state.shape[1] != 5:
            raise ValueError("state must have shape (N,5)")
        if not np.isfinite(state).all() or np.any(state[:, :3] <= 0):
            raise ValueError("states must be finite with strictly positive actions")
        values = []
        base = -energy(state) / 5.0 - self.log_gibbs_partition
        for model, logz in zip(self.models, self.log_ratio_normalizers):
            ratio = np.concatenate([
                model.log_ratio(tf.constant(state[start:start + batch], tf.float64)).numpy()
                for start in range(0, len(state), batch)
            ])
            values.append(base + ratio - logz)
        return np.asarray(values).T

    def log_density(self, state, batch=8192):
        components = self.component_log_density(state, batch)
        return logsumexp(components, axis=1) - math.log(components.shape[1])

    def density(self, state, batch=8192):
        return np.exp(self.log_density(state, batch))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="CSV with I1,I2,I3,theta1,theta3")
    parser.add_argument("output", type=Path)
    parser.add_argument("--timestep", choices=("fine", "coarse"), default="fine")
    args = parser.parse_args()
    data = np.genfromtxt(args.input, delimiter=",", names=True)
    names = ("I1", "I2", "I3", "theta1", "theta3")
    state = np.column_stack([data[name] for name in names])
    model = FinalNESSDensity(args.timestep)
    logp = model.log_density(state)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle); writer.writerow((*names, "log_rho_ss", "rho_ss"))
        writer.writerows(np.column_stack((state, logp, np.exp(logp))))


if __name__ == "__main__":
    main()
