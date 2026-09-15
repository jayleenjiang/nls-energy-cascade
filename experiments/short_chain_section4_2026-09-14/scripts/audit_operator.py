#!/usr/bin/env python3
"""Independent Cartesian-to-reduced generator and Gibbs-null audit."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE / "stationary_density"))

from reduced_operator import (  # noqa: E402
    diffusion_covariance,
    drift,
    energy,
    gibbs_adjoint_relative_residual,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cartesian_force(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    action = x*x + y*y
    sqr = x*x - y*y
    sqi = 2*x*y
    mass = action.sum(axis=-1)
    fr = np.zeros_like(x)
    fi = np.zeros_like(y)
    for j in range(3):
        nr = np.zeros_like(mass)
        ni = np.zeros_like(mass)
        if j > 0:
            nr += sqr[..., j-1]
            ni += sqi[..., j-1]
        if j < 2:
            nr += sqr[..., j+1]
            ni += sqi[..., j+1]
        onsite = 2*mass - action[..., j]
        fr[..., j] = onsite*x[..., j] + 2*(nr*x[..., j] + ni*y[..., j])
        fi[..., j] = onsite*y[..., j] + 2*(ni*x[..., j] - nr*y[..., j])
    return fr, fi


def cartesian_induced_coefficients(state: np.ndarray, t1: float, t3: float,
                                   gamma: float) -> tuple[np.ndarray, np.ndarray]:
    i1, i2, i3, th1, th3 = np.moveaxis(state, -1, 0)
    phi = np.stack((0.5*th1, np.zeros_like(th1), 0.5*th3), axis=-1)
    action = np.stack((i1, i2, i3), axis=-1)
    radius = np.sqrt(action)
    x = radius*np.cos(phi)
    y = radius*np.sin(phi)
    fr, fi = cartesian_force(x, y)
    ax = -fi
    ay = fr
    ax[..., 0] -= gamma*fr[..., 0]
    ay[..., 0] -= gamma*fi[..., 0]
    ax[..., 2] -= gamma*fr[..., 2]
    ay[..., 2] -= gamma*fi[..., 2]

    i_drift = 2*x*ax + 2*y*ay
    i_drift[..., 0] += 4*gamma*t1
    i_drift[..., 2] += 4*gamma*t3
    phi_drift = (-y*ax+x*ay)/action
    reduced_drift = np.stack((
        i_drift[..., 0], i_drift[..., 1], i_drift[..., 2],
        2*(phi_drift[..., 0]-phi_drift[..., 1]),
        2*(phi_drift[..., 2]-phi_drift[..., 1]),
    ), axis=-1)

    # Jacobian of reduced variables with respect to bath coordinates
    # (x1,y1,x3,y3), multiplied by the Cartesian bath covariance.
    jac = np.zeros(state.shape[:-1] + (5, 4), dtype=np.float64)
    jac[..., 0, 0] = 2*x[..., 0]
    jac[..., 0, 1] = 2*y[..., 0]
    jac[..., 2, 2] = 2*x[..., 2]
    jac[..., 2, 3] = 2*y[..., 2]
    jac[..., 3, 0] = -2*y[..., 0]/i1
    jac[..., 3, 1] = 2*x[..., 0]/i1
    jac[..., 4, 2] = -2*y[..., 2]/i3
    jac[..., 4, 3] = 2*x[..., 2]/i3
    cart_cov = np.diag((2*gamma*t1, 2*gamma*t1, 2*gamma*t3, 2*gamma*t3))
    reduced_cov = np.einsum("...ik,kl,...jl->...ij", jac, cart_cov, jac)
    return reduced_drift, reduced_cov


def main() -> None:
    rng = np.random.default_rng(2026091401)
    n = 20_000
    states = np.empty((n, 5), dtype=np.float64)
    states[:, :3] = np.exp(rng.uniform(np.log(1e-4), np.log(12.0), size=(n, 3)))
    states[:, 3:] = rng.uniform(-np.pi, np.pi, size=(n, 2))

    comparisons = {}
    for label, t1, t3 in (("equilibrium", 5.0, 5.0), ("driven", 2.0, 8.0)):
        induced_b, induced_cov = cartesian_induced_coefficients(states, t1, t3, 0.1)
        analytic_b = drift(states, t1, t3, 0.1)
        analytic_cov = diffusion_covariance(states, t1, t3, 0.1)
        offdiagonal = induced_cov.copy()
        diagonal_index = np.arange(5)
        offdiagonal[..., diagonal_index, diagonal_index] = 0.0
        comparisons[label] = {
            "max_abs_drift_error": float(np.max(np.abs(induced_b-analytic_b))),
            "max_rel_drift_error": float(np.max(np.abs(induced_b-analytic_b)/(1+np.abs(analytic_b)))),
            "max_abs_covariance_error": float(np.max(np.abs(induced_cov-analytic_cov))),
            "max_offdiagonal_covariance": float(np.max(np.abs(offdiagonal))),
        }

    gibbs_residual = gibbs_adjoint_relative_residual(states, 5.0, 0.1)
    shifted = states.copy()
    shifted[:, 3] += 2*np.pi
    periodic_energy_error = float(np.max(np.abs(energy(shifted)-energy(states))))
    shifted[:, 4] -= 4*np.pi
    periodic_energy_error = max(periodic_energy_error,
                                float(np.max(np.abs(energy(shifted)-energy(states)))))

    # At I1=0 or I3=0, b_I=4 gamma T and d_I(4 gamma T I rho)=4 gamma T rho
    # for finite rho, so the normal stationary probability current vanishes.
    boundary_current_coefficients = {
        "I1_b_minus_da_at_zero": float(4*0.1*5.0-4*0.1*5.0),
        "I3_b_minus_da_at_zero": float(4*0.1*5.0-4*0.1*5.0),
        "I2_drift_at_zero": 0.0,
    }

    source = REPO / "experiments/spectral_gap_controlled_2026-09-07/source/NLS_stationary_autocorr_controlled.cpp"
    notebook = REPO / "KDE/4:15_NN/FKE_5d_NLS.ipynb"
    result = {
        "protocol_version": "section4-v1",
        "points": n,
        "state_action_range": [float(states[:, :3].min()), float(states[:, :3].max())],
        "comparisons": comparisons,
        "gibbs_relative_residual": {
            "max_abs": float(np.max(np.abs(gibbs_residual))),
            "rms": float(np.sqrt(np.mean(gibbs_residual*gibbs_residual))),
        },
        "periodic_energy_max_abs_error": periodic_energy_error,
        "boundary_current_coefficients": boundary_current_coefficients,
        "source_sha256": sha256(source),
        "legacy_notebook_sha256": sha256(notebook),
    }
    tolerance = 2e-10
    result["pass"] = bool(
        all(v["max_abs_drift_error"] <= tolerance and
            v["max_abs_covariance_error"] <= tolerance and
            v["max_offdiagonal_covariance"] <= tolerance
            for v in comparisons.values())
        and result["gibbs_relative_residual"]["max_abs"] <= tolerance
        and periodic_energy_error <= tolerance
        and all(abs(x) <= tolerance for x in boundary_current_coefficients.values())
    )
    out = HERE / "stationary_density" / "OPERATOR_AUDIT.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
