#!/usr/bin/env python3
"""Audited reduced Ito generator for the three-mode Cartesian NLS dynamics."""

from __future__ import annotations

import numpy as np


def energy(state: np.ndarray) -> np.ndarray:
    """Physical energy E=H/2 for (...,5) reduced states."""
    state = np.asarray(state, dtype=np.float64)
    i1, i2, i3, th1, th3 = np.moveaxis(state, -1, 0)
    mass = i1 + i2 + i3
    return (0.5 * mass * mass - 0.25 * (i1*i1 + i2*i2 + i3*i3)
            + i1*i2*np.cos(th1) + i2*i3*np.cos(th3))


def drift(state: np.ndarray, t1: float, t3: float, gamma: float = 0.1) -> np.ndarray:
    """Reduced Ito drift in (I1,I2,I3,theta1,theta3)."""
    state = np.asarray(state, dtype=np.float64)
    i1, i2, i3, th1, th3 = np.moveaxis(state, -1, 0)
    s1, c1 = np.sin(th1), np.cos(th1)
    s3, c3 = np.sin(th3), np.cos(th3)
    b1 = (4*i1*i2*s1 - 4*gamma*i1*i2*(1+c1) + 4*gamma*t1
          - 2*gamma*(i1*i1 + 2*i1*i3))
    b2 = -4*i2*(i1*s1 + i3*s3)
    b3 = (4*i3*i2*s3 - 4*gamma*i3*i2*(1+c3) + 4*gamma*t3
          - 2*gamma*(i3*i3 + 2*i1*i3))
    bt1 = 2*(i2-i1) + 4*(i2-i1)*c1 + 4*gamma*i2*s1 - 4*i3*c3
    bt3 = 2*(i2-i3) + 4*(i2-i3)*c3 + 4*gamma*i2*s3 - 4*i1*c1
    return np.stack((b1, b2, b3, bt1, bt3), axis=-1)


def diffusion_covariance(state: np.ndarray, t1: float, t3: float,
                         gamma: float = 0.1) -> np.ndarray:
    """Ito covariance B B^T of the reduced state increments per unit time."""
    state = np.asarray(state, dtype=np.float64)
    i1, _, i3, _, _ = np.moveaxis(state, -1, 0)
    out = np.zeros(state.shape[:-1] + (5, 5), dtype=np.float64)
    out[..., 0, 0] = 8*gamma*t1*i1
    out[..., 2, 2] = 8*gamma*t3*i3
    out[..., 3, 3] = 8*gamma*t1/i1
    out[..., 4, 4] = 8*gamma*t3/i3
    return out


def drift_divergence(state: np.ndarray, t1: float, t3: float,
                     gamma: float = 0.1) -> np.ndarray:
    """Coordinate divergence of the reduced drift."""
    del t1, t3
    state = np.asarray(state, dtype=np.float64)
    i1, i2, i3, th1, th3 = np.moveaxis(state, -1, 0)
    s1, c1 = np.sin(th1), np.cos(th1)
    s3, c3 = np.sin(th3), np.cos(th3)
    db1 = 4*i2*s1 - 4*gamma*i2*(1+c1) - 4*gamma*(i1+i3)
    db2 = -4*(i1*s1+i3*s3)
    db3 = 4*i2*s3 - 4*gamma*i2*(1+c3) - 4*gamma*(i3+i1)
    db4 = -4*(i2-i1)*s1 + 4*gamma*i2*c1
    db5 = -4*(i2-i3)*s3 + 4*gamma*i2*c3
    return db1 + db2 + db3 + db4 + db5


def gibbs_log_derivatives(state: np.ndarray, temperature: float) -> tuple[np.ndarray, np.ndarray]:
    """Gradient and diagonal Hessian of log rho_eq=-E/T+constant."""
    state = np.asarray(state, dtype=np.float64)
    i1, i2, i3, th1, th3 = np.moveaxis(state, -1, 0)
    mass = i1 + i2 + i3
    s1, c1 = np.sin(th1), np.cos(th1)
    s3, c3 = np.sin(th3), np.cos(th3)
    e_grad = np.stack((
        mass - 0.5*i1 + i2*c1,
        mass - 0.5*i2 + i1*c1 + i3*c3,
        mass - 0.5*i3 + i2*c3,
        -i1*i2*s1,
        -i2*i3*s3,
    ), axis=-1)
    e_hdiag = np.stack((
        np.full_like(i1, 0.5),
        np.full_like(i2, 0.5),
        np.full_like(i3, 0.5),
        -i1*i2*c1,
        -i2*i3*c3,
    ), axis=-1)
    return -e_grad/temperature, -e_hdiag/temperature


def adjoint_relative_residual_from_log_derivatives(
    state: np.ndarray,
    log_grad: np.ndarray,
    log_hdiag: np.ndarray,
    t1: float,
    t3: float,
    gamma: float = 0.1,
) -> np.ndarray:
    """Return (L^dagger rho)/rho from derivatives of log rho."""
    state = np.asarray(state, dtype=np.float64)
    log_grad = np.asarray(log_grad, dtype=np.float64)
    log_hdiag = np.asarray(log_hdiag, dtype=np.float64)
    b = drift(state, t1, t3, gamma)
    div_b = drift_divergence(state, t1, t3, gamma)
    i1, _, i3, _, _ = np.moveaxis(state, -1, 0)
    a = np.stack((
        4*gamma*t1*i1,
        np.zeros_like(i1),
        4*gamma*t3*i3,
        4*gamma*t1/i1,
        4*gamma*t3/i3,
    ), axis=-1)
    diffusion = np.sum(a*(log_hdiag + log_grad*log_grad), axis=-1)
    diffusion += 8*gamma*t1*log_grad[..., 0]
    diffusion += 8*gamma*t3*log_grad[..., 2]
    return -div_b - np.sum(b*log_grad, axis=-1) + diffusion


def gibbs_adjoint_relative_residual(state: np.ndarray, temperature: float,
                                    gamma: float = 0.1) -> np.ndarray:
    grad, hdiag = gibbs_log_derivatives(state, temperature)
    return adjoint_relative_residual_from_log_derivatives(
        state, grad, hdiag, temperature, temperature, gamma
    )


WEAK_FUNCTION_NAMES = (
    "I1", "I2", "I3", "I1_sq", "I2_sq", "I3_sq",
    "I1_I2", "I2_I3", "I1_I3",
    "sin_theta1", "cos_theta1", "sin_theta3", "cos_theta3",
    "I1_sin_theta1", "I2_sin_theta1", "I3_sin_theta3", "I2_sin_theta3",
    "I1_I2_sin_theta1", "I1_I2_cos_theta1",
    "I3_I2_sin_theta3", "I3_I2_cos_theta3",
    "M", "M_sq", "E",
)


def weak_generator_values(state: np.ndarray, t1: float, t3: float,
                          gamma: float = 0.1) -> np.ndarray:
    """Evaluate Lf for the frozen weak-identity family."""
    state = np.asarray(state, dtype=np.float64)
    i1, i2, i3, th1, th3 = np.moveaxis(state, -1, 0)
    s1, c1 = np.sin(th1), np.cos(th1)
    s3, c3 = np.sin(th3), np.cos(th3)
    b = drift(state, t1, t3, gamma)
    b1, b2, b3, bt1, bt3 = np.moveaxis(b, -1, 0)
    a1 = 4*gamma*t1*i1
    a3 = 4*gamma*t3*i3
    at1 = 4*gamma*t1/i1
    at3 = 4*gamma*t3/i3

    lsin1 = bt1*c1-at1*s1
    lcos1 = -bt1*s1-at1*c1
    lsin3 = bt3*c3-at3*s3
    lcos3 = -bt3*s3-at3*c3
    li1i2 = i2*b1+i1*b2
    li2i3 = i3*b2+i2*b3
    li1i3 = i3*b1+i1*b3
    mass = i1+i2+i3
    bmass = b1+b2+b3

    e_grad, e_hdiag = gibbs_log_derivatives(state, -1.0)
    # gibbs_log_derivatives(T=-1) returns +grad(E), +diag Hess(E).
    lenergy = np.sum(b*e_grad, axis=-1)
    lenergy += a1*e_hdiag[..., 0] + a3*e_hdiag[..., 2]
    lenergy += at1*e_hdiag[..., 3] + at3*e_hdiag[..., 4]

    return np.stack((
        b1, b2, b3,
        2*i1*b1+2*a1, 2*i2*b2, 2*i3*b3+2*a3,
        li1i2, li2i3, li1i3,
        lsin1, lcos1, lsin3, lcos3,
        s1*b1+i1*lsin1,
        s1*b2+i2*lsin1,
        s3*b3+i3*lsin3,
        s3*b2+i2*lsin3,
        s1*li1i2+i1*i2*lsin1,
        c1*li1i2+i1*i2*lcos1,
        s3*li2i3+i3*i2*lsin3,
        c3*li2i3+i3*i2*lcos3,
        bmass,
        2*mass*bmass+2*(a1+a3),
        lenergy,
    ), axis=-1)
