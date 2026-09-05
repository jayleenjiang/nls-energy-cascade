#!/usr/bin/env python3
"""Symbolic checks for n2_ft_theorem.tex.

The script verifies only model-specific algebra.  The recurrence, path-space,
and positive-operator arguments remain mathematical proofs, not computer
algebra claims.
"""

import sympy as sp


def main() -> None:
    x1, y1, x2, y2 = sp.symbols("x1 y1 x2 y2", real=True)
    variables = (x1, y1, x2, y2)
    i1 = x1**2 + y1**2
    i2 = x2**2 + y2**2
    mass = i1 + i2
    c1 = x1 + sp.I * y1
    c2 = x2 + sp.I * y2
    energy = sp.expand(
        sp.Rational(1, 2) * mass**2
        - sp.Rational(1, 4) * (i1**2 + i2**2)
        + sp.re(c2**2 * sp.conjugate(c1) ** 2)
    )

    euler = sum(value * sp.diff(energy, value) for value in variables)
    assert sp.simplify(euler - 4 * energy) == 0

    for pair in ((x1, y1), (x2, y2)):
        laplacian = sum(sp.diff(energy, value, 2) for value in pair)
        assert sp.simplify(laplacian - 4 * mass) == 0

    # Coefficient of |F_r|^2 in L_k exp(aE) / exp(aE).
    a, k, temperature = sp.symbols("a k T", real=True)
    coefficient = (
        temperature * a**2
        + (2 * k - 1) * a
        + (k**2 - k) / temperature
    )
    factorized = (temperature * a + k) * (
        a - (1 - k) / temperature
    )
    assert sp.simplify(coefficient - factorized) == 0

    # The un-tilted exponential Lyapunov coefficient is the k=0 case.
    assert sp.simplify(coefficient.subs(k, 0) + a * (1 - a * temperature)) == 0

    print("PASS: n=2 homogeneity, Laplacians, and tilted Lyapunov factorization")


if __name__ == "__main__":
    main()
