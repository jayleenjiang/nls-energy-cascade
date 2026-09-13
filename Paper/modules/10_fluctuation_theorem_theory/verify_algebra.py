#!/usr/bin/env python3
"""Symbolic checks supporting finite_n_ft_proof.tex.

This script is not the proof.  It independently checks the model-specific
polynomial identities for n=3 and the scalar coefficient algebra used in the
tilted-generator adjoint calculation.
"""

import sympy as sp


def main() -> None:
    xs = sp.symbols("x1:4", real=True)
    ys = sp.symbols("y1:4", real=True)
    variables = [value for pair in zip(xs, ys) for value in pair]
    cs = [sp.expand_complex(x + sp.I * y) for x, y in zip(xs, ys)]
    actions = [sp.expand(x * x + y * y) for x, y in zip(xs, ys)]
    mass = sum(actions)

    interaction = 0
    for j in range(1, 3):
        interaction += sp.re(cs[j] ** 2 * sp.conjugate(cs[j - 1]) ** 2)
    energy = sp.expand(
        sp.Rational(1, 2) * mass**2
        - sp.Rational(1, 4) * sum(value**2 for value in actions)
        + interaction
    )

    reversed_energy = sp.expand(energy.subs(
        {value: -value for value in ys}, simultaneous=True
    ))
    assert sp.simplify(reversed_energy - energy) == 0

    for site in (0, 2):
        boundary_laplacian = (
            sp.diff(energy, xs[site], 2) + sp.diff(energy, ys[site], 2)
        )
        assert sp.simplify(boundary_laplacian - 4 * mass) == 0

    # The code force is grad E.  Check the explicit formula at each site.
    square_real = [x * x - y * y for x, y in zip(xs, ys)]
    square_imag = [2 * x * y for x, y in zip(xs, ys)]
    for j in range(3):
        neighbor_real = sum(
            square_real[m] for m in (j - 1, j + 1) if 0 <= m < 3
        )
        neighbor_imag = sum(
            square_imag[m] for m in (j - 1, j + 1) if 0 <= m < 3
        )
        onsite = 2 * mass - actions[j]
        force_x = onsite * xs[j] + 2 * (
            neighbor_real * xs[j] + neighbor_imag * ys[j]
        )
        force_y = onsite * ys[j] + 2 * (
            neighbor_imag * xs[j] - neighbor_real * ys[j]
        )
        assert sp.simplify(sp.diff(energy, xs[j]) - force_x) == 0
        assert sp.simplify(sp.diff(energy, ys[j]) - force_y) == 0

    # Adjoint at k: drift coefficient alpha=gamma(2k-1).  Taking the
    # Lebesgue adjoint contributes -alpha*Delta E.  The resulting coefficients
    # must equal those of the tilt 1-k.
    k, gamma = sp.symbols("k gamma", real=True)
    alpha = gamma * (2 * k - 1)
    adjoint_drift = -alpha
    target_drift = gamma * (2 * (1 - k) - 1)
    assert sp.simplify(adjoint_drift - target_drift) == 0

    adjoint_laplacian_potential = gamma * k - alpha
    target_laplacian_potential = gamma * (1 - k)
    assert sp.simplify(
        adjoint_laplacian_potential - target_laplacian_potential
    ) == 0
    assert sp.simplify(k**2 - k - ((1 - k) ** 2 - (1 - k))) == 0

    print("PASS: energy, force, time reversal, boundary Laplacian, and tilt algebra")


if __name__ == "__main__":
    main()
