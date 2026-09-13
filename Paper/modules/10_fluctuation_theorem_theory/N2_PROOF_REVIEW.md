# Internal adversarial review: two-site FT theorem

Date: 2026-08-30

## Verdict

The `n=2` working proof closes the analytic assumptions that remained open in
the arbitrary finite-chain ledger for interior tilts `0<k<1`.  No numerical
fit is used in the proof.  The result is suitable for advisor review, but it
should not be inserted into the paper as an established theorem until a
probability/functional-analysis reader checks the cited theorem mappings and
operator-core details.

## Dependency audit

| Step | Evidence | Internal status |
|---|---|---|
| quartic coercivity | positive action-angle decomposition | closed |
| gradient coercivity | Euler identity plus Cauchy--Schwarz | closed |
| nonexplosion and NESS | exponential Foster--Lyapunov drift | closed |
| uniqueness and positive density | uniform ellipticity, irreducibility, strong Feller | closed |
| finite-time total-entropy FT | finite Fisher information plus diffusion time reversal | closed via cited theorem |
| entropy tilted generator | direct Ito product calculation | closed |
| principal SCGF | multiplicative Lyapunov, majorizer, positive kernel, minorization | closed via cited theorem |
| common adjoint domain | coercive `L^2` form and compact embedding | closed at working-note level |
| GC symmetry | Hilbert adjoint plus unitary time reversal | closed for `0<k<1` |

## Stress tests and exclusions

1. The proof cannot be copied to `n>=3`: uniform ellipticity is lost.
2. Endpoint moments `k=0,1` are excluded; the Lyapunov interval degenerates at
   `k=1` and needs a separate domain argument.
3. The finite-time identity is forward--reverse for total entropy, not a
   same-forward histogram identity.
4. The long-time theorem concerns medium entropy, not action current.
5. The theorem is for the continuous Cartesian SDE, not the finite-step
   splitting integrator.
6. No conclusion is drawn about uniformity in `n` or `n -> infinity`.

## Items an external reader should check

- the cutoff passage in the stationary Fisher-information identity;
- the precise finite-entropy hypotheses in the Haussmann--Pardoux/Chetrite--
  Gawedzki time-reversal theorem;
- that the stochastic-integral Feynman--Kac kernel satisfies the exact
  regularity formulation used by Ferre--Rousset--Stoltz;
- the core/closure statement for the non-selfadjoint `L^2` form.

These are requests for theorem-level verification, not known counterexamples
or currently open algebraic gaps.
