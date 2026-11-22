"""4.3.3: phase average of the I2 growth rate at FROZEN endpoint actions.

With I2 = 0 and (I1, I3) held fixed, the two bond angles (advisor-note convention th1 = 2(phi1-phi2),
th3 = 2(phi3-phi2)) form a closed diffusion on the torus:
  canonical bath, paper normalisation:
     dth1 = [-2 I1 (1 + 2 cos th1) - 4 I3 cos th3] dt + sqrt(8 gamma T1 / I1) dW1
     dth3 = [-2 I3 (1 + 2 cos th3) - 4 I1 cos th1] dt + sqrt(8 gamma T3 / I3) dW3
  HLNS reduced dynamics (their normalisation, angle noise sqrt(gamma) dB as on the advisor's slide):
     dth1 = -[I1 (1 + 2 cos th1) + 2 I3 cos th3] dt + sqrt(gamma) dW1
     dth3 = -[I3 (1 + 2 cos th3) + 2 I1 cos th1] dt + sqrt(gamma) dW3
F(I1, I3) = E[I1 sin th1 + I3 sin th3] under the stationary law of this diffusion (the advisor's quantity).
The stationary Fokker--Planck equation is solved by finite volumes with Scharfetter--Gummel fluxes
(positivity-preserving, exact for constant drift), N x N cells, sparse direct solve.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

TWO_PI = 2 * np.pi


def bern(x):
    x = np.asarray(x, float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-8
    out[small] = 1 - x[small] / 2
    xs = x[~small]
    out[~small] = xs / np.expm1(xs)
    return out


def drift(model, I1, I3):
    if model == "canonical":
        return (lambda a, b: -2 * I1 * (1 + 2 * np.cos(a)) - 4 * I3 * np.cos(b),
                lambda a, b: -2 * I3 * (1 + 2 * np.cos(b)) - 4 * I1 * np.cos(a))
    if model == "hlns":
        return (lambda a, b: -(I1 * (1 + 2 * np.cos(a)) + 2 * I3 * np.cos(b)),
                lambda a, b: -(I3 * (1 + 2 * np.cos(b)) + 2 * I1 * np.cos(a)))
    raise ValueError(model)


def half_var(model, I1, I3, gamma, T1, T3):
    """d = sigma^2 / 2 for th1, th3"""
    if model == "canonical":
        return 4 * gamma * T1 / I1, 4 * gamma * T3 / I3
    return gamma / 2, gamma / 2


def solve(model, I1, I3, gamma=0.1, T1=2.0, T3=8.0, N=256):
    h = TWO_PI / N
    c = -np.pi + (np.arange(N) + 0.5) * h          # cell centres
    f = -np.pi + (np.arange(N) + 1.0) * h          # face to the right of cell i
    b1f, b3f = drift(model, I1, I3)
    d1, d3 = half_var(model, I1, I3, gamma, T1, T3)
    idx = np.arange(N * N).reshape(N, N)           # [i (th1), j (th3)]
    rows, cols, vals = [], [], []

    def add_faces(bface, d, L, R):
        # flux from cell L to cell R: J = (d/h) [B(-Pe) p_L - B(Pe) p_R], Pe = b h / d
        Pe = bface * h / d
        a, bb = (d / h) * bern(-Pe), (d / h) * bern(Pe)
        L, R, a, bb = L.ravel(), R.ravel(), a.ravel(), bb.ravel()
        # cell L loses J/h, cell R gains J/h
        rows.extend([L, L, R, R]); cols.extend([L, R, L, R])
        vals.extend([-a / h, bb / h, a / h, -bb / h])

    A1, A3 = np.meshgrid(f, c, indexing="ij")      # faces in th1 direction: th1 = f_i, th3 = c_j
    add_faces(b1f(A1, A3), d1, idx, np.roll(idx, -1, axis=0))
    B1, B3 = np.meshgrid(c, f, indexing="ij")      # faces in th3 direction
    add_faces(b3f(B1, B3), d3, idx, np.roll(idx, -1, axis=1))
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N * N, N * N))
    A = A.tolil(); A[0, :] = h * h; A = A.tocsr()
    rhs = np.zeros(N * N); rhs[0] = 1.0
    p = spla.spsolve(A, rhs).reshape(N, N)
    T1g, T3g = np.meshgrid(c, c, indexing="ij")
    w = p * h * h
    out = {"F": float(np.sum(w * (I1 * np.sin(T1g) + I3 * np.sin(T3g)))),
           "E_sin1": float(np.sum(w * np.sin(T1g))), "E_sin3": float(np.sum(w * np.sin(T3g))),
           "pmin": float(p.min()), "mass": float(w.sum())}
    return out, p, c


def one_d_exact_sin(model, I3, gamma, T3, M=4000):
    """I1 -> 0 limit: th3 alone, dth = b(th) dt + sqrt(2 d) dW on the circle.  Exact stationary density
    (Risken):  p(th) ∝ int_0^{2pi} exp(Phi(th) - Phi(th + s)) ds,  Phi(th) = int_0^th b/d
    (from J = b p - d p' = const and periodicity).
    Evaluated in log space (logsumexp) to avoid overflow."""
    from scipy.special import logsumexp
    if model == "canonical":
        b = lambda t: -2 * I3 * (1 + 2 * np.cos(t)); d = 4 * gamma * T3 / I3
    else:
        b = lambda t: -I3 * (1 + 2 * np.cos(t)); d = gamma / 2
    ds = TWO_PI / M
    t = np.arange(2 * M) * ds
    # exact antiderivative of b/d: b = -k (1 + 2 cos t)  ->  int_0^t b/d = -(k/d)(t + 2 sin t)
    k = 2 * I3 if model == "canonical" else I3
    Phi = -(k / d) * (t + 2 * np.sin(t))
    th = t[:M]
    lp = np.empty(M)
    for a in range(0, M, 200):
        kk = np.arange(a, min(a + 200, M))
        S = Phi[kk][:, None] - Phi[kk[:, None] + np.arange(M)[None, :]]
        lp[kk] = logsumexp(S, axis=1) + np.log(ds)
    lp -= lp.max()
    p = np.exp(lp); p /= p.sum() * ds
    return float(np.sum(p * np.sin(th)) * ds)
