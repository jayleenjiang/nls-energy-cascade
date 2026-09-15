"""
Generate Latin Hypercube initial points for NLS backward solver.
5D: (I1, I2, I3, theta1, theta3)
I_j in [0.1, I_max], theta_j in [-pi, pi]
"""
import numpy as np
from scipy.stats.qmc import LatinHypercube

N_train = 64000 # must be multiple of 16
I_max = 10.0
pi = np.pi

sampler = LatinHypercube(d=5, seed=42)
lhs_raw = sampler.random(n=N_train).astype(np.float32)

X_train = np.column_stack([
    (lhs_raw[:, 0] * (I_max - 0.1) + 0.1),   # I1 in [0.1, I_max]
    (lhs_raw[:, 1] * (I_max - 0.1) + 0.1),   # I2
    (lhs_raw[:, 2] * (I_max - 0.1) + 0.1),   # I3
    (lhs_raw[:, 3] * 2 * pi - pi),             # theta1 in [-pi, pi]
    (lhs_raw[:, 4] * 2 * pi - pi),             # theta3 in [-pi, pi]
])

np.savetxt('NLS_backward_LHS_X_train.txt', X_train, fmt='%.8f')
print(f"Generated {N_train} LHS points in 5D")
print(f"I range: [{X_train[:,:3].min():.2f}, {X_train[:,:3].max():.2f}]")
print(f"θ range: [{X_train[:,3:].min():.4f}, {X_train[:,3:].max():.4f}]")
print(f"Saved to NLS_backward_LHS_X_train.txt")
