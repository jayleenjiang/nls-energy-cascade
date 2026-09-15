import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.01
tt = np.arange(Y.shape[1]) * dt_eff

def model_osc(t, A, lamR, lamI, phi, c):
    return A * np.exp(lamR * t) * np.sin(lamI * t + phi) + c

lamR_all = []
lamI_all = []
for i in range(len(Y)):
    try:
        popt, _ = curve_fit(model_osc, tt, Y[i,:],
                           p0=[0.5, -1.0, 1.5, 0.0, Y[i,-50:].mean()],
                           maxfev=10000,
                           bounds=([-5, -5, 0, -np.pi, -2], [5, 0, 5, np.pi, 2]))
        if abs(popt[1]) < 4.9 and popt[2] < 4.9:
            lamR_all.append(popt[1])
            lamI_all.append(popt[2])
    except:
        pass

lamR_all = np.array(lamR_all)
lamI_all = np.array(lamI_all)
print(f"Good fits: {len(lamR_all)}/{len(Y)}")
print(f"λ_R = {lamR_all.mean():.3f} ± {lamR_all.std():.3f}")
print(f"λ_I = {lamI_all.mean():.3f} ± {lamI_all.std():.3f}")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Scatter: λ_R vs λ_I
axes[0].scatter(lamR_all, lamI_all, s=1, alpha=0.3)
axes[0].set_xlabel('λ_R'); axes[0].set_ylabel('λ_I')
axes[0].set_title('Eigenvalue estimates (each point = one initial condition)')

# Histogram: λ_R
axes[1].hist(lamR_all, bins=80, edgecolor='white')
axes[1].set_xlabel('λ_R'); axes[1].set_title(f'λ_R distribution (mean={lamR_all.mean():.2f})')

# Histogram: λ_I
axes[2].hist(lamI_all, bins=80, edgecolor='white')
axes[2].set_xlabel('λ_I'); axes[2].set_title(f'λ_I distribution (mean={lamI_all.mean():.2f})')

plt.tight_layout(); plt.savefig('eigenvalue_scatter.png', dpi=150); plt.show()
