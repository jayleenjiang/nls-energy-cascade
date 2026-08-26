import numpy as np
import matplotlib.pyplot as plt

Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.001
tt = np.arange(Y.shape[1]) * dt_eff

# Find initial points with strongest signal (furthest from steady state)
ss = Y[:, -10:].mean(axis=1)
signal_strength = np.abs(Y[:, 0] - ss)
top_idx = np.argsort(signal_strength)[-20:]  # top 20 strongest

fig, axes = plt.subplots(2, 1, figsize=(14, 8))
for i in top_idx:
    signal = Y[i,:] - ss[i]
    axes[0].plot(tt, signal, linewidth=1)
    # Normalize to see shape
    axes[1].plot(tt, signal / signal[0], linewidth=1)

axes[0].set_xlabel('time'); axes[0].set_ylabel('signal')
axes[0].set_title('Top 20 strongest raw signals')
axes[0].axhline(0, color='k', linewidth=0.5)

axes[1].set_xlabel('time'); axes[1].set_ylabel('signal / signal(0)')
axes[1].set_title('Normalized')
axes[1].axhline(0, color='k', linewidth=0.5)
axes[1].set_ylim(-0.5, 1.2)

plt.tight_layout(); plt.show()