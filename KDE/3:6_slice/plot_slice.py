"""
Plot 2D slice test results: conditional PDF f(theta1, theta3 | I1=I2=I3=2)
Run after NLS5D_slice produces output files.
"""
import numpy as np
import matplotlib.pyplot as plt

G_bin = 50
G_kde = 30

# ============================================================
# Test 1: Simple binning 50x50
# ============================================================
data1 = np.loadtxt("slice_test1_binning50.txt")
th1_1 = data1[:, 0].reshape(G_bin, G_bin)
th3_1 = data1[:, 1].reshape(G_bin, G_bin)
den_1 = data1[:, 2].reshape(G_bin, G_bin)

# ============================================================
# Test 2: KDE results
# ============================================================
files = {
    "local_h0.05":  "slice_test2_kde_local_h0.05.txt",
    "local_h0.10":  "slice_test2_kde_local_h0.10.txt",
    "global_h0.05": "slice_test2_kde_global_h0.05.txt",
    "global_h0.10": "slice_test2_kde_global_h0.10.txt",
}

kde_data = {}
for key, fname in files.items():
    d = np.loadtxt(fname)
    kde_data[key] = {
        "th1": d[:, 0].reshape(G_kde, G_kde),
        "th3": d[:, 1].reshape(G_kde, G_kde),
        "den": d[:, 2].reshape(G_kde, G_kde),
    }

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle(r"Conditional PDF $f(\theta_1, \theta_3 \mid I_1=I_2=I_3=2)$", fontsize=16)

# (1) Simple binning 50x50
ax = axes[0, 0]
im = ax.pcolormesh(th1_1, th3_1, den_1, cmap="viridis", shading="auto")
ax.set_title("Test 1: Simple binning 50x50")
ax.set_xlabel(r"$\theta_1$")
ax.set_ylabel(r"$\theta_3$")
plt.colorbar(im, ax=ax, label="density")

# (2a) KDE local, h=0.05
ax = axes[0, 1]
d = kde_data["local_h0.05"]
im = ax.pcolormesh(d["th1"], d["th3"], d["den"], cmap="viridis", shading="auto")
ax.set_title("Test 2a: KDE local bin, h=0.05")
ax.set_xlabel(r"$\theta_1$")
ax.set_ylabel(r"$\theta_3$")
plt.colorbar(im, ax=ax, label="density")

# (2a) KDE local, h=0.10
ax = axes[0, 2]
d = kde_data["local_h0.10"]
im = ax.pcolormesh(d["th1"], d["th3"], d["den"], cmap="viridis", shading="auto")
ax.set_title("Test 2a: KDE local bin, h=0.10")
ax.set_xlabel(r"$\theta_1$")
ax.set_ylabel(r"$\theta_3$")
plt.colorbar(im, ax=ax, label="density")

# (2b) KDE global, h=0.05
ax = axes[1, 0]
d = kde_data["global_h0.05"]
im = ax.pcolormesh(d["th1"], d["th3"], d["den"], cmap="viridis", shading="auto")
ax.set_title("Test 2b: KDE global, h=0.05")
ax.set_xlabel(r"$\theta_1$")
ax.set_ylabel(r"$\theta_3$")
plt.colorbar(im, ax=ax, label="density")

# (2b) KDE global, h=0.10
ax = axes[1, 1]
d = kde_data["global_h0.10"]
im = ax.pcolormesh(d["th1"], d["th3"], d["den"], cmap="viridis", shading="auto")
ax.set_title("Test 2b: KDE global, h=0.10")
ax.set_xlabel(r"$\theta_1$")
ax.set_ylabel(r"$\theta_3$")
plt.colorbar(im, ax=ax, label="density")

# Comparison: 1D cross-section at theta3 ~ 0
ax = axes[1, 2]
# Simple binning: find row closest to theta3=0
idx3_bin = np.argmin(np.abs(th3_1[0, :]))
ax.plot(th1_1[:, idx3_bin], den_1[:, idx3_bin], 'k-', lw=2, label="Binning 50x50", alpha=0.7)

# KDE global
for h_str, color, ls in [("global_h0.05", "blue", "--"), ("global_h0.10", "red", "-.")]:
    d = kde_data[h_str]
    idx3_kde = np.argmin(np.abs(d["th3"][0, :]))
    ax.plot(d["th1"][:, idx3_kde], d["den"][:, idx3_kde], color=color, ls=ls, lw=1.5,
            label=f"KDE global {h_str.split('_')[1]}")

# KDE local
for h_str, color, ls in [("local_h0.05", "cyan", ":"), ("local_h0.10", "orange", ":")]:
    d = kde_data[h_str]
    idx3_kde = np.argmin(np.abs(d["th3"][0, :]))
    ax.plot(d["th1"][:, idx3_kde], d["den"][:, idx3_kde], color=color, ls=ls, lw=1.5,
            label=f"KDE local {h_str.split('_')[1]}")

ax.set_title(r"Cross-section at $\theta_3 \approx 0$")
ax.set_xlabel(r"$\theta_1$")
ax.set_ylabel("density")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("slice_test_comparison.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: slice_test_comparison.png")
