"""
Plot slice test v2 results: simulated vs theoretical Gibbs density
Run in the same directory as the output txt files.
"""
import numpy as np
import matplotlib.pyplot as plt

G_bin = 50
G_kde = 30

# Load data
bin50 = np.loadtxt("slice_test1_binning50.txt")
kde05 = np.loadtxt("slice_test2_kde_global_h0.05.txt")
kde10 = np.loadtxt("slice_test2_kde_global_h0.10.txt")

has_theory = True
try:
    theory30 = np.loadtxt("slice_theory_30x30.txt")
    theory50 = np.loadtxt("slice_theory_50x50.txt")
except:
    has_theory = False

# Reshape
d_bin50 = bin50[:, 2].reshape(G_bin, G_bin)
th1_bin = bin50[:G_bin*G_bin:G_bin, 0]
th3_bin = bin50[:G_bin, 1]

d_kde05 = kde05[:, 2].reshape(G_kde, G_kde)
d_kde10 = kde10[:, 2].reshape(G_kde, G_kde)
th1_kde = kde05[:G_kde*G_kde:G_kde, 0]
th3_kde = kde05[:G_kde, 1]

if has_theory:
    d_th30 = theory30[:, 2].reshape(G_kde, G_kde)
    d_th50 = theory50[:, 2].reshape(G_bin, G_bin)

# ---- Determine common color scale from theory ----
if has_theory:
    vmin = d_th30.min()
    vmax = d_th30.max()
else:
    vmin = min(d_bin50.min(), d_kde05.min(), d_kde10.min())
    vmax = max(d_bin50.max(), d_kde05.max(), d_kde10.max())

# ============================================================
# Figure 1: 2D heatmaps comparison
# ============================================================
if has_theory:
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    fig.suptitle(r"Conditional PDF $f(\theta_1, \theta_3 \mid I_1=I_2=I_3=2)$, $T_1=T_3$", fontsize=15)

    # Theory
    ax = axes[0, 0]
    im = ax.pcolormesh(th1_kde, th3_kde, d_th30, cmap="viridis", shading="auto", vmin=vmin, vmax=vmax)
    ax.set_title("Theoretical Gibbs")
    ax.set_xlabel(r"$\theta_1$"); ax.set_ylabel(r"$\theta_3$")
    plt.colorbar(im, ax=ax)

    # Binning 50x50
    ax = axes[0, 1]
    im = ax.pcolormesh(th1_bin, th3_bin, d_bin50, cmap="viridis", shading="auto", vmin=vmin, vmax=vmax)
    ax.set_title("Simple binning 50x50")
    ax.set_xlabel(r"$\theta_1$"); ax.set_ylabel(r"$\theta_3$")
    plt.colorbar(im, ax=ax)

    # KDE h=0.05
    ax = axes[1, 0]
    im = ax.pcolormesh(th1_kde, th3_kde, d_kde05, cmap="viridis", shading="auto", vmin=vmin, vmax=vmax)
    ax.set_title("KDE global h=0.05")
    ax.set_xlabel(r"$\theta_1$"); ax.set_ylabel(r"$\theta_3$")
    plt.colorbar(im, ax=ax)

    # KDE h=0.10
    ax = axes[1, 1]
    im = ax.pcolormesh(th1_kde, th3_kde, d_kde10, cmap="viridis", shading="auto", vmin=vmin, vmax=vmax)
    ax.set_title("KDE global h=0.10")
    ax.set_xlabel(r"$\theta_1$"); ax.set_ylabel(r"$\theta_3$")
    plt.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig("slice2_heatmaps.png", dpi=150, bbox_inches="tight")
    print("Saved: slice2_heatmaps.png")

# ============================================================
# Figure 2: Cross-sections at theta3 ~ 0 and theta3 ~ pi
# ============================================================
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig2.suptitle(r"Cross-sections: simulated vs theoretical", fontsize=14)

for ax, target_th3, label in [(ax1, 0.0, r"$\theta_3 \approx 0$"),
                                (ax2, 3.0, r"$\theta_3 \approx \pi$")]:
    # Find closest index
    j3_kde = np.argmin(np.abs(th3_kde - target_th3))
    j3_bin = np.argmin(np.abs(th3_bin - target_th3))

    ax.plot(th1_bin, d_bin50[:, j3_bin], 'k-', lw=1.5, alpha=0.5, label="Binning 50x50")
    ax.plot(th1_kde, d_kde05[:, j3_kde], 'b--', lw=1.5, label="KDE h=0.05")
    ax.plot(th1_kde, d_kde10[:, j3_kde], 'r-.', lw=1.5, label="KDE h=0.10")
    if has_theory:
        ax.plot(th1_kde, d_th30[:, j3_kde], 'g-', lw=2.5, label="Theory (Gibbs)")

    ax.set_title(label)
    ax.set_xlabel(r"$\theta_1$")
    ax.set_ylabel("density")
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig("slice2_crosssections.png", dpi=150, bbox_inches="tight")
print("Saved: slice2_crosssections.png")

# ============================================================
# Figure 3: Error map (simulated - theory) if theory available
# ============================================================
if has_theory:
    fig3, axes3 = plt.subplots(1, 3, figsize=(16, 4.5))
    fig3.suptitle("Error: simulated - theoretical", fontsize=14)

    for ax, data, G, th1, th3, title in [
        (axes3[0], d_bin50, G_bin, th1_bin, th3_bin, "Binning 50x50"),
        (axes3[1], d_kde05, G_kde, th1_kde, th3_kde, "KDE h=0.05"),
        (axes3[2], d_kde10, G_kde, th1_kde, th3_kde, "KDE h=0.10"),
    ]:
        if G == G_bin:
            err = data - d_th50
        else:
            err = data - d_th30
        emax = max(abs(err.min()), abs(err.max()))
        im = ax.pcolormesh(th1, th3, err, cmap="RdBu_r", shading="auto",
                           vmin=-emax, vmax=emax)
        ax.set_title(title)
        ax.set_xlabel(r"$\theta_1$"); ax.set_ylabel(r"$\theta_3$")
        plt.colorbar(im, ax=ax, label="error")

    plt.tight_layout()
    plt.savefig("slice2_errors.png", dpi=150, bbox_inches="tight")
    print("Saved: slice2_errors.png")

plt.show()