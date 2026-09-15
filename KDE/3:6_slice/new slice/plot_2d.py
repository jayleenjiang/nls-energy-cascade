#!/usr/bin/env python3
"""Plot 2D density of (theta1, theta3) from NLS5D_slice3 output."""

import numpy as np
import matplotlib.pyplot as plt
import os, glob

def load_grid(fname):
    d = np.loadtxt(fname, delimiter=',', skiprows=1)
    th1 = np.unique(d[:,0]); th3 = np.unique(d[:,1])
    Z = d[:,2].reshape(len(th1), len(th3))
    return th1, th3, Z

# Detect which files exist
files = {
    'bin':    'density_bin.csv',
    'kde05':  'density_kde_h0.05.csv',
    'kde10':  'density_kde_h0.10.csv',
    'theory': 'density_theory.csv',
}
avail = {k: v for k, v in files.items() if os.path.isfile(v)}
has_theory = 'theory' in avail

# Decide layout
panels = []
if 'bin' in avail:     panels.append(('bin',    'Binning 50×50'))
if 'kde05' in avail:   panels.append(('kde05',  'KDE h=0.05'))
if 'kde10' in avail:   panels.append(('kde10',  'KDE h=0.10'))
if has_theory:          panels.append(('theory', 'Theory'))

n = len(panels)
if n == 0:
    print("No density CSV files found. Run NLS5D_slice3 first.")
    exit(1)

ncols = min(n, 2)
nrows = (n + ncols - 1) // ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False)

# Global colorbar range from all data
vmin, vmax = 1e30, -1e30
grids = {}
for key, label in panels:
    th1, th3, Z = load_grid(avail[key])
    grids[key] = (th1, th3, Z)
    vmin = min(vmin, Z.min())
    vmax = max(vmax, Z.max())

for idx, (key, label) in enumerate(panels):
    r, c = divmod(idx, ncols)
    ax = axes[r][c]
    th1, th3, Z = grids[key]
    im = ax.pcolormesh(th3, th1, Z, shading='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    ax.set_xlabel(r'$\theta_3 = 2(\phi_3 - \phi_2)$')
    ax.set_ylabel(r'$\theta_1 = 2(\phi_1 - \phi_2)$')
    ax.set_title(label)
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, shrink=0.8)

# Hide unused axes
for idx in range(n, nrows*ncols):
    r, c = divmod(idx, ncols)
    axes[r][c].set_visible(False)

fig.suptitle(r'Invariant density on $(\theta_1, \theta_3)$ slice', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('density_2d.png', dpi=200, bbox_inches='tight')
print(f"Saved density_2d.png  ({n} panels)")
plt.close()

# Also make a difference plot if theory is available
if has_theory and 'kde05' in avail:
    th1_t, th3_t, Z_t = grids['theory']
    th1_k, th3_k, Z_k = grids['kde05']
    if Z_t.shape == Z_k.shape:
        diff = Z_k - Z_t
        fig2, ax2 = plt.subplots(figsize=(6,5))
        lim = max(abs(diff.min()), abs(diff.max()))
        im2 = ax2.pcolormesh(th3_t, th1_t, diff, shading='auto', cmap='RdBu_r',
                             vmin=-lim, vmax=lim)
        ax2.set_xlabel(r'$\theta_3$'); ax2.set_ylabel(r'$\theta_1$')
        ax2.set_title('KDE (h=0.05) − Theory')
        ax2.set_aspect('equal')
        fig2.colorbar(im2, ax=ax2, shrink=0.8)
        plt.tight_layout()
        plt.savefig('density_diff.png', dpi=200, bbox_inches='tight')
        print("Saved density_diff.png")
        plt.close()