import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.signal import periodogram

# ============================================================
# Load data
# ============================================================
Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.1  # gap=100, dt=0.001
tt = np.arange(Y.shape[1]) * dt_eff

print("="*60)
print("EIGENFUNCTION ANALYSIS — cos θ₁ observable")
print("="*60)
print(f"Data: {Y.shape[0]} initial points × {Y.shape[1]} time steps")
print(f"Time range: [0, {tt[-1]}] (dt_eff = {dt_eff})")

# ============================================================
# 1. Raw signal and decaying signal
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(14, 8))
for i in range(10):
    axes[0].plot(tt, Y[i,:], linewidth=0.8)
axes[0].set_xlabel('time'); axes[0].set_ylabel('E_x[cos θ₁]')
axes[0].set_title('MC-averaged observable (first 10 initial points)')

for i in range(10):
    ss = Y[i, -10:].mean()
    axes[1].plot(tt, Y[i,:] - ss, linewidth=0.8)
axes[1].set_xlabel('time'); axes[1].set_ylabel('E_x[cos θ₁] − steady state')
axes[1].set_title('After subtracting steady state')
axes[1].axhline(0, color='k', linewidth=0.5)
plt.tight_layout(); plt.savefig('eigen_1_raw_signal.png', dpi=150); plt.show()

# ============================================================
# 2. Zoom + log scale
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 5))
t_cut = 30; idx = tt <= t_cut
for i in range(20):
    ss = Y[i, -10:].mean()
    signal = Y[i,:] - ss
    axes[0].plot(tt[idx], signal[idx], 'o-', linewidth=1, markersize=3)
axes[0].set_xlabel('time'); axes[0].set_ylabel('signal')
axes[0].set_title('Early time (zoom)'); axes[0].axhline(0, color='k', linewidth=0.5)

for i in range(20):
    ss = Y[i, -10:].mean()
    signal = np.abs(Y[i,:] - ss)
    signal[signal < 1e-6] = np.nan
    axes[1].plot(tt[idx], signal[idx], 'o-', linewidth=1, markersize=3)
axes[1].set_yscale('log'); axes[1].set_xlabel('time'); axes[1].set_ylabel('|signal|')
axes[1].set_title('Log scale (slope = λ_R)')
plt.tight_layout(); plt.savefig('eigen_2_zoom.png', dpi=150); plt.show()

# ============================================================
# 3. Nonlinear regression: A * exp(λR*t) * sin(λI*t + φ) + c
# ============================================================
print("\n" + "="*60)
print("NONLINEAR REGRESSION")
print("="*60)

def model_osc(t, A, lamR, lamI, phi, c):
    return A * np.exp(lamR * t) * np.sin(lamI * t + phi) + c

lamR_all = []
lamI_all = []
for i in range(len(Y)):
    try:
        popt, _ = curve_fit(model_osc, tt, Y[i,:],
                           p0=[0.5, -1.0, 1.5, 0.0, -0.25],
                           maxfev=10000,
                           bounds=([-5, -10, 0, -np.pi, -2], [5, 0, 10, np.pi, 2]))
        A, lamR, lamI, phi, c = popt
        if abs(lamR) < 5 and lamI < 8:
            lamR_all.append(lamR)
            lamI_all.append(lamI)
    except:
        pass

lamR_all = np.array(lamR_all)
lamI_all = np.array(lamI_all)
print(f"Good fits: {len(lamR_all)}/{len(Y)}")
print(f"λ_R (all): {lamR_all.mean():.3f} ± {lamR_all.std():.3f}")
print(f"λ_I (all): {lamI_all.mean():.3f} ± {lamI_all.std():.3f}")

# ============================================================
# 4. Histograms of eigenvalues
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].hist(lamR_all, bins=50, color='steelblue', edgecolor='white')
axes[0].set_xlabel('λ_R'); axes[0].set_ylabel('count'); axes[0].set_title('Decay rate λ_R')
axes[0].axvline(lamR_all.mean(), color='red', linestyle='--', label=f'mean = {lamR_all.mean():.2f}')
axes[0].legend()

axes[1].hist(lamI_all, bins=50, color='steelblue', edgecolor='white')
axes[1].set_xlabel('λ_I'); axes[1].set_ylabel('count'); axes[1].set_title('Oscillation frequency λ_I')
axes[1].legend()
plt.tight_layout(); plt.savefig('eigen_3_histograms.png', dpi=150); plt.show()

# ============================================================
# 5. Separate two modes
# ============================================================
print("\n" + "="*60)
print("TWO EIGENMODES DETECTED")
print("="*60)

mask_osc = (lamI_all > 0.5) & (lamI_all < 2.0) & (lamR_all > -3)
mask_dec = lamI_all < 0.3

print(f"\nMode 1 (pure decay):")
print(f"  Fits: {mask_dec.sum()}")
print(f"  λ_R = {lamR_all[mask_dec].mean():.3f} ± {lamR_all[mask_dec].std():.3f}")
print(f"  λ_I ≈ 0")
print(f"  Relaxation time τ = {1/abs(lamR_all[mask_dec].mean()):.2f}")

print(f"\nMode 2 (decaying oscillation):")
print(f"  Fits: {mask_osc.sum()}")
print(f"  λ_R = {lamR_all[mask_osc].mean():.3f} ± {lamR_all[mask_osc].std():.3f}")
print(f"  λ_I = {lamI_all[mask_osc].mean():.3f} ± {lamI_all[mask_osc].std():.3f}")
print(f"  Relaxation time τ = {1/abs(lamR_all[mask_osc].mean()):.2f}")
print(f"  Oscillation period T = {2*np.pi/lamI_all[mask_osc].mean():.2f}")

# ============================================================
# 6. Example fits
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 8))
count = 0
for i in range(len(Y)):
    if count >= 6: break
    try:
        popt, _ = curve_fit(model_osc, tt, Y[i,:],
                           p0=[0.5, -1.0, 1.5, 0.0, -0.25],
                           maxfev=10000,
                           bounds=([-5, -10, 0, -np.pi, -2], [5, 0, 10, np.pi, 2]))
        A, lamR, lamI, phi, c = popt
        if 0.5 < lamI < 2.0 and -3 < lamR < 0:
            ax = axes[count//3, count%3]
            ax.plot(tt, Y[i,:], 'b-', linewidth=1.5, label='MC data')
            ax.plot(tt, model_osc(tt, *popt), 'r--', linewidth=1.5, label='fit')
            ax.set_title(f'λR={lamR:.2f}, λI={lamI:.2f}')
            ax.set_xlabel('time'); ax.legend(fontsize=8)
            count += 1
    except:
        pass
plt.suptitle('Example nonlinear fits (decaying oscillation mode)', fontsize=13)
plt.tight_layout(); plt.savefig('eigen_4_example_fits.png', dpi=150); plt.show()

fig, ax = plt.subplots(figsize=(10, 5))
ax.hist(lamI_all, bins=100, color='steelblue', edgecolor='white')
ax.axvline(0.3, color='red', linestyle='--', linewidth=2, label='cutoff = 0.3')
ax.set_xlabel('λ_I', fontsize=14); ax.set_ylabel('count', fontsize=14)
ax.set_title('λ_I distribution — two clusters visible', fontsize=14)
ax.set_xlim(0, 3)
ax.legend(fontsize=12)
plt.tight_layout(); plt.savefig('eigen_lamI_zoom.png', dpi=150); plt.show()


# ============================================================
# 7. Summary
# ============================================================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"""
Observable: f(x) = cos(θ₁), θ₁ = 2(φ₁ − φ₂)
Parameters: γ = 0.1, T₁ = T₃ = 5.0 (equilibrium)
MC: {len(Y)} initial points, N_sample = 10000

Two eigenmodes detected:
  Mode 1: λ = {lamR_all[mask_dec].mean():.2f}          (pure decay, action relaxation)
  Mode 2: λ = {lamR_all[mask_osc].mean():.2f} ± {lamI_all[mask_osc].mean():.2f}i  (decaying oscillation, phase dynamics)

Physical interpretation:
  - Relaxation time τ ≈ {1/abs(lamR_all[mask_osc].mean()):.1f} time units (much faster than 1/γ = 10)
  - Oscillation period T ≈ {2*np.pi/lamI_all[mask_osc].mean():.1f} time units (energy sloshing between modes)
  - Hamiltonian dynamics accelerate convergence to steady state
  - Oscillatory mode confirms Hamiltonian structure (not pure diffusion)
""")