import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.1
tt = np.arange(Y.shape[1]) * dt_eff

# Find strongest signals
ss = Y[:, -10:].mean(axis=1)
signal_strength = np.abs(Y[:, 0] - ss)
top_idx = np.argsort(signal_strength)[-20:]

def model_osc(t, A, lamR, lamI, phi, c):
    return A * np.exp(lamR * t) * np.sin(lamI * t + phi) + c

def model_pure(t, A, lamR, c):
    return A * np.exp(lamR * t) + c

fig, axes = plt.subplots(4, 5, figsize=(20, 12))
for k, i in enumerate(top_idx):
    ax = axes[k//5, k%5]
    ax.plot(tt, Y[i,:], 'b-', linewidth=1.5)
    
    # Fit both models
    try:
        popt_osc, _ = curve_fit(model_osc, tt, Y[i,:],
                                p0=[0.5, -1.0, 1.5, 0.0, ss[i]],
                                maxfev=10000,
                                bounds=([-10, -10, 0, -np.pi, -2], [10, 0, 10, np.pi, 2]))
        ax.plot(tt, model_osc(tt, *popt_osc), 'r--', linewidth=1.5)
        
        popt_pure, _ = curve_fit(model_pure, tt, Y[i,:],
                                 p0=[0.5, -1.0, ss[i]], maxfev=10000)
        ax.plot(tt, model_pure(tt, *popt_pure), 'g:', linewidth=1.5)
        
        res_osc = np.mean((Y[i,:] - model_osc(tt, *popt_osc))**2)
        res_pure = np.mean((Y[i,:] - model_pure(tt, *popt_pure))**2)
        
        ax.set_title(f'λI={popt_osc[2]:.2f} osc/pure={res_osc/res_pure:.2f}', fontsize=8)
    except:
        ax.set_title('fit failed', fontsize=8)

plt.suptitle('Red=oscillatory fit, Green=pure decay fit', fontsize=13)
plt.tight_layout(); plt.show()


lamR_top = []
lamI_top = []
for i in top_idx:
    try:
        popt_osc, _ = curve_fit(model_osc, tt, Y[i,:],
                                p0=[0.5, -1.0, 1.5, 0.0, ss[i]],
                                maxfev=10000,
                                bounds=([-10, -10, 0, -np.pi, -2], [10, 0, 10, np.pi, 2]))
        lamR_top.append(popt_osc[1])
        lamI_top.append(popt_osc[2])
        print(f"Point {i}: λR={popt_osc[1]:.3f}, λI={popt_osc[2]:.3f}")
    except:
        print(f"Point {i}: fit failed")

lamR_top = np.array(lamR_top)
lamI_top = np.array(lamI_top)
print(f"\nMean λR = {lamR_top.mean():.3f} ± {lamR_top.std():.3f}")
print(f"Mean λI = {lamI_top.mean():.3f} ± {lamI_top.std():.3f}")

for i in range(len(lamR_top)):
    if abs(lamR_top[i]) < 9 and lamI_top[i] < 9:
        print(f"λR={lamR_top[i]:.3f}, λI={lamI_top[i]:.3f}")
        
mask = (np.abs(lamR_top) < 9) & (lamI_top < 9)
print(f"\nGood fits: {mask.sum()}/{len(lamR_top)}")
print(f"λR = {lamR_top[mask].mean():.3f} ± {lamR_top[mask].std():.3f}")
print(f"λI = {lamI_top[mask].mean():.3f} ± {lamI_top[mask].std():.3f}")


# Use medium-signal points instead of strongest
mid_idx = np.argsort(signal_strength)[len(Y)//2 : len(Y)//2 + 100]

lamR_mid = []
lamI_mid = []
for i in mid_idx:
    try:
        popt, _ = curve_fit(model_osc, tt, Y[i,:],
                            p0=[0.1, -1.0, 1.0, 0.0, ss[i]],
                            maxfev=10000,
                            bounds=([-5, -5, 0, -np.pi, -2], [5, 0, 5, np.pi, 2]))
        if abs(popt[1]) < 4.9 and popt[2] < 4.9:
            lamR_mid.append(popt[1])
            lamI_mid.append(popt[2])
    except:
        pass

lamR_mid = np.array(lamR_mid)
lamI_mid = np.array(lamI_mid)
print(f"Good fits: {len(lamR_mid)}/100")
print(f"λR = {lamR_mid.mean():.3f} ± {lamR_mid.std():.3f}")
print(f"λI = {lamI_mid.mean():.3f} ± {lamI_mid.std():.3f}")

for i in top_idx:
    try:
        popt_osc, _ = curve_fit(model_osc, tt, Y[i,:],
                                p0=[0.5, -1.0, 1.5, 0.0, ss[i]],
                                maxfev=10000,
                                bounds=([-10, -10, 0, -np.pi, -2], [10, 0, 10, np.pi, 2]))
        popt_pure, _ = curve_fit(model_pure, tt, Y[i,:],
                                 p0=[0.5, -1.0, ss[i]], maxfev=10000)
        res_osc = np.mean((Y[i,:] - model_osc(tt, *popt_osc))**2)
        res_pure = np.mean((Y[i,:] - model_pure(tt, *popt_pure))**2)
        print(f"Point {i}: ratio={res_osc/res_pure:.4f}  λR={popt_osc[1]:.2f}  λI={popt_osc[2]:.2f}")
    except:
        pass