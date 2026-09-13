import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.01  # 改：gap=10, dt=0.001
tt = np.arange(Y.shape[1]) * dt_eff

ss = Y[:, -50:].mean(axis=1)
signal_strength = np.abs(Y[:, 0] - ss)

def model_osc(t, A, lamR, lamI, phi, c):
    return A * np.exp(lamR * t) * np.sin(lamI * t + phi) + c

def model_pure(t, A, lamR, c):
    return A * np.exp(lamR * t) + c

# Use medium-signal points
mid_idx = np.argsort(signal_strength)[len(Y)//2 : len(Y)//2 + 100]

lamR_list = []
lamI_list = []
ratios = []
for i in mid_idx:
    try:
        popt_osc, _ = curve_fit(model_osc, tt, Y[i,:],
                                p0=[0.1, -1.0, 1.0, 0.0, ss[i]],
                                maxfev=10000,
                                bounds=([-5, -5, 0, -np.pi, -2], [5, 0, 5, np.pi, 2]))
        popt_pure, _ = curve_fit(model_pure, tt, Y[i,:],
                                 p0=[0.1, -1.0, ss[i]], maxfev=10000)
        res_osc = np.mean((Y[i,:] - model_osc(tt, *popt_osc))**2)
        res_pure = np.mean((Y[i,:] - model_pure(tt, *popt_pure))**2)
        if abs(popt_osc[1]) < 4.9 and popt_osc[2] < 4.9:
            lamR_list.append(popt_osc[1])
            lamI_list.append(popt_osc[2])
            ratios.append(res_osc / res_pure)
    except:
        pass

lamR_list = np.array(lamR_list)
lamI_list = np.array(lamI_list)
ratios = np.array(ratios)
print(f"Good fits: {len(lamR_list)}/100")
print(f"λR = {lamR_list.mean():.3f} ± {lamR_list.std():.3f}")
print(f"λI = {lamI_list.mean():.3f} ± {lamI_list.std():.3f}")
print(f"Mean residual ratio (osc/pure): {ratios.mean():.3f}")


from scipy.signal import find_peaks

count_with_peaks = 0
count_without = 0
periods = []

for i in range(len(Y)):
    ss_i = Y[i, -50:].mean()
    signal = Y[i,:] - ss_i
    peaks, _ = find_peaks(signal, distance=20)
    troughs, _ = find_peaks(-signal, distance=20)
    
    # Need at least 2 peaks or 2 troughs for a period
    if len(peaks) >= 2:
        per = np.mean(np.diff(tt[peaks]))
        periods.append(per)
        count_with_peaks += 1
    else:
        count_without += 1

print(f"Trajectories with ≥2 peaks: {count_with_peaks}/{len(Y)}")
print(f"Trajectories without: {count_without}/{len(Y)}")
if len(periods) > 0:
    periods = np.array(periods)
    print(f"Period: {periods.mean():.3f} ± {periods.std():.3f}")
    print(f"→ λI = 2π/period = {2*np.pi/periods.mean():.3f}")


from scipy.signal import find_peaks

count_with = 0
periods = []

# Only look at signal region (first 1.0 time units = first 100 points)
t_cut = 100

for i in range(len(Y)):
    ss_i = Y[i, -50:].mean()
    signal = Y[i,:t_cut] - ss_i
    amp = np.abs(signal).max()
    
    # Require peaks to be at least 20% of max amplitude
    peaks, _ = find_peaks(signal, distance=10, prominence=0.2*amp)
    troughs, _ = find_peaks(-signal, distance=10, prominence=0.2*amp)
    
    all_extrema = np.sort(np.concatenate([peaks, troughs]))
    if len(all_extrema) >= 3:  # need 3 extrema for one full oscillation
        per = 2 * np.mean(np.diff(tt[all_extrema]))
        periods.append(per)
        count_with += 1

print(f"Trajectories with real oscillation: {count_with}/{len(Y)}")
if len(periods) > 0:
    periods = np.array(periods)
    print(f"Period: {periods.mean():.3f} ± {periods.std():.3f}")
    print(f"→ λI = 2π/period = {2*np.pi/periods.mean():.3f}")
else:
    print("No oscillation detected — λI = 0 confirmed")

count_crossings = []

for i in range(len(Y)):
    ss_i = Y[i, -50:].mean()
    signal = Y[i,:] - ss_i
    
    # Only count zero crossings in first 2 time units (first 200 points)
    sig = signal[:200]
    crossings = np.sum(np.diff(np.sign(sig)) != 0)
    count_crossings.append(crossings)

count_crossings = np.array(count_crossings)
print(f"Zero crossings in first 2 time units:")
print(f"  Mean: {count_crossings.mean():.1f}")
print(f"  Median: {np.median(count_crossings):.0f}")
print(f"  0 crossings: {(count_crossings==0).sum()}")
print(f"  1 crossing:  {(count_crossings==1).sum()}")
print(f"  2 crossings: {(count_crossings==2).sum()}")
print(f"  3+ crossings: {(count_crossings>=3).sum()}")

mask_0 = count_crossings == 0
mask_3 = count_crossings >= 3

print("0 crossings — initial conditions:")
print(f"  mean I1={X[mask_0,0].mean():.2f}, I2={X[mask_0,1].mean():.2f}, I3={X[mask_0,2].mean():.2f}")

print("3+ crossings — initial conditions:")
print(f"  mean I1={X[mask_3,0].mean():.2f}, I2={X[mask_3,1].mean():.2f}, I3={X[mask_3,2].mean():.2f}")