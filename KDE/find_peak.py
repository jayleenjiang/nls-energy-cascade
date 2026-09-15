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