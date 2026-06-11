import numpy as np
from scipy.optimize import curve_fit

Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.01  # gap=10, dt=0.001
tt = np.arange(Y.shape[1]) * dt_eff

def model_pure(t, A, lamR, c):
    return A * np.exp(lamR * t) + c

rates = []
for i in range(len(Y)):
    try:
        popt, _ = curve_fit(model_pure, tt, Y[i,:],
                           p0=[0.5, -1.0, Y[i,-50:].mean()], maxfev=5000)
        if -5 < popt[1] < 0:
            rates.append(popt[1])
    except:
        pass

rates = np.array(rates)
lamR = np.median(rates)
print(f"λ_R = {lamR:.3f} (median of {len(rates)} fits)")

Tstart = 10; Tend = 200
t_fit = tt[Tstart:Tend]
A_mat = np.column_stack([np.exp(lamR * t_fit), np.ones(len(t_fit))])

Q1_vals = []
X_good = []
for i in range(len(Y)):
    y_fit = Y[i, Tstart:Tend]
    sol, res, _, _ = np.linalg.lstsq(A_mat, y_fit, rcond=None)
    Q1 = sol[0]
    fitted = A_mat @ sol
    SSR = np.mean((y_fit - fitted)**2)
    if SSR < 0.01:
        Q1_vals.append(Q1)
        X_good.append(X[i,:])

Q1_vals = np.array(Q1_vals)
X_good = np.array(X_good)
print(f"Good fits: {len(Q1_vals)}/{len(Y)}")
print(f"Q1 range: [{Q1_vals.min():.3f}, {Q1_vals.max():.3f}]")

np.savetxt('backward_NLS_X.txt', X_good)
np.savetxt('backward_NLS_Q1.txt', Q1_vals)