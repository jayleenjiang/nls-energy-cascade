import numpy as np

Y = np.loadtxt('NLS_backward_Y_train.txt')
X = np.loadtxt('NLS_backward_X_train.txt')
dt_eff = 0.01
tt = np.arange(Y.shape[1]) * dt_eff

# Estimate λ_R from pure decay fit
from scipy.optimize import curve_fit

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
lamR = np.median(rates)  # median more robust than mean
print(f"λ_R = {lamR:.3f} (median of {len(rates)} fits)")

# Extract Q1 using fixed λ_R (Kuramoto method: linear regression A\y)
# signal(t) = Q1 * exp(λR * t) + c
# → linear system: [exp(λR*t), 1] @ [Q1, c] = signal
Tstart = 10   # skip initial transient (0.1 time units)
Tend = 200    # before noise floor (2.0 time units)
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
    if SSR < 0.01:  # good fit
        Q1_vals.append(Q1)
        X_good.append(X[i,:])

Q1_vals = np.array(Q1_vals)
X_good = np.array(X_good)
print(f"Good fits: {len(Q1_vals)}/{len(Y)}")
print(f"Q1 range: [{Q1_vals.min():.3f}, {Q1_vals.max():.3f}]")

np.savetxt('backward_NLS_X.txt', X_good)
np.savetxt('backward_NLS_Q1.txt', Q1_vals)