import numpy as np
import matplotlib.pyplot as plt

gamma = 0.1; T1 = 5.0; T3 = 5.0
dt = 0.001; N_steps = 500000
I1, I2, I3 = 2.0, 2.0, 2.0
p1, p2, p3 = 1.0, 0.0, -0.5
gap = 100; N_record = N_steps // gap
obs_bump = np.zeros(N_record)
t_arr = np.arange(N_record) * dt * gap

for step in range(N_steps):
    if step % gap == 0:
        idx = step // gap
        obs_bump[idx] = max(0, 3-I1) + max(0, 3-I2) + max(0, 3-I3)
    M = I1 + I2 + I3
    s12 = np.sin(2*(p1-p2)); c12 = np.cos(2*(p1-p2))
    s32 = np.sin(2*(p3-p2)); c32 = np.cos(2*(p3-p2))
    dI1 = 4*I1*I2*s12 + 2*gamma*(2*T1 - (2*M*I1 - I1**2 + 2*I2*I1*c12))
    dI2 = 4*I2*(-I1*s12 - I3*s32)
    dI3 = 4*I3*I2*s32 + 2*gamma*(2*T3 - (2*M*I3 - I3**2 + 2*I2*I3*c32))
    dp1 = 2*M - I1 + 2*I2*c12 + gamma*2*I2*s12
    dp2 = 2*M - I2 + 2*I1*c12 + 2*I3*c32
    dp3 = 2*M - I3 + 2*I2*c32 + gamma*2*I2*s32
    nI1 = np.random.randn() * 2*np.sqrt(2*gamma*T1*max(I1,1e-14))
    nI3 = np.random.randn() * 2*np.sqrt(2*gamma*T3*max(I3,1e-14))
    np1 = np.random.randn() * np.sqrt(2*gamma*T1/max(I1,1e-14))
    np3 = np.random.randn() * np.sqrt(2*gamma*T3/max(I3,1e-14))
    I1 = max(I1 + dI1*dt + nI1*np.sqrt(dt), 1e-14)
    I2 = max(I2 + dI2*dt, 1e-14)
    I3 = max(I3 + dI3*dt + nI3*np.sqrt(dt), 1e-14)
    p1 = p1 + dp1*dt + np1*np.sqrt(dt)
    p2 = p2 + dp2*dt
    p3 = p3 + dp3*dt + np3*np.sqrt(dt)

fig, axes = plt.subplots(2, 1, figsize=(14, 7))
axes[0].plot(t_arr, obs_bump, linewidth=0.5)
axes[0].set_xlabel('time'); axes[0].set_ylabel('f(x)')
axes[0].set_title('max(0,3−I₁)+max(0,3−I₂)+max(0,3−I₃)')

axes[1].plot(t_arr[:500], obs_bump[:500], linewidth=0.8)
axes[1].set_xlabel('time'); axes[1].set_ylabel('f(x)')
axes[1].set_title('first 50 time units')

plt.tight_layout(); plt.savefig('trajectory_bump.png', dpi=150); plt.show()