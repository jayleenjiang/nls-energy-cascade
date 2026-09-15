import numpy as np
import matplotlib.pyplot as plt

gamma = 0.001; T1 = 5.0; T3 = 5.0
dt = 0.001; N_steps = 1000000
gap = 10; N_record = N_steps // gap
t_arr = np.arange(N_record) * dt * gap

fig, ax = plt.subplots(figsize=(12, 5))

for trial, (i1, i2, i3, pp1, pp3) in enumerate([
    (2.0, 2.0, 2.0, 1.0, -0.5),
    (1.0, 3.0, 1.0, 0.5, 0.5),
    (4.0, 1.0, 3.0, -1.0, 2.0),
    (0.5, 0.5, 0.5, 2.0, -1.0),
    (3.0, 3.0, 3.0, 0.0, 1.0)]):
    
    I1, I2, I3 = i1, i2, i3
    p1, p2, p3 = pp1, 0.0, pp3
    obs = np.zeros(N_record)
    
    for step in range(N_steps):
        if step % gap == 0:
            obs[step//gap] = I1
        M = I1+I2+I3
        s12=np.sin(2*(p1-p2));c12=np.cos(2*(p1-p2))
        s32=np.sin(2*(p3-p2));c32=np.cos(2*(p3-p2))
        dI1=4*I1*I2*s12+2*gamma*(2*T1-(2*M*I1-I1**2+2*I2*I1*c12))
        dI2=4*I2*(-I1*s12-I3*s32)
        dI3=4*I3*I2*s32+2*gamma*(2*T3-(2*M*I3-I3**2+2*I2*I3*c32))
        dp1=2*M-I1+2*I2*c12+gamma*2*I2*s12
        dp2=2*M-I2+2*I1*c12+2*I3*c32
        dp3=2*M-I3+2*I2*c32+gamma*2*I2*s32
        nI1=np.random.randn()*2*np.sqrt(2*gamma*T1*max(I1,1e-14))
        nI3=np.random.randn()*2*np.sqrt(2*gamma*T3*max(I3,1e-14))
        np1=np.random.randn()*np.sqrt(2*gamma*T1/max(I1,1e-14))
        np3=np.random.randn()*np.sqrt(2*gamma*T3/max(I3,1e-14))
        I1=max(I1+dI1*dt+nI1*np.sqrt(dt),1e-14)
        I2=max(I2+dI2*dt,1e-14)
        I3=max(I3+dI3*dt+nI3*np.sqrt(dt),1e-14)
        p1+=dp1*dt+np1*np.sqrt(dt);p2+=dp2*dt;p3+=dp3*dt+np3*np.sqrt(dt)
    
    obs_c = obs - obs.mean()
    acf = np.correlate(obs_c, obs_c, mode='full')
    acf = acf[len(acf)//2:]; acf /= acf[0]
    ax.plot(t_arr[:500], acf[:500], linewidth=1, label=f'IC {trial+1}')

ax.set_xlabel('lag (time units)'); ax.set_ylabel('autocorrelation')
ax.set_title('I₁ autocorrelation at γ=0.001 — 5 different initial conditions')
ax.axhline(0, color='k', linewidth=0.5)
ax.legend(); ax.set_xlim(0, 5)
plt.tight_layout(); plt.savefig('autocorr_multi_IC.png', dpi=150); plt.show()