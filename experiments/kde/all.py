
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for gamma_val in [0.1, 0.01, 0.001]:
    gamma = gamma_val; T1 = 5.0; T3 = 5.0
    dt = 0.001; N_steps = 1000000
    I1, I2, I3 = 2.0, 2.0, 2.0
    p1, p2, p3 = 1.0, 0.0, -0.5
    gap = 10; N_record = N_steps // gap
    obs_I1 = np.zeros(N_record)
    obs_I2 = np.zeros(N_record)
    obs_cos = np.zeros(N_record)
    
    for step in range(N_steps):
        if step % gap == 0:
            idx = step // gap
            obs_I1[idx] = I1; obs_I2[idx] = I2
            obs_cos[idx] = np.cos(2*(p1-p2))
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
    
    t_arr = np.arange(N_record)*dt*gap
    for col, (obs, name) in enumerate([(obs_I1,'I₁'),(obs_I2,'I₂'),(obs_cos,'cos θ₁')]):
        obs_c = obs - obs.mean()
        acf = np.correlate(obs_c, obs_c, mode='full')
        acf = acf[len(acf)//2:]; acf /= acf[0]
        axes[col].plot(t_arr[:5000], acf[:5000], linewidth=1, label=f'γ={gamma_val}')

for col, name in enumerate(['I₁','I₂','cos θ₁']):
    axes[col].set_title(f'Autocorrelation: {name}')
    axes[col].set_xlabel('lag (time units)')
    axes[col].axhline(0, color='k', linewidth=0.5)
    axes[col].legend()
plt.tight_layout(); plt.savefig('autocorr_gamma_comparison.png', dpi=150); plt.show()