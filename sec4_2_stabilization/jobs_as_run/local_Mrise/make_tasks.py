"""batch_07 task list: worst-case start for n = 3..6.
Initial actions: I_1 = I_n = 0.01, the middle modes share M0 equally (I_2 = ... = I_{n-1} = M0/(n-2)).
Initial angles (report convention): theta_1 = 2(phi_1 - phi_2) and theta_n = 2(phi_n - phi_{n-1}) (endpoint bonds, measured
outward); interior bonds theta_j = 2(phi_j - phi_{j-1}), j = 3..n-1 (left to right).
  worst  : all angles 4pi/3 (both endpoints at the unstable equilibrium of I(1 + 2 cos theta))
  anti   : endpoint angles 4pi/3, interior angles pi (neighbouring middle modes in anti-phase); n >= 4 only
  random : random phases (no --phases)
(T1, Tn) = (2, 8), gamma = 0.1, 128 trajectories, seed 701 for all (common random numbers across the three phase types).
M0 = 50 / 100 / 200: dt = 1e-4 / 1e-4 / 5e-5, run length 300 / 800 / 1500, sampled every 0.1.
Plus high-resolution first-time-unit runs (folder hires, M0 = 100, every 1e-3) and one n = 6 NESS run.  Usage: python3 make_tasks.py > tasks.tsv"""
import numpy as np


def phases(n, interior):
    th = [4 * np.pi / 3] + [interior] * (n - 3) + [4 * np.pi / 3]   # theta_1, theta_3..theta_{n-1}, theta_n
    p = np.zeros(n); p[0] = th[0] / 2                                  # phi_2 = 0, phi_1 = theta_1 / 2
    for j in range(2, n):
        p[j] = p[j - 1] + th[j - 1] / 2
    return ",".join(f"{x:.12f}" for x in p)


RUN = {50: ("1e-4", 300), 100: ("1e-4", 800), 200: ("5e-5", 1500)}
for n in (3, 4, 5, 6):
    for M0 in (50, 100, 200):
        dt, T = RUN[M0]
        init = ",".join(["0.01"] + [f"{M0 / (n - 2):.10g}"] * (n - 2) + ["0.01"])
        base = f"--n {n} --T1 2 --Tn 8 --gamma 0.1 --init {init} --dt {dt} --burn 0 --measure {T} --interval 0.1 --streams 128 --seed 701"
        print(f"worst\tn{n}_M{M0}_worst\t{base} --phases {phases(n, 4 * np.pi / 3)}")
        if n >= 4:
            print(f"worst\tn{n}_M{M0}_anti\t{base} --phases {phases(n, np.pi)}")
        print(f"worst\tn{n}_M{M0}_random\t{base}")
for n in (3, 4, 5, 6):                      # first time unit at high resolution (angle escape from 4pi/3), M0 = 100
    init = ",".join(["0.01"] + [f"{100 / (n - 2):.10g}"] * (n - 2) + ["0.01"])
    base = f"--n {n} --T1 2 --Tn 8 --gamma 0.1 --init {init} --dt 1e-4 --burn 0 --measure 1 --interval 1e-3 --streams 128 --seed 701"
    print(f"hires\tn{n}_M100_worst\t{base} --phases {phases(n, 4 * np.pi / 3)}")
    print(f"hires\tn{n}_M100_random\t{base}")
print("ness\tn6_T2_8\t--n 6 --T1 2 --Tn 8 --gamma 0.1 --dt 1e-3 --burn 50 --measure 2000 --interval 0.1 --streams 128 --seed 702")
