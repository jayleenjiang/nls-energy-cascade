"""batch_08 task lists: worst-case start at small gamma, run until M has decayed.
Same start as batch_07 and the local small-gamma test: I_1 = I_n = 0.01, middle modes share M0 = 100 equally; (T1, Tn) = (2, 8).
Phases (report convention, as batch_07): worst = all angles 4pi/3; anti = endpoint angles 4pi/3, interior angles pi (n >= 4);
random = random phases. Seed 801 for all (common random numbers across phase types). 32 trajectories, dt 1e-4.
  tasks_A.tsv: gamma = 1e-3, run length 3e4, sampled every 1
  tasks_B.tsv: gamma = 1e-4, run length 3e5, sampled every 10
Run length = 30 / gamma: local test n = 4 (16 trajectories) gave last t with M >= 30 of about 3000-6000 at gamma = 1e-3
(about 300 at 1e-2), i.e. time ~ 1/gamma for the slow trajectories.  Usage: python3 make_tasks.py A > tasks_A.tsv"""
import sys
import numpy as np


def phases(n, interior):
    th = [4 * np.pi / 3] + [interior] * (n - 3) + [4 * np.pi / 3]   # theta_1, theta_3..theta_{n-1}, theta_n
    p = np.zeros(n); p[0] = th[0] / 2                                  # phi_2 = 0, phi_1 = theta_1 / 2
    for j in range(2, n):
        p[j] = p[j - 1] + th[j - 1] / 2
    return ",".join(f"{x:.12f}" for x in p)


part = sys.argv[1]
g, T, iv = {"A": ("1e-3", "30000", "1"), "B": ("1e-4", "300000", "10")}[part]
for n in (3, 4, 5, 6):
    init = ",".join(["0.01"] + [f"{100 / (n - 2):.10g}"] * (n - 2) + ["0.01"])
    base = f"--n {n} --T1 2 --Tn 8 --gamma {g} --init {init} --dt 1e-4 --burn 0 --measure {T} --interval {iv} --streams 32 --seed 801"
    print(f"g{g}\tn{n}_worst\t{base} --phases {phases(n, 4 * np.pi / 3)}")
    if n >= 4:
        print(f"g{g}\tn{n}_anti\t{base} --phases {phases(n, np.pi)}")
    print(f"g{g}\tn{n}_random\t{base}")
