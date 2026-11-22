"""batch_09: Sec. 3.3 (action-current conductivity) with nls_cartesian_v3 (= v2 + per-stream bond currents).
Reproduces the two figures of Sec. 3.3 with the current integrator:
  Fig. 6 (linear response): n = 15, 25, 35, 45; T1 = 5 - dT/2, Tn = 5 + dT/2, dT = 0.25, 0.5, 1.0, ..., 4.0;   -> group resp
  Fig. 7 (chain length):    (T1, Tn) = (4, 6); n = 20, 30, ..., 250                                          -> group len
All four couplings BC1, BC2, BC3, BC3b; open chain; gamma 0.1; dt 1e-3; 160 trajectories (10 batches of 16);
burn-in max(2000, 3.2 n^2) (Sec. 3.1 rule; the old Sec. 3.3 runs used 2000 for every n); measurement 20000
(as the old runs); sampling interval 0.1; 4 checkpoints.
  tasks_A.csv: resp (144) + len with n <= 150 (56)   tasks_B.csv: len with n = 160..250 (40)
Usage: python3 make_tasks.py   (writes both files and prints the cost estimate)"""
import csv
BCS = ("BC1", "BC2", "BC3", "BC3b"); MEAS = 20000; BATCH = 10; SEED0 = 202611090000
SEC_PER_SITE_TIME = 4.8e-4      # seconds per batch (16 lanes) per site per unit time, from dense n = 150 on Unity
burn = lambda n: max(2000, round(3.2 * n * n))
rows = []
for bc in BCS:
    for n in (15, 25, 35, 45):
        for dT in (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
            T1, Tn = 5 - dT / 2, 5 + dT / 2
            rows.append(("resp", f"resp_{bc}_dT{dT:g}_n{n}", bc, f"{T1:g}", f"{Tn:g}", n))
for bc in BCS:
    for n in range(20, 251, 10):
        rows.append(("len", f"len_{bc}_T4_T6_n{n}", bc, "4", "6", n))
A = [r for r in rows if not (r[0] == "len" and r[5] > 150)]; Bl = [r for r in rows if r[0] == "len" and r[5] > 150]
seed = SEED0
for name, rs in (("tasks_A.csv", A), ("tasks_B.csv", Bl)):
    with open(name, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["task", "group", "label", "bc", "T1", "Tn", "n", "closed", "burn", "measure", "batches", "dt", "seed"])
        core = 0; wall = 0
        for t, (g, lab, bc, T1, Tn, n) in enumerate(rs):
            w.writerow([t, g, lab, bc, T1, Tn, n, 0, burn(n), MEAS, BATCH, "1e-3", seed]); seed += 1
            per_batch = SEC_PER_SITE_TIME * n * (burn(n) + MEAS)
            core += per_batch * BATCH; wall = max(wall, per_batch)
    print(f"{name}: {len(rs)} tasks, about {core/3600:.0f} core-hours, longest task about {wall/3600:.1f} h on {BATCH} cores")
