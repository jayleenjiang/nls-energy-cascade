"""Mean of the last t with M(t) >= 30 per run (first table of the worst-case section of report_15), and the number of
trajectories still >= 30 at the end of the run (counted with the run length, so the mean is then a lower bound).
Usage (from the batch_07 folder): python3 last_passage_mean.py   -> prints one line per run"""
import numpy as np
for l in open("tasks.tsv"):
    d, nm, a = l.rstrip("\n").split("\t")
    if d != "worst":
        continue
    t_ = a.split(); g = lambda k: t_[t_.index(k) + 1]; n = int(g("--n")); D = 2 * n - 1
    x = np.fromfile(f"out/worst/{nm}_dump.f32", dtype=np.float32).reshape(8, -1, 16, D)
    M = x[..., :n].sum(-1).transpose(0, 2, 1).reshape(128, -1)
    t = (np.arange(M.shape[1]) + 1) * 0.1; above = M >= 30
    last = np.array([t[np.where(r)[0][-1]] if r.any() else 0.0 for r in above]); nf = above[:, -1]
    print(nm, {"mean_last30": float(np.where(nf, t[-1], last).mean()), "not_finished": int(nf.sum()),
               "frac_last30_gt5": float((np.where(nf, np.inf, last) > 5).mean())})
