"""Figure + table for 4.3.3 from the part files frozen_map_part_*.json.  Usage: python3 fig_frozen_map.py OUTDIR"""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
O = sys.argv[1]
d = {"maps": {}}
for p in sorted(glob.glob(f"{O}/frozen_map_part_*.json")):
    q = json.load(open(p)); d["grid"] = q["grid"]; d["maps"].update(q.get("maps", {}))
    if "averaging_T5" in q:
        d["averaging_T5"] = q["averaging_T5"]
g = np.array(d["grid"]); keys = [k for k in ["canonical_T2_8", "canonical_T5_5", "hlns"] if k in d["maps"]]
titles = {"canonical_T2_8": "canonical bath, (T1,T3)=(2,8)", "canonical_T5_5": "canonical bath, (T1,T3)=(5,5)", "hlns": "HLNS reduced dynamics"}
fig, axes = plt.subplots(2, len(keys), figsize=(5 * len(keys), 8.5), squeeze=False)
lines = ["F(I1,I3)/I3 at the smallest I1 = 0.05 (the I1 << I3 regime of the HLNS bound F <= -0.48 I3)", "",
         "| I3 | " + " | ".join(titles[k] for k in keys) + " |", "|---" * (len(keys) + 1) + "|"]
for c, k in enumerate(keys):
    F = np.array(d["maps"][k]["F"]); F2 = np.array(d["maps"][k]["F_N256"])
    ax = axes[0, c]
    im = ax.imshow(F / (g[:, None] + g[None, :]), origin="lower", cmap="viridis", vmin=-0.9, vmax=0)
    for i in range(len(g)):
        for j in range(len(g)):
            ax.text(j, i, f"{F[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color="w")
    ax.set_xticks(range(len(g)), [f"{x:g}" for x in g], fontsize=7); ax.set_yticks(range(len(g)), [f"{x:g}" for x in g], fontsize=7)
    ax.set_xlabel("$I_3$"); ax.set_ylabel("$I_1$"); ax.set_title(f"{titles[k]}\ncolour F/(I1+I3) (perfect locking: -0.87), numbers F", fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax = axes[1, c]
    for i in [0, 3, 4, 6, 7]:
        ax.plot(g / 1.0, F[i] / g, "o-", ms=3, label=f"I1={g[i]:g}")
    ax.axhline(-0.48, color="r", ls="--", lw=1, label="HLNS bound -0.48")
    ax.axhline(-np.sqrt(3) / 2, color="grey", ls=":", lw=1, label="perfect locking -0.87")
    ax.set_xscale("log"); ax.set_xlabel("$I_3$"); ax.set_ylabel("$F/I_3$"); ax.set_ylim(-1.2, 0.05); ax.legend(fontsize=6)
    ax.set_title("F/I3 vs I3", fontsize=9)
    d["maps"][k]["max_grid_diff_N256_vs_N512"] = float(np.max(np.abs(F - F2)))
    d["maps"][k]["max_F"] = float(F.max())
for j, I3 in enumerate(g):
    lines.append(f"| {I3:g} | " + " | ".join(f"{np.array(d['maps'][k]['F'])[0, j] / I3:.3f}" for k in keys) + " |")
lines += ["", "max F over the grid (all must be < 0 for the sign claim) and max |F(N=512) - F(N=256)|:", ""]
for k in keys:
    lines.append(f"- {titles[k]}: max F = {d['maps'][k]['max_F']:.4f}; grid difference {d['maps'][k]['max_grid_diff_N256_vs_N512']:.3f}")
if "averaging_T5" in d:
    a = d["averaging_T5"]
    lines += ["", f"Averaging check (T1 = T3 = 5): mean of F over (I1,I3) ~ restricted Gibbs = {a['mean_F']:.3f} ± {a['se']:.3f} ({a['draws']} draws); "
              "full boundary dynamics (Sec. A): Q = -2.153 ± 0.002"]
open(f"{O}/frozen_map_table.md", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
fig.tight_layout(); fig.savefig(f"{O}/fig_frozen_map.png", dpi=150)
json.dump(d, open(f"{O}/frozen_map_merged.json", "w"), indent=1)
