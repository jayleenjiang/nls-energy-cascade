"""Fig. 1: the CKSTT resonant-family embedding of the first five generations.

Reads the four rectangle files ../geometry/rectangles_N5_R100_f<k>_to_f<k+1>.csv (made by ../geometry/main.cpp)
and writes cascade_embedding.pdf / .png into the folder given as argument (default: this folder).
Each file lists, for every nuclear family of one generation step, the two parents and the two children as complex
frequencies (exact fractions); a family is drawn as the rectangle parent, child, parent, child.
This is the routine make_embedding_figure of the original generate_manuscript_figures.py, with the file paths changed.
"""
import csv
import os
import sys
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEOMETRY = os.path.join(HERE, "..", "geometry")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def read_transition(path):
    """rows of one rectangle file: family id, parent/child, real and imaginary part"""
    with open(path, newline="") as f:
        return [{"family_id": int(r["family_id"]), "point_type": r["point_type"],
                 "real": float(Fraction(r["real"])), "imag": float(Fraction(r["imag"]))} for r in csv.DictReader(f)]


plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10, "legend.fontsize": 9,
                     "figure.dpi": 140, "savefig.dpi": 220})
files = [("$\\Lambda_%d\\to\\Lambda_%d$" % (k, k + 1), os.path.join(GEOMETRY, "rectangles_N5_R100_f%d_to_f%d.csv" % (k, k + 1)))
         for k in range(1, 5)]
fig, ax = plt.subplots(figsize=(6.0, 6.0))
colors = plt.cm.viridis(np.linspace(0.08, 0.88, len(files)))
for (label, path), color in zip(files, colors):
    by_family = {}
    for row in read_transition(path):
        by_family.setdefault(row["family_id"], []).append(row)
    for rows in by_family.values():
        pts = [r for r in rows if r["point_type"] == "parent"] + [r for r in rows if r["point_type"] == "child"]
        xs = [p["real"] for p in pts]; ys = [p["imag"] for p in pts]
        if len(xs) == 4:                       # rectangle parent1 -> child1 -> parent2 -> child2 -> parent1
            order = [0, 2, 1, 3, 0]
            ax.plot([xs[i] for i in order], [ys[i] for i in order], color=color, lw=1.0, alpha=0.55)
        ax.scatter(xs, ys, s=16, color=color, edgecolor="black", linewidth=0.25, label=label)
handles, labels = ax.get_legend_handles_labels()       # one legend entry per generation step
seen = set()
uniq = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
ax.legend([h for h, _ in uniq], [l for _, l in uniq], frameon=False, fontsize=9, loc="upper right")
ax.axhline(0, color="0.65", lw=0.8); ax.axvline(0, color="0.65", lw=0.8)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel(r"$\operatorname{Re} K$"); ax.set_ylabel(r"$\operatorname{Im} K$")
ax.set_title("CKSTT resonant-family embedding, first five generations")
ax.grid(alpha=0.22)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cascade_embedding.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "cascade_embedding.png"), dpi=220, bbox_inches="tight")
print("written", os.path.join(OUT, "cascade_embedding.png"))
