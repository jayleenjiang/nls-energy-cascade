"""Fig. 8 (heat_ft_grid.pdf): the fluctuation-relation test R_tau(q) against Delta_beta * q for three temperature pairs.

Data: fig8_data_extracted_from_pdf.json. The per-q histogram ratios behind this figure were computed from raw heat
blocks that are no longer stored, so the points, error bars and fit lines were read back from the vector graphics of
the original figure file (heat_ft_paper_grid.pdf) and are redrawn here in the paper style.
The a_inf intervals quoted in the caption are in ../results/affinity_sweep_n10/extrapolation_summary.csv.
Usage: python3 make_figures.py [output folder]   (default: this folder)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "common"))
import paperstyle as S          # shared paper style: fonts, colors, panel titles
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
D = json.load(open(os.path.join(HERE, "fig8_data_extracted_from_pdf.json")))
PAIRS = ["(6.5,5.5)", "(7,5)", "(8,4)"]

fig = plt.figure(figsize=(S.TEXTWIDTH, 2.3))
gs = fig.add_gridspec(1, 3)
L = S.letters()
for col in range(3):
    d = D[str(col)]["a"]
    a = fig.add_subplot(gs[0, col])
    lim = d["lim"]
    for i, tau in enumerate((5, 20, 80)):          # one color per window length tau
        x, y = d["pts"][i]
        c = S.MATLAB[i]
        a.plot(x, y, "o", ms=2, mfc="none", mec=c, mew=0.5, label=rf"$\tau={tau}$")
        for ex, lo, hi in d["err"][str(i)]:
            a.plot([ex, ex], [lo, hi], color=c, lw=0.5)
        (x0, y0), (x1, y1) = d["fit"][str(i)]
        a.plot([x0, x1], [y0, y1], color=c, lw=0.8)
    a.plot([0, lim], [0, lim], "k:", lw=0.8, label=r"$\Delta\beta\,q$")   # the fluctuation-relation prediction
    a.set_xlim(0, lim); a.set_ylim(0, lim)
    a.set_xlabel(r"$\Delta\beta\,q$"); a.set_ylabel(r"$R_\tau(q)$")
    S.panel_title(a, L[col], f"$(T_1,T_n)={PAIRS[col]}$")
    if col == 0:
        a.legend(loc="upper left")
fig.tight_layout(w_pad=0.5, h_pad=0.6)
fig.savefig(os.path.join(OUT, "heat_ft_grid.pdf"))
print("written", os.path.join(OUT, "heat_ft_grid.pdf"))
