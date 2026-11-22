"""Common figure style (v3: LaTeX Latin Modern text to match the paper, thin gray axes).
Old docstring: for the NLS paper, modelled on the MATLAB-default figures
in Li & Sun (arXiv:2511.12583) and Jiang, Chu & Li (arXiv:2406.12002):
Helvetica-type font, boxed axes with inward ticks, no grid, MATLAB colour order,
thin lines without markers, boxed legends, jet colormap for every 2D map,
panel titles '(a) ...' above each panel, no figure titles."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from cycler import cycler

for f in ["/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-regular.otf",
          "/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-italic.otf",
          "/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-bold.otf",
          "/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-bolditalic.otf"]:
    try: font_manager.fontManager.addfont(f)
    except Exception: pass

MATLAB = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]
CMAP = "jet"
TEXTWIDTH = 6.27         # inches, \textwidth of the paper (A4, 1in margins)

AXC = "0.15"   # MATLAB default axes colour [0.15 0.15 0.15]
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath,amssymb}",
    "font.family": "serif",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.4,
    "axes.edgecolor": AXC, "axes.labelcolor": AXC,
    "xtick.color": AXC, "ytick.color": AXC,
    "axes.prop_cycle": cycler(color=MATLAB),
    "axes.grid": False,
    "axes.titlepad": 4,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.minor.size": 1.2, "ytick.minor.size": 1.2,
    "xtick.major.width": 0.4, "ytick.major.width": 0.4,
    "xtick.minor.width": 0.3, "ytick.minor.width": 0.3,
    "lines.linewidth": 0.8,
    "lines.markersize": 3,
    "patch.linewidth": 0.4,
    "legend.frameon": True, "legend.fancybox": False, "legend.framealpha": 1.0,
    "legend.edgecolor": AXC, "legend.borderpad": 0.3, "legend.handlelength": 1.8,
    "legend.labelspacing": 0.25,
    "figure.dpi": 100, "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "image.cmap": CMAP,
    "pdf.fonttype": 42,
})

def panel_title(ax, letter, text=""):
    ax.set_title(f"({letter}) {text}".rstrip())

def colorbar(im, ax, **kw):
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03, **kw)
    cb.outline.set_linewidth(0.4); cb.outline.set_edgecolor(AXC)
    cb.ax.tick_params(direction="in", width=0.4, length=2, labelsize=8, color=AXC)
    from matplotlib.ticker import MaxNLocator
    cb.locator = MaxNLocator(nbins=3, steps=[1, 2, 5, 10]); cb.update_ticks()
    return cb

def letters():
    return "abcdefghijklmnopqrstuvwxyz"

import numpy as _np
def angle_axis(ax, which="xy"):
    t = [-_np.pi, 0, _np.pi]; lab = [r"$-\pi$", r"$0$", r"$\pi$"]
    if "x" in which: ax.set_xticks(t); ax.set_xticklabels(lab); ax.set_xlim(-_np.pi, _np.pi)
    if "y" in which: ax.set_yticks(t); ax.set_yticklabels(lab); ax.set_ylim(-_np.pi, _np.pi)

def read_states(path):
    """density-code dump: magic(8) n(u64) d(u32) b(u32) float32 data -> (n, d) float64"""
    import struct
    with open(path, "rb") as f:
        magic = f.read(8); n, d, b = struct.unpack("<QII", f.read(16)); x = _np.fromfile(f, "<f4")
    assert magic.startswith(b"NLS") and b == 4 and x.size == n * d, path
    return x.reshape(n, d).astype(_np.float64)

def top_legend(fig, ax, ncol, y=1.0):
    """one framed legend above the panels, so that it never covers a curve"""
    h, l = ax.get_legend_handles_labels()
    return fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, y), ncol=ncol, columnspacing=1.2)
