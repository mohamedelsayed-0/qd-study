"""Figure: Han 2018 absorber-level validation.

(a) Active-layer absorptance for pure MAPbI3 versus the Han-loading blend,
    with the 800 nm MAPbI3 collection edge marked. The dots add absorption only
    on the long-wavelength side of that edge.
(b) Optical current and its above-gap fraction versus dot fill, spanning Han's
    loadings (0.5-1 mg/mL ~ 3-7e-4) up to the 20 vol% used elsewhere in this
    work.

Regenerable via `python scripts/make_fig_han_validation.py` once
scripts/sweeps/han_validation.py has produced results/fdtd/han_validation.csv.
"""
from __future__ import annotations

import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CD1, CD2, CACC, CANNOT, CGRID, apply_style, panel_label, grid
apply_style()

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(REPO, "figures")
RES = os.path.join(REPO, "results", "fdtd")

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5))

# ---- panel (a): absorptance, pure vs Han vs 20% -----------------------------
# Pure and the Han-loading blend are visually indistinguishable (the dots are
# negligible at 3e-4); the 20 vol% curve shows what a real dot contribution
# looks like -- and that all the added absorption lands beyond the 800 nm edge.
ax = axes[0]
ax.axvspan(800, 1200, color=CGRID, alpha=0.18, zorder=0)   # sub-gap: not collected
# pure drawn thick underneath, Han thin on top: where they coincide the orange
# tracks inside the grey band, so "indistinguishable" reads as such rather than
# as a missing curve.
curves = [
    ("han_absorptance_f1e-06.csv", CACC, 2.9, "pure MAPbI$_3$"),
    ("han_absorptance_f3e-04.csv", CD2,  1.3, "blend, 0.5 mg mL$^{-1}$ (Han)"),
    ("han_absorptance_f2e-01.csv", CD1,  1.9, "blend, 20 vol.%"),
]
for fname, col, lwv, lab in curves:
    p = os.path.join(RES, fname)
    if os.path.isfile(p):
        d = np.loadtxt(p, delimiter=",", skiprows=1)
        ax.plot(d[:, 0], d[:, 1], color=col, lw=lwv, label=lab)
ax.axvline(800, color=CANNOT, lw=0.9, ls=":")
ax.text(806, 0.06, "sub-gap\n(not collected)", fontsize=7.4, color=CANNOT,
        ha="left", va="bottom")
ax.set_xlim(300, 1200); ax.set_ylim(0, 1.08)
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Active-layer absorptance")
ax.legend(loc="upper right", fontsize=8.0)
grid(ax); panel_label(ax, "(a)")

# ---- panel (b): Jopt and above-gap fraction vs fill -------------------------
ax = axes[1]
if os.path.isfile(os.path.join(RES, "han_validation.csv")):
    v = np.loadtxt(os.path.join(RES, "han_validation.csv"), delimiter=",", skiprows=1)
    v = np.atleast_2d(v)
    fill, _mg, jopt, agf = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
    fplot = np.where(fill <= 0, 1e-5, fill)   # place "pure" at the axis floor
    ax.semilogx(fplot, jopt, "-o", color=CD2, lw=1.8, ms=6, mec="white", mew=0.7,
                label="$J_{\\mathrm{opt}}$")
    ax.set_ylim(jopt.min() - 1.2, jopt.max() + 1.2)   # keep points off the axis
    ax.set_xlabel("PbS dot volume fraction")
    ax.set_ylabel("$J_{\\mathrm{opt}}$ (mA cm$^{-2}$)", color=CD2)
    ax.tick_params(axis="y", labelcolor=CD2)
    ax.axvspan(3.3e-4, 6.6e-4, color="#d9d2c5", alpha=0.8, zorder=0)
    ax.text(4.7e-4, ax.get_ylim()[1], "Han\nloading", fontsize=7.6,
            color=CANNOT, ha="center", va="top")
    axr = ax.twinx()
    axr.semilogx(fplot, 100*agf, "-s", color=CD1, lw=1.6, ms=5, mec="white", mew=0.7,
                 label="above-gap %")
    axr.set_ylabel("above-gap fraction (%)", color=CD1)
    axr.tick_params(axis="y", labelcolor=CD1)
    axr.set_ylim(0, 100)
grid(ax); panel_label(ax, "(b)")

plt.tight_layout()
out = os.path.join(FIG, "fig_han_validation.png")
plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print("wrote", out)
