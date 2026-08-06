"""Graphical abstract: same optics, ranking set by the collection assumption.

(a) Measured active-layer absorptance of the planar (D1) and embedded (D2)
    cells, with the 1.55 eV / 800 nm host edge marked; beyond it the embedded
    blend absorbs more, and that absorption is what the strict single-gap
    reading discards.
(b) Power conversion efficiency of the two morphologies under the two
    collection readings: parity under the intermediate-band reading, a clear
    planar lead under the strict single-gap reading.

Real data (results/fdtd/*.mat) in the shared muted style. Exports PNG + PDF.
"""
from __future__ import annotations

import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CD1, CD2, CANNOT, CGRID, apply_style, panel_label, grid
apply_style()

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, "results", "fdtd")

with h5py.File(os.path.join(RES, "planar_results.mat"), "r") as h:
    lam1 = np.array(h["lam"]).ravel() * 1e9
    A_d1 = np.clip(np.array(h["Apvk"]).ravel() + np.array(h["Aqd"]).ravel(), 0, 1)
with h5py.File(os.path.join(RES, "embedded_5nm_results.mat"), "r") as h:
    lam2 = np.array(h["lam"]).ravel() * 1e9
    A_d2 = np.clip(np.array(h["Aactive"]).ravel(), 0, 1)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.0, 3.5),
                               gridspec_kw={"width_ratios": [1.35, 1.0]})

# ---- (a) absorptance -------------------------------------------------------
axA.axvspan(800, 1300, color=CGRID, alpha=0.16, zorder=0)
axA.plot(lam1, A_d1, color=CD1, lw=1.9, label="D1 planar")
axA.plot(lam2, A_d2, color=CD2, lw=1.9, label="D2 embedded")
axA.axvline(800, color=CANNOT, lw=0.9, ls=":")
axA.text(812, 0.09, "sub-gap NIR\ndiscarded under\nsingle-gap reading",
         fontsize=7.4, color=CANNOT, va="bottom", ha="left", linespacing=1.25)
axA.set_xlim(320, 1300); axA.set_ylim(0, 1.05)
axA.set_xlabel("Wavelength (nm)")
axA.set_ylabel("Active-layer absorptance")
axA.legend(loc="upper right", fontsize=8.6, frameon=False,
           handlelength=1.4, borderaxespad=0.4)
grid(axA); panel_label(axA, "(a)")

# ---- (b) PCE under the two readings ---------------------------------------
regimes = ["Intermediate\nband", "Strict\nsingle gap"]
d1 = np.array([24.00, 18.80]); d2 = np.array([24.21, 14.23])
x = np.arange(2); w = 0.36
b1 = axB.bar(x - w / 2, d1, w, color=CD1, label="D1 planar", zorder=3)
b2 = axB.bar(x + w / 2, d2, w, color=CD2, label="D2 embedded", zorder=3)
for bars, vals in ((b1, d1), (b2, d2)):
    for r, v in zip(bars, vals):
        axB.text(r.get_x() + r.get_width() / 2, v + 0.35, f"{v:.1f}",
                 ha="center", va="bottom", fontsize=8.2)
axB.set_xticks(x); axB.set_xticklabels(regimes)
axB.set_ylabel("Power conversion efficiency (%)")
axB.set_ylim(0, 28.5)
# parity bracket over the intermediate-band cluster
axB.plot([-w / 2, w / 2], [25.6, 25.6], color=CANNOT, lw=0.8)
axB.text(0, 26.0, "parity", ha="center", va="bottom", fontsize=8, color=CANNOT)
# gap arrow over the single-gap cluster
axB.annotate("", xy=(1 + w / 2 + 0.02, 18.80), xytext=(1 + w / 2 + 0.02, 14.23),
             arrowprops=dict(arrowstyle="<->", color=CANNOT, lw=0.9))
axB.text(1 + w / 2 + 0.12, 16.5, "4.6\npts", ha="left", va="center",
         fontsize=8, color=CANNOT, linespacing=1.05)
grid(axB); panel_label(axB, "(b)")

plt.tight_layout(w_pad=1.6)
out = os.path.join(REPO, "figures", "fig_graphical_abstract.png")
plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(out.replace(".png", ".pdf"), bbox_inches="tight", facecolor="white")
print("wrote", out)
