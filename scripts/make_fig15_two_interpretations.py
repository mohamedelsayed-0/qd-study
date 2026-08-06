"""fig15: side-by-side unfiltered vs strict single-gap headlines.

Left panel: PCE vs Rs for both interpretations of D1 optimized and D2
optimized. Right panel: bar chart at Rs=0 and Rs=2.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (CD1, CD2, CACC, CPVK, CQD, CGRID, CANNOT,
                      CSHADE, apply_style, panel_label, grid)
apply_style()

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, "results", "charge")
FIG = os.path.join(REPO, "figures")




def pce_rs(jv_csv, rs):
    d = np.loadtxt(jv_csv, delimiter=",", skiprows=1)
    V, J = d[:, 0], d[:, 1]
    return float(np.max((V - J * 1e-3 * rs) * J) / 100.0 * 100.0)


rs = np.linspace(0, 3, 61)
curves = {
    "D1 optimized (unfiltered)": ("final_D1_max_jv.csv",       CD1, "-"),
    "D2 optimized (unfiltered)": ("final_D2_max_jv.csv",       CD2, "-"),
    "D1 optimized (single-gap)": ("d1_max_filtered_jv.csv",    CD1, ":"),
    "D2 optimized (single-gap)": ("d2_max_filtered_jv.csv",    CD2, ":"),
}

fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.8), gridspec_kw={"width_ratios": [1.2, 1]})

ax = axes[0]
for label, (fn, color, ls) in curves.items():
    p = os.path.join(RES, fn)
    y = [pce_rs(p, r) for r in rs]
    ax.plot(rs, y, color=color, ls=ls, lw=2.0, label=label)
ax.set_xlabel("Series resistance $R_\\mathrm{s}$ ($\\Omega$ cm$^{2}$)")
ax.set_ylabel("PCE (%)")
ax.set_xlim(0, 3); ax.set_ylim(10, 26.5)
ax.legend(loc="lower left", fontsize=8.2, ncol=2, columnspacing=1.0,
          handlelength=1.9, borderaxespad=0.4)
grid(ax)
# Regime annotation sits in the empty band above the curves, not over them.
ax.text(0.985, 0.965, "solid: unfiltered\ndotted: strict single-gap",
        transform=ax.transAxes, fontsize=8.2, va="top", ha="right",
        color=CANNOT)
panel_label(ax, "(a)")

# Right panel: bar chart
ax = axes[1]
labels = ["D1 opt", "D2 opt", "D2 hi-V"]
ib_vals = [24.00, 24.21, 24.91]
sg_vals = [18.80, 14.23, 14.78]
x = np.arange(3); w = 0.35
b1 = ax.bar(x - w/2, ib_vals, w, color=[CD1, CD2, CD2], alpha=0.9,
            edgecolor="none", label="unfiltered")
b2 = ax.bar(x + w/2, sg_vals, w, color=[CD1, CD2, CD2], alpha=0.5,
            edgecolor="none", label="strict single-gap")
for xi, v in zip(x - w/2, ib_vals):
    ax.text(xi, v + 0.3, f"{v:.2f}", ha="center", fontsize=9)
for xi, v in zip(x + w/2, sg_vals):
    ax.text(xi, v + 0.3, f"{v:.2f}", ha="center", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("PCE at $R_\\mathrm{s}=0$ (%)")
# Headroom so the legend clears the tallest bar label (24.91).
ax.set_ylim(0, 34)
ax.set_yticks([0, 5, 10, 15, 20, 25])
ax.legend(loc="upper center", fontsize=8.8, ncol=2, columnspacing=1.0,
          handlelength=1.6, borderaxespad=0.2)
grid(ax, axis="y", alpha=0.14)
panel_label(ax, "(b)")

plt.tight_layout()
out = os.path.join(FIG, "fig15_two_interpretations.png")
plt.savefig(out)
print(f"wrote {out}")
