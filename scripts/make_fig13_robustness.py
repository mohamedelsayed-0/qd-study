"""fig13: D2 robustness — mobility sensitivity (a) and doping map (b).

Uses the committed CSV files. Regenerable via
`python scripts/make_fig13_robustness.py`.
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
FIG = os.path.join(REPO, "figures")
RES = os.path.join(REPO, "results")





# --- panel a: mobility sensitivity ---
mu = np.loadtxt(os.path.join(RES, "charge", "d2_mobility_sensitivity.csv"),
                delimiter=",", skiprows=1)
mu_vals, jsc, voc, ff, pce, pce_rs1 = mu.T

fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))

ax = axes[0]
ax.plot(mu_vals, pce, "-o", color=CD2, lw=2.0, ms=6, mec="white", mew=0.7,
        label="Ideal")
ax.plot(mu_vals, pce_rs1, "-s", color=CACC, lw=1.7, ms=5, mec="white", mew=0.7,
        alpha=0.75, label="Rs = 1 Ω cm²")
ax.axhline(20.27, color=CGRID, ls=":", lw=0.9)
ax.text(0.022, 20.36, "D2 planar baseline", fontsize=8.0,
        color=CANNOT, ha="left", va="bottom")
ax.set_xscale("log")
ax.set_xlabel("D2 blend mobility (cm² V⁻¹ s⁻¹)")
ax.set_ylabel("PCE (%)")
ax.set_ylim(20, 25)
ax.legend(loc="lower right")
grid(ax)
panel_label(ax, "(a)")

# --- panel b: doping map (dense) ---
# Dense CSV columns: dop_cm3,Jsc,Voc,FF,PCE_Rs0,PCE_Rs0p5,PCE_Rs1,PCE_Rs2,notes
# Missing entries (non-convergent points) leave fields blank — skip.
rows = []
with open(os.path.join(RES, "charge", "d2_doping_map_dense.csv")) as f:
    next(f)
    for line in f:
        p = line.strip().split(",")
        if len(p) < 8 or not p[1] or not p[4]:
            continue
        try:
            rs2 = float(p[7]) if p[7] else None
            rows.append((float(p[0]), float(p[1]), float(p[2]), float(p[3]),
                         float(p[4]), rs2))
        except ValueError:
            continue
dops = np.array([r[0] for r in rows])
jsc_d = np.array([r[1] for r in rows])
voc_d = np.array([r[2] for r in rows])
pce_d = np.array([r[4] for r in rows])
pce_rs2 = np.array([r[5] if r[5] is not None else np.nan for r in rows])

ax = axes[1]
color_d = CD2
ax.plot(dops, pce_d, "-o", color=color_d, lw=2.0, ms=6, mec="white", mew=0.7,
        label="Ideal")
ax.plot(dops, pce_rs2, "-s", color=CACC, lw=1.7, ms=5, mec="white", mew=0.7,
        alpha=0.75, label="Rs = 2 Ω cm²")
ax.set_xscale("log")
ax.set_xlabel("D2 blend p-doping (cm⁻³)")
ax.set_ylabel("PCE (%)")
ax.set_ylim(19, 29)
# annotate high-voltage regime (5e15 - 1e17: Voc keeps rising with doping)
ax.axvspan(5e15, 1.5e17, color="#f2e6d9", alpha=0.6, zorder=0)
ax.text(1e16, 28.3, "high-voltage regime\n(current → voltage tradeoff)",
        fontsize=8.5, color="#666", ha="center", va="top")
ax.legend(loc="lower left")
grid(ax)
panel_label(ax, "(b)")

plt.tight_layout()
out = os.path.join(FIG, "fig13_d2_robustness.png")
plt.savefig(out)
print(f"wrote {out}")
