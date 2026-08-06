"""Rs crossover figure (fig12): PCE vs lumped series resistance for the three
optimized devices. Purely post-processes committed JV curves.

Central paper claim: the D1/D2 ranking at zero Rs is reversed above 0.76
Ohm cm2. The heavily-doped D2 variant beats both across the full Rs range.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os as _os
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from figstyle import CD1, CD2, CACC, apply_style
apply_style()

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, "results", "charge")
FIG = os.path.join(REPO, "figures")


plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5, "axes.labelsize": 11.5,
    "axes.linewidth": 0.9, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def pce_rs(jv_csv, rs):
    d = np.loadtxt(jv_csv, delimiter=",", skiprows=1)
    V, J = d[:, 0], d[:, 1]
    return np.max((V - J * 1e-3 * rs) * J) / 100.0 * 100.0


rs = np.linspace(0, 3, 601)
curves = {
    "D1 optimized (42.7 mA cm$^{-2}$)":     ("final_D1_max_jv.csv", CD1, "-"),
    "D2 optimized (50.2 mA cm$^{-2}$)":     ("final_D2_max_jv.csv", CD2, "-"),
    "D2 high-voltage (36.9 mA cm$^{-2}$)":  ("final_D2_hiV_jv.csv", CD2, "--"),
}
fig, ax = plt.subplots(figsize=(5.2, 3.7))
vals = {}
for label, (fn, color, ls) in curves.items():
    p = os.path.join(RES, fn)
    if not os.path.isfile(p):
        print(f"missing {fn}, skipping")
        continue
    y = np.array([pce_rs(p, r) for r in rs])
    vals[label] = y
    ax.plot(rs, y, color=color, ls=ls, lw=2.2, label=label)
ax.set_xlabel("Series resistance $R_\\mathrm{s}$ ($\\Omega$ cm$^{2}$)")
ax.set_ylabel("PCE (%)")
ax.set_xlim(0, 3)
ax.set_ylim(17, 26.9)
ax.legend(loc="lower left", fontsize=8.6, borderaxespad=0.4,
          labelspacing=0.35)
ax.grid(alpha=0.16, lw=0.5, color="#B8B8B8")

# annotate crossover
labels = list(vals)
if len(labels) >= 2:
    d1v = vals[labels[0]]; d2v = vals[labels[1]]
    ix = np.where(np.diff(np.sign(d1v - d2v)))[0]
    if len(ix):
        cx = rs[ix[0]] + (rs[ix[0]+1] - rs[ix[0]]) * abs(d1v[ix[0]] - d2v[ix[0]]) / (abs(d1v[ix[0]] - d2v[ix[0]]) + abs(d1v[ix[0]+1] - d2v[ix[0]+1]) + 1e-12)
        ax.axvline(cx, color="#888", lw=0.9, ls=":")
        # Callout placed above every curve so it cannot collide with the
        # legend in the lower-left corner.
        ax.annotate(f"crossover  {cx:.2f} $\\Omega$ cm$^{{2}}$",
                    xy=(cx, 22.78), xytext=(cx + 0.16, 26.35),
                    fontsize=8.4, color="#6E6E6E", ha="left", va="top",
                    arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.7))
        print(f"D1/D2 crossover at Rs = {cx:.4f}")

plt.tight_layout()
out = os.path.join(FIG, "fig12_rs_crossover.png")
plt.savefig(out)
print(f"wrote {out}")
