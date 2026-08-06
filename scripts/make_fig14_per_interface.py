"""fig14: D1 per-interface SRV sensitivity.

Reveals the MAPbI3/PbS-QD internal heterojunction (pvqd) as the dominant
D1 interface — at S=10^4 cm/s the internal junction collapses D1 to 11.1%,
while the outer QD/EDT interface barely moves and the SnO2/MAPbI3
sensitivity is carrier-asymmetric (hole velocity dominates).
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
CQEDT, CPVQD, CSNPV = CQD, CD2, CPVK

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(REPO, "figures")
RES = os.path.join(REPO, "results")





# Read the interface CSV manually (mixed types)
rows = []
with open(os.path.join(RES, "charge", "d1_interface_map.csv")) as f:
    next(f)
    for line in f:
        p = line.strip().split(",")
        if len(p) < 7:
            continue
        try:
            rows.append((p[0], float(p[1]), p[2], float(p[3]), float(p[4]),
                         float(p[5]), float(p[6])))
        except ValueError:
            continue

# Extract per-interface curves
def extract(iface, mode):
    S, P = [], []
    for r in rows:
        if r[0] == iface and r[2] == mode:
            S.append(r[1]); P.append(r[6])
    return np.array(S), np.array(P)

fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.7))

# --- panel a: e+h symmetric ---
ax = axes[0]
baseline_pce = 20.59
ax.axhline(baseline_pce, color="#999", ls=":", lw=0.9)
ax.text(1e6, baseline_pce + 0.15, "baseline", fontsize=8.5, color="#666", ha="right")
for iface, color, label in (("qedt", CQEDT, "QD / EDT (outer)"),
                             ("pvqd", CPVQD, "MAPbI3 / QD (internal)"),
                             ("snpv", CSNPV, "SnO2 / MAPbI3 (outer)")):
    S, P = extract(iface, "e+h")
    # add baseline at S=0
    S = np.concatenate([[1e1], S])
    P = np.concatenate([[baseline_pce], P])
    ax.plot(S, P, "-o", color=color, lw=2.0, ms=6, mec="white", mew=0.7, label=label)
ax.set_xscale("log")
ax.set_xlabel("Interface recombination velocity S (cm s⁻¹, common e = h)")
ax.set_ylabel("D1 PCE (%)")
ax.set_ylim(7, 22)
ax.legend(loc="lower left", fontsize=8.5)
grid(ax)
panel_label(ax, "(a)")

# --- panel b: carrier-asymmetric at S=1e4 and 1e6 ---
ax = axes[1]
# For each interface show e+h vs e_only bar comparison at S=1e4
S_TARGET = 1e4
data = {}
for iface, color, label in (("qedt", CQEDT, "QD/EDT"),
                             ("pvqd", CPVQD, "MAPbI3/QD"),
                             ("snpv", CSNPV, "SnO2/MAPbI3")):
    for mode in ("e+h", "e_only"):
        for r in rows:
            if r[0] == iface and r[2] == mode and abs(r[1] - S_TARGET) < 1e-3:
                data[(iface, mode)] = r[6]
                break
labels = ["QD/EDT", "MAPbI3/QD", "SnO2/MAPbI3"]
colors = [CQEDT, CPVQD, CSNPV]
xpos = np.arange(3)
w = 0.35
eh = [data[("qedt", "e+h")], data[("pvqd", "e+h")], data[("snpv", "e+h")]]
eo = [data[("qedt", "e_only")], data[("pvqd", "e_only")], data[("snpv", "e_only")]]
b1 = ax.bar(xpos - w/2, eh, w, color=colors, alpha=0.9, edgecolor="none", label="e = h symmetric")
b2 = ax.bar(xpos + w/2, eo, w, color=colors, alpha=0.5, edgecolor="none", label="e only (h ~ 0)")
for xi, v in zip(xpos - w/2, eh):
    ax.text(xi, v + 0.15, f"{v:.1f}", ha="center", fontsize=8.5)
for xi, v in zip(xpos + w/2, eo):
    ax.text(xi, v + 0.15, f"{v:.1f}", ha="center", fontsize=8.5)
ax.axhline(baseline_pce, color="#999", ls=":", lw=0.9)
ax.text(2.5, baseline_pce + 0.15, "baseline 20.59%", fontsize=8.5, color="#666")
ax.set_xticks(xpos)
ax.set_xticklabels(labels)
ax.set_ylabel("D1 PCE (%)")
# Headroom above the tallest bar so the legend clears the data.
ax.set_ylim(7, 25.5)
ax.set_yticks([8, 10, 12, 14, 16, 18, 20, 22])
ax.set_title(f"S = 10⁴ cm s⁻¹", fontsize=10, pad=8, color="#555")
ax.legend(loc="upper left", fontsize=8.5, ncol=2, columnspacing=1.0,
          handlelength=1.5, borderaxespad=0.3)
grid(ax, axis="y")
panel_label(ax, "(b)")

plt.tight_layout()
out = os.path.join(FIG, "fig14_per_interface.png")
plt.savefig(out)
print(f"wrote {out}")
