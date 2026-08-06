"""fig0: device architecture schematic for the two morphologies.

A conventional device-stack figure: the planar (D1) and embedded (D2) cells
side by side on an identical layer stack, with quantum-dot placement as the
only variable. Layer thicknesses are drawn to scale.

Regenerable via `python scripts/make_fig0_architecture.py`.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (CD1, CACC, CANNOT, LAYERS, EDGE,
                      FACE_TOP, FACE_SIDE, apply_style, shade)

apply_style()

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(REPO, "figures")

# ---- canvas geometry ------------------------------------------------------
# The stack is intrinsically tall and thin while the labelled layout is wide,
# so an equal aspect ratio is not usable. Instead the axes limits are fixed
# here and the isometric depth offset is derived from the resulting
# inch-per-unit scales, so the perspective angle stays correct regardless.
FIGW, FIGH = 7.1, 3.55
XLIM = (-0.30, 8.70)
YLIM = (-2.35, 6.55)

XSPAN = XLIM[1] - XLIM[0]
YSPAN = YLIM[1] - YLIM[0]
KX = FIGW / XSPAN          # inches per x-unit
KY = FIGH / YSPAN          # inches per y-unit

DEPTH_IN = 0.20            # visual length of the isometric offset, inches
ANGLE = np.deg2rad(30.0)
DX = DEPTH_IN * np.cos(ANGLE) / KX
DY = DEPTH_IN * np.sin(ANGLE) / KY

LW = 0.6
NM_PER_UNIT = 330.0        # nm per y-unit -> 1650 nm stack = 5.0 units
LBL = 7.0
MIN_SEP = 0.30             # minimum vertical gap between fanned labels

# bottom -> top: (label, thickness_nm, palette key)
D1_STACK = [
    ("Au 100",        100, "Au"),
    ("MoO$_x$ 10",     10, "MoOx"),
    ("p-PbS-EDT 50",   50, "p-PbS-EDT"),
    ("PbS-QD 300",    300, "PbS-QD"),
    ("MAPbI$_3$ 550", 550, "MAPbI3"),
    ("SnO$_2$ 40",     40, "SnO2"),
    ("FTO 300",       300, "FTO"),
    ("Glass 300",     300, "Glass"),
]

D2_STACK = [
    ("Au 100",        100, "Au"),
    ("MoO$_x$ 10",     10, "MoOx"),
    ("p-PbS-EDT 50",   50, "p-PbS-EDT"),
    ("MAPbI$_3$:PbS-QD\nblend 850", 850, "Blend"),
    ("SnO$_2$ 40",     40, "SnO2"),
    ("FTO 300",       300, "FTO"),
    ("Glass 300",     300, "Glass"),
]


def draw_layer(ax, x0, y0, w, h, colour, dots=False):
    faces = [
        ([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)], colour),
        ([(x0 + w, y0), (x0 + w + DX, y0 + DY),
          (x0 + w + DX, y0 + h + DY), (x0 + w, y0 + h)],
         shade(colour, FACE_SIDE)),
        ([(x0, y0 + h), (x0 + w, y0 + h),
          (x0 + w + DX, y0 + h + DY), (x0 + DX, y0 + h + DY)],
         shade(colour, FACE_TOP)),
    ]
    for pts, fc in faces:
        ax.add_patch(Polygon(pts, closed=True, facecolor=fc, edgecolor=EDGE,
                             linewidth=LW, joinstyle="miter", zorder=2))
    if dots:
        rng = np.random.default_rng(4)
        n = 300
        xs = rng.uniform(x0 + 0.04, x0 + w - 0.04, n)
        ys = rng.uniform(y0 + 0.04, y0 + h - 0.04, n)
        ax.scatter(xs, ys, s=1.7, c=shade(LAYERS["PbS-QD"], 0.85),
                   alpha=0.85, linewidths=0, zorder=3)


def fan(y_desired, min_sep=MIN_SEP):
    y = list(y_desired)
    for i in range(1, len(y)):
        if y[i] - y[i - 1] < min_sep:
            y[i] = y[i - 1] + min_sep
    return y


def draw_stack(ax, stack, x0, width, panel, name, subtitle, dot_layer=None):
    y, anchors = 0.0, []
    for label, thick_nm, key in stack:
        h = thick_nm / NM_PER_UNIT
        draw_layer(ax, x0, y, width, h, LAYERS[key], dots=(key == dot_layer))
        anchors.append((y + h / 2 + DY / 2, label))
        y += h

    label_x = x0 + width + DX + 0.34
    for (y_anchor, label), y_text in zip(anchors, fan([a for a, _ in anchors])):
        ax.plot([x0 + width + DX, label_x - 0.20, label_x - 0.05],
                [y_anchor, y_text, y_text],
                color=CANNOT, lw=0.5, zorder=1)
        ax.text(label_x, y_text, label, fontsize=LBL, va="center", ha="left",
                color="#1A1A1A", linespacing=1.2)

    # Illumination arrow
    top = y + DY
    ax.add_patch(FancyArrowPatch((x0 + width / 2, top + 0.80),
                                 (x0 + width / 2, top + 0.16),
                                 arrowstyle="-|>", mutation_scale=9,
                                 color=CD1, lw=1.2, zorder=4))
    ax.text(x0 + width / 2, top + 0.90, "AM1.5G", fontsize=7.6,
            ha="center", va="bottom", color=CD1)

    # Panel label + name on one line, descriptor well below it.
    cx = x0 + width / 2 + DX / 2
    ax.text(cx, -0.45, f"{panel}  {name}", fontsize=9.0, ha="center",
            va="top", fontweight="bold")
    ax.text(cx, -1.15, subtitle, fontsize=7.0, ha="center", va="top",
            color=CANNOT, linespacing=1.4)


fig, ax = plt.subplots(figsize=(FIGW, FIGH))
ax.set_axis_off()
ax.set_xlim(*XLIM)
ax.set_ylim(*YLIM)

STACK_W = 0.62
draw_stack(ax, D1_STACK, x0=0.62, width=STACK_W,
           panel="(a)", name="D1 — planar",
           subtitle="discrete 300 nm PbS-QD film\nbeneath the perovskite",
           dot_layer=None)

draw_stack(ax, D2_STACK, x0=4.85, width=STACK_W,
           panel="(b)", name="D2 — embedded",
           subtitle="5 nm PbS QDs dispersed in MAPbI$_3$\n(Maxwell–Garnett medium, 20 vol.%)",
           dot_layer="Blend")

# Scale bar
sb = 500 / NM_PER_UNIT
ax.plot([0.30, 0.30], [0.0, sb], color=CACC, lw=1.1)
for yy in (0.0, sb):
    ax.plot([0.22, 0.38], [yy, yy], color=CACC, lw=1.1)
ax.text(0.14, sb / 2, "500 nm", fontsize=7.0, rotation=90,
        va="center", ha="right", color=CACC)

plt.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
out = os.path.join(FIG, "fig0_architecture.png")
plt.savefig(out, dpi=400, bbox_inches="tight", facecolor="white")
print("wrote %s" % out)
