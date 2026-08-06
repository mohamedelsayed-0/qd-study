"""Shared figure style for the manuscript.

Single source of truth for the palette and rcParams so every figure in the
paper reads as one system. The earlier per-script palette used saturated
screen colours (#1f6feb / #d95f02) that are too vibrant for print; the values
below are desaturated equivalents that keep hue separation, survive greyscale
conversion, and sit comfortably next to body text.

Usage:
    from figstyle import CD1, CD2, CACC, apply_style, panel_label, grid
    apply_style()
"""
from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt

# --- primary series palette ------------------------------------------------
# D1 and D2 are the two morphologies compared throughout the paper. They are
# separated in both hue (cool/warm) and lightness, so the distinction holds in
# greyscale print.
CD1 = "#2E5C8A"   # D1 planar        - slate blue   (L* ~ 38)
CD2 = "#B4652A"   # D2 embedded      - terracotta   (L* ~ 51)
CACC = "#3A3A3A"  # accent / reference curves
CPVK = "#6A5B8C"  # MAPbI3 perovskite - muted violet
CQD = "#3E7C7C"   # PbS quantum dots  - muted teal

# --- supporting greys ------------------------------------------------------
CGRID = "#B8B8B8"
CANNOT = "#6E6E6E"
CSHADE = "#ECE7DF"   # neutral band fill (e.g. "passivated" range)
CSHADE2 = "#E4E9EE"  # secondary band fill

# --- device-stack layer palette (architecture schematic) -------------------
# Muted, low-chroma fills so the schematic does not overpower the data
# figures it sits next to.
LAYERS = {
    "Glass":     "#DFE6EA",
    "FTO":       "#8FA9B8",
    "SnO2":      "#B6C3CA",
    "MAPbI3":    "#C6A96B",
    "PbS-QD":    "#8C6070",
    "Blend":     "#B79A72",
    "p-PbS-EDT": "#B4652A",
    "MoOx":      "#9A8C79",
    "Au":        "#C7A85C",
}
EDGE = "#3A3A3A"      # layer outline
FACE_TOP = 1.10       # lightness multiplier for the isometric top face
FACE_SIDE = 0.82      # lightness multiplier for the isometric side face

RC = {
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.labelsize": 11.5,
    "axes.linewidth": 0.9,
    "axes.edgecolor": "#444444",
    "axes.labelcolor": "#1A1A1A",
    "text.color": "#1A1A1A",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.color": "#444444", "ytick.color": "#444444",
    "xtick.labelcolor": "#1A1A1A", "ytick.labelcolor": "#1A1A1A",
    "xtick.major.size": 4, "ytick.major.size": 4,
    "xtick.minor.size": 2.5, "ytick.minor.size": 2.5,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "legend.frameon": False,
    "legend.fontsize": 9.5,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "figure.dpi": 130,
}


def apply_style() -> None:
    """Apply the shared rcParams. Call once at import time in each script."""
    matplotlib.use("Agg")
    plt.rcParams.update(RC)


def shade(hex_colour: str, factor: float) -> str:
    """Lighten (factor>1) or darken (factor<1) a hex colour."""
    h = hex_colour.lstrip("#")
    rgb = [int(h[i:i + 2], 16) for i in (0, 2, 4)]
    out = [max(0, min(255, int(round(c * factor)))) for c in rgb]
    return "#%02x%02x%02x" % tuple(out)


def panel_label(ax, s, xoff=-0.16, yoff=1.03, fontsize=12):
    ax.text(xoff, yoff, s, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom")


def grid(ax, axis="both", alpha=0.16):
    ax.grid(alpha=alpha, lw=0.5, axis=axis, color=CGRID)
