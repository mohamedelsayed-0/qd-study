"""Rebuild all paper figures in a single, journal-consistent style. v2."""
from __future__ import annotations
import os
import numpy as np
import scipy.io as sio
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (CD1, CD2, CACC, CPVK, CQD, CGRID, CANNOT,
                      CSHADE, apply_style, panel_label, grid)
apply_style()

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(REPO, "figures")
RES  = os.path.join(REPO, "results")


def save(name):
    plt.savefig(os.path.join(FIG, name)); plt.close(); print("wrote", name)

# ---- loaders -----------------------------------------------------------------
def load_jv(v):
    d = np.loadtxt(os.path.join(RES,"charge",f"final_{v}_jv.csv"),
                   delimiter=",", skiprows=1)
    return d[:,0], d[:,1]

def mpp(V, J):
    i = int(np.argmax(V*J)); return V[i], J[i]

def load_fdtd(v):
    p = os.path.join(RES,"fdtd", "planar_results.mat" if v=="D1" else "embedded_5nm_results.mat")
    with h5py.File(p,"r") as h:
        lam = np.array(h["lam"]).ravel()*1e9
        out = dict(lam=lam,
                   A=np.array(h["Atot"]).ravel(),
                   R=np.array(h["Rspec"]).ravel(),
                   T=np.array(h["Tspec"]).ravel())
        if v=="D1":
            out["Apv"] = np.array(h["Apvk"]).ravel()
            out["Aqd"] = np.array(h["Aqd"]).ravel()
        return out

# =============================================================================
def fig_jv():
    V1,J1 = load_jv("D1"); V2,J2 = load_jv("D2")
    d1m = np.loadtxt(os.path.join(RES,"charge","final_D1_max_jv.csv"), delimiter=",", skiprows=1)
    d2m = np.loadtxt(os.path.join(RES,"charge","final_D2_max_jv.csv"), delimiter=",", skiprows=1)
    fig,ax = plt.subplots(figsize=(5.2,3.8))
    ax.plot(V1,J1,color=CD1,lw=1.8,alpha=0.55,label="D1 planar")
    ax.plot(V2,J2,color=CD2,lw=1.8,alpha=0.55,label="D2 planar")
    ax.plot(d1m[:,0],d1m[:,1],color=CD1,lw=2.3,ls="--",label="D1 optimized")
    ax.plot(d2m[:,0],d2m[:,1],color=CD2,lw=2.3,ls="--",label="D2 optimized")
    ax.set_xlabel("Voltage (V)"); ax.set_ylabel("Current density (mA cm$^{-2}$)")
    ax.set_xlim(0,0.8); ax.set_ylim(0,53)
    tbl  = ("optimized      PCE      Jsc     Voc      FF\n"
            "D1          24.00 %   42.7   0.741   0.76\n"
            "D2          24.21 %   50.2   0.719   0.67")
    ax.text(0.03, 0.03, tbl, transform=ax.transAxes,
            fontsize=8.6, family="monospace", va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#ccc", lw=0.6))
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.98), fontsize=8.8)
    grid(ax)
    save("fig1_jv_final.png")

def fig_optics():
    d1 = load_fdtd("D1"); d2 = load_fdtd("D2")
    fig, axes = plt.subplots(1,2, figsize=(8.8,3.4))
    axA, axR = axes
    axA.plot(d1["lam"], d1["A"], color=CD1, lw=1.9, label="D1 — planar QD film")
    axA.plot(d2["lam"], d2["A"], color=CD2, lw=1.9, label="D2 — embedded QD")
    axA.fill_between(d1["lam"], 0, d1["A"], color=CD1, alpha=0.06)
    axA.set_xlim(300,1700); axA.set_ylim(0,1.02)
    axA.set_xlabel("Wavelength (nm)"); axA.set_ylabel("Absorptance")
    axA.legend(loc="lower left")
    panel_label(axA, "(a)"); grid(axA)

    axR.plot(d1["lam"], d1["R"], color=CD1, lw=1.9, label="D1")
    axR.plot(d2["lam"], d2["R"], color=CD2, lw=1.9, label="D2")
    axR.set_xlim(300,1700); axR.set_ylim(0,1.02)
    axR.set_xlabel("Wavelength (nm)"); axR.set_ylabel("Reflectance")
    axR.legend(loc="upper left", framealpha=0.0)
    panel_label(axR, "(b)"); grid(axR)
    plt.tight_layout()
    save("fig2_optics_RA.png")

def fig_interface():
    d = np.loadtxt(os.path.join(RES,"charge","interface_sensitivity.csv"),
                   delimiter=",", skiprows=1)
    srv = d[:,0]; p1 = d[:,1]; p2 = d[:,2]
    fig,ax = plt.subplots(figsize=(5.0,3.6))
    ax.axvspan(1e2, 1e3, color="#eaeaea", zorder=0)
    ax.text(3.16e2, 21.7, "passivated\n(literature)", fontsize=8.5, color="#666",
            ha="center", va="top")
    ax.axhline(20.59, color=CD1, ls=":", lw=0.9, alpha=0.7)
    ax.axhline(20.27, color=CD2, ls=":", lw=0.9, alpha=0.7)
    ax.text(1.03e6, 21.15, "ideal ($S=0$)", fontsize=8.5, color=CANNOT,
        ha="right", va="bottom")
    ax.plot(srv, p1, "-o", color=CD1, lw=2.0, ms=6, mec="white", mew=0.8, label="D1 — planar")
    ax.plot(srv, p2, "-s", color=CD2, lw=2.0, ms=6, mec="white", mew=0.8, label="D2 — embedded")
    ax.set_xscale("log")
    ax.set_xlabel("Interface recombination velocity $S$ (cm s$^{-1}$)")
    ax.set_ylabel("PCE (%)")
    ax.set_ylim(13, 22)
    ax.legend(loc="lower left")
    grid(ax)
    save("fig4_interface_sensitivity.png")

def fig_bandgap():
    d = np.loadtxt(os.path.join(RES,"charge","qd_bandgap_sweep.csv"),
                   delimiter=",", skiprows=1, usecols=(1,3,4,5,6))
    order = np.argsort(d[:,0]); d = d[order]
    Eg, Jsc, Voc, FF, PCE = d.T
    ib = int(np.argmax(PCE))

    fig, axes = plt.subplots(1,2, figsize=(8.4,3.5))
    ax = axes[0]
    ax.plot(Eg, PCE, "-o", color=CD1, lw=2.0, ms=6.5, mec="white", mew=0.8, zorder=2)
    ax.plot([Eg[ib]], [PCE[ib]], "*", color=CD2, ms=15, mec="white", mew=0.8, zorder=3)
    ax.annotate("optimum  1.03 eV\n20.4 %",
                xy=(Eg[ib], PCE[ib]),
                xytext=(1.11, 21.5), textcoords="data",
                fontsize=9, color=CD2, ha="left",
                arrowprops=dict(arrowstyle="-", color="#888", lw=0.7))
    ax.set_xlabel("QD bandgap (eV)"); ax.set_ylabel("PCE (%)")
    ax.set_ylim(12, 23); grid(ax); panel_label(ax, "(a)")

    ax = axes[1]
    ax.plot(Eg, Jsc/Jsc.max(), "-o", color=CQD, lw=2.0, ms=6, mec="white", mew=0.8, label="$J_\\mathrm{sc}$ (norm.)")
    ax.plot(Eg, Voc/Voc.max(), "-s", color=CD2, lw=2.0, ms=6, mec="white", mew=0.8, label="$V_\\mathrm{oc}$ (norm.)")
    ax.set_xlabel("QD bandgap (eV)"); ax.set_ylabel("Normalized to maximum")
    ax.set_ylim(0.55, 1.05); grid(ax)
    ax.legend(loc="lower center")
    panel_label(ax, "(b)")
    plt.tight_layout()
    save("fig5_qd_bandgap.png")

def fig_eqe():
    d = load_fdtd("D1")
    lam, Apv, Aqd, Atot = d["lam"], d["Apv"], d["Aqd"], d["A"]
    fig,ax = plt.subplots(figsize=(5.4,3.6))
    ax.fill_between(lam, 0, Apv, color=CPVK, alpha=0.22)
    ax.fill_between(lam, 0, Aqd, color=CQD, alpha=0.22)
    ax.plot(lam, Apv, color=CPVK, lw=2.0, label="MAPbI$_3$")
    ax.plot(lam, Aqd, color=CQD, lw=2.0, label="PbS-QD")
    ax.plot(lam, Atot, color=CACC, lw=1.3, ls="--", label="Total absorbed")
    # MAPbI3 cutoff from the modelled gap: 1239.842 / 1.55 eV = 800 nm.
    # (Must match scripts/spectral_filter.py, which uses the same relation.)
    ax.axvline(799.9, color=CGRID, lw=0.8, ls=":")
    ax.text(799.9, 1.02, " MAPbI$_3$ cutoff (800 nm)", fontsize=8.5, color=CANNOT,
            ha="left", va="bottom")
    ax.set_xlim(300,1700); ax.set_ylim(0,1.08)
    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Fraction absorbed per layer")
    ax.legend(loc="upper right")
    grid(ax)
    save("fig8_absorptance.png")

def fig_tornado():
    # Scans whose endpoints are bit-identical to the base value did not
    # propagate to the solver and are excluded rather than plotted as a
    # zero-width bar. This currently removes the D2 blend-gap row
    # (1.40-1.50 eV -> 20.27/20.27), which cannot be a physical result: the
    # bandgap sets ni and therefore Voc. Recorded in the CSV for provenance
    # and documented in the manuscript's excluded-configurations paragraph.
    rows = []
    with open(os.path.join(RES,"charge","sensitivity_oat.csv")) as f:
        next(f)
        for line in f:
            p = line.strip().split(",")
            if len(p) < 5 or p[3] in ("","None") or p[4] in ("","None"): continue
            base, lo, hi = float(p[1]), float(p[3]), float(p[4])
            if lo == hi == base:
                print("  excluded (no solver response): %s / %s" % (p[0], p[2]))
                continue
            rows.append((p[0], base, p[2], lo, hi))
    D1 = [r for r in rows if r[0]=="D1"]
    D2 = [r for r in rows if r[0]=="D2"]

    def draw(ax, data, color, ref_val, xlim, title):
        # most sensitive at top -> sort descending span, then reverse for barh order
        data = sorted(data, key=lambda r: abs(r[4]-r[3]))
        labels = [r[2] for r in data]
        lo = np.array([r[3] for r in data]); hi = np.array([r[4] for r in data])
        y = np.arange(len(data))
        left = np.minimum(lo,hi); width = np.abs(hi-lo)
        # The plotted references are the conservative baselines (D1 20.59,
        # D2 20.27), not the optimized devices.
        ax.axvline(ref_val, color="#555", lw=1.1, ls="--", zorder=1,
                   label="reference device")
        ax.barh(y, width, left=left, color=color, alpha=0.78,
                edgecolor=color, height=0.55, zorder=2)
        ax.plot(lo,y,"o",color=color,mec="white",mew=0.7,ms=5, zorder=3)
        ax.plot(hi,y,"o",color=color,mec="white",mew=0.7,ms=5, zorder=3)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
        ax.set_ylim(-0.7, len(data) - 0.15)
        ax.set_xlabel("PCE (%)")
        ax.set_xlim(*xlim)
        ax.set_title(title, fontsize=11, color=color, pad=8)
        ax.legend(loc="lower right", fontsize=8.5,
                  borderaxespad=0.3)
        grid(ax, axis="x")

    fig, axes = plt.subplots(1,2, figsize=(9.8,3.6))
    draw(axes[0], D1, CD1, D1[0][1], (16.5, 23), "D1 — planar QD film")
    draw(axes[1], D2, CD2, D2[0][1], (16.5, 23), "D2 — embedded QD")
    panel_label(axes[0], "(a)", xoff=-0.32, yoff=1.08)
    panel_label(axes[1], "(b)", xoff=-0.32, yoff=1.08)
    plt.tight_layout()
    save("fig10_tornado.png")

def fig_grating():
    rows = {}
    with open(os.path.join(RES, "charge", "grating_lighttrapping_full.csv")) as f:
        next(f)
        for line in f:
            p = line.strip().split(",")
            if len(p) == 6:
                rows[p[0]] = [float(v) for v in p[1:]]
    # cols: Jopt, Jsc, Voc, FF, PCE
    def prof(gp):
        d = sio.loadmat(gp); x=d["G_x"].ravel(); y=d["G_y"].ravel(); z=d["G_z"].ravel()
        raw=np.asarray(d["G_3d"],float); exp=(len(x),len(y),len(z))
        G = raw if raw.shape==exp else raw.transpose(2,1,0)
        p = np.trapezoid(np.trapezoid(G, x=x, axis=0), x=y, axis=0)/((x[-1]-x[0])*(y[-1]-y[0]))
        return z*1e9, p
    # Use committed reference depth profiles (rebuildable via
    # scripts/build_grating_fdtd.py). If either file is missing, refuse to
    # draw the panel with a placeholder — regenerate the FDTD first.
    flat_ref = os.path.join(RES, "fdtd", "G_d1_flat_ref.mat")
    grat_ref = os.path.join(RES, "fdtd", "G_d1_max.mat")
    if not (os.path.isfile(flat_ref) and os.path.isfile(grat_ref)):
        raise SystemExit(
            "fig11 requires committed FDTD profiles; run "
            "`python scripts/build_grating_fdtd.py d1_flat_ref` and "
            "`python scripts/build_grating_fdtd.py d1_max` first")
    zf, pf = prof(flat_ref)
    zg, pg = prof(grat_ref)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    ax.plot(zf, pf/1e27, color="#666", lw=2.0, label="Flat Au back")
    ax.plot(zg, pg/1e27, color=CD2, lw=2.0, label="Textured Au grating")
    ax.set_xlabel("Depth from back contact (nm)")
    ax.set_ylabel("Generation rate ($10^{27}$ m$^{-3}$ s$^{-1}$)")
    ax.legend(loc="upper left"); grid(ax)
    panel_label(ax, "(a)")

    ax = axes[1]
    cfgs = ["Flat", "Grating"]
    jsc = [40.03, rows["P600_t200"][1]]
    pce = [21.14, rows["P600_t200"][4]]
    x = np.arange(2); w = 0.32
    ax.bar(x-w/2, jsc, w, color=CD1, edgecolor="none")
    ax.set_ylabel("$J_\\mathrm{sc}$ (mA cm$^{-2}$)", color=CD1)
    ax.tick_params(axis='y', colors=CD1)
    ax.set_ylim(38, 44); ax.set_xticks(x); ax.set_xticklabels(cfgs)
    for xi, v in zip(x-w/2, jsc):
        ax.text(xi, v+0.08, f"{v:.1f}", ha="center", fontsize=9, color=CD1)
    ax2 = ax.twinx()
    ax2.spines['right'].set_visible(True); ax2.spines['top'].set_visible(False)
    ax2.bar(x+w/2, pce, w, color=CD2, edgecolor="none")
    ax2.set_ylabel("PCE (%)", color=CD2)
    ax2.tick_params(axis='y', colors=CD2)
    ax2.set_ylim(20, 23)
    for xi, v in zip(x+w/2, pce):
        ax2.text(xi, v+0.04, f"{v:.2f}", ha="center", fontsize=9, color=CD2)
    ax.grid(alpha=0.14, axis="y", lw=0.5)
    panel_label(ax, "(b)")
    plt.tight_layout()
    save("fig11_grating.png")

# =============================================================================
def clean_supplementary():
    for name in ("fig3_generation.png","fig6_band_diagram.png","fig7_d2_fill.png","fig9_energy_check.png"):
        p = os.path.join(FIG, name)
        if os.path.isfile(p): os.remove(p); print("removed", name)

if __name__ == "__main__":
    clean_supplementary()
    fig_jv(); fig_optics(); fig_interface(); fig_bandgap(); fig_eqe(); fig_tornado(); fig_grating()
    print("done")
