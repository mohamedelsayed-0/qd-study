"""Produce spectrally filtered generation profiles for a strict single-gap CHARGE rerun.

The committed FDTD generation profiles integrate absorption at all wavelengths,
including photons below each absorber's electrical bandgap. Under the strict
single-gap assumption, sub-gap absorption cannot create collectable carriers.

This script:
  1) Reads the layer-resolved absorptance spectra from the committed FDTD
     results (planar_results.mat, embedded_5nm_results.mat).
  2) Computes the per-layer above-gap fraction of Jopt using an ASTM G173
     photon flux (300-1700 nm, integrated to 945.62 W/m^2 per NREL).
  3) Reads each committed CHARGE generation file, splits G(y) by depth
     using the known stack boundaries, and rescales each layer segment by
     its above-gap fraction.
  4) Writes filtered generation files data/generation/charge_<name>_filtered.mat
     which can be consumed by the standard build_planar.py / build_embedded.py.

Layer boundaries (nm) used by both FDTD and CHARGE:
  Au         0 - 100     (back electrode)
  MoOx     100 - 110
  EDT      110 - 160
  QD-Pbs   160 - 460     (D1)      -> Eg = 1.03 eV
  MAPbI3   460 - 1010    (D1)      -> Eg = 1.55 eV
  blend    160 - 1010    (D2)      -> Eg = 1.44 eV
  SnO2    1010 - 1050
  FTO     1050 - 1350
"""
from __future__ import annotations
import os
import numpy as np
import scipy.io as sio
import h5py

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q = 1.602176634e-19
HC = 1.98644586e-25
AM15G_300_1700_W = 945.62

_G173 = np.array([
    [300, 0.00], [350, 0.30], [400, 0.65], [450, 1.65], [500, 1.85],
    [550, 1.80], [600, 1.60], [650, 1.55], [700, 1.42], [750, 1.20],
    [800, 1.08], [850, 0.94], [900, 0.79], [950, 0.60], [1000, 0.75],
    [1050, 0.72], [1100, 0.60], [1150, 0.36], [1200, 0.50], [1250, 0.42],
    [1300, 0.25], [1350, 0.05], [1400, 0.03], [1450, 0.35], [1500, 0.32],
    [1550, 0.27], [1600, 0.24], [1650, 0.20], [1700, 0.16],
])
_S_lam = _G173[:, 0]
_S = _G173[:, 1] * (AM15G_300_1700_W / np.trapezoid(_G173[:, 1], _G173[:, 0]))


def photon_flux(lam_nm):
    return np.interp(lam_nm, _S_lam, _S) * (lam_nm * 1e-9) / HC


def above_gap_fraction(lam_nm, A, gap_eV):
    """Fraction of layer's Jopt that comes from photons with energy > gap."""
    order = np.argsort(lam_nm)
    lam = lam_nm[order]
    A = np.clip(A[order], 0, 1)
    flux = photon_flux(lam)
    integrand = A * flux
    tot = np.trapezoid(integrand, lam)
    lam_gap = 1239.842 / gap_eV
    above = lam <= lam_gap
    top = np.trapezoid(integrand[above], lam[above])
    return float(top / tot) if tot > 0 else 0.0


def per_layer_fractions():
    """Return dict of layer -> above-gap fraction for D1 and D2 blends."""
    with h5py.File(os.path.join(REPO, "results", "fdtd", "planar_results.mat"), "r") as h:
        lam = np.array(h["lam"]).ravel() * 1e9
        Apvk = np.array(h["Apvk"]).ravel()
        Aqd = np.array(h["Aqd"]).ravel()
    with h5py.File(os.path.join(REPO, "results", "fdtd", "embedded_5nm_results.mat"), "r") as h:
        lam2 = np.array(h["lam"]).ravel() * 1e9
        Aa = np.array(h["Aactive"]).ravel()
    return {
        "D1_QD": above_gap_fraction(lam, Aqd, 1.03),
        "D1_MAPbI3": above_gap_fraction(lam, Apvk, 1.55),
        "D2_blend": above_gap_fraction(lam2, Aa, 1.44),
    }


# Depth boundaries (y in meters) matching the CHARGE geometry
BOUNDARIES = {
    "D1": {
        "QD": (160e-9, 460e-9, "D1_QD"),
        "MAPbI3": (460e-9, 1010e-9, "D1_MAPbI3"),
    },
    "D2": {
        "blend": (160e-9, 1010e-9, "D2_blend"),
    },
}


def filter_charge_mat(src_path, out_path, design, fractions):
    d = sio.loadmat(src_path)
    y = d["y"].ravel()
    G = np.asarray(d["G"], dtype=float)     # shape (2, Ny, 2)
    # Grab the middle-depth profile (identical across x and z by construction)
    prof = G[0, :, 0].copy()
    layers = BOUNDARIES[design]
    scale = np.ones_like(prof)
    for name, (y0, y1, frac_key) in layers.items():
        f = fractions[frac_key]
        mask = (y >= y0) & (y <= y1)
        scale[mask] = f
    prof_filtered = prof * scale
    # write back into the same 3D shape CHARGE expects
    Gc = np.broadcast_to(prof_filtered[None, :, None], (2, len(y), 2)).copy()
    sio.savemat(out_path,
                {"x": d["x"], "y": d["y"], "z": d["z"], "G": Gc},
                do_compression=True)
    jopt_before = Q * np.trapezoid(prof, y) / 10.0
    jopt_after = Q * np.trapezoid(prof_filtered, y) / 10.0
    return jopt_before, jopt_after


if __name__ == "__main__":
    fractions = per_layer_fractions()
    print("Above-gap fractions:")
    for k, v in fractions.items():
        print(f"  {k}: {v:.3f}")
    print()

    plan = [
        # (source generation file, design, filtered output tag)
        ("charge_hedge.mat",       "D1", "charge_hedge_filtered.mat"),
        ("charge_embedded.mat",    "D2", "charge_embedded_filtered.mat"),
        ("charge_d1_max.mat",      "D1", "charge_d1_max_filtered.mat"),
        ("charge_d2_max.mat",      "D2", "charge_d2_max_filtered.mat"),
        ("charge_d1_max_unpol.mat","D1", "charge_d1_max_unpol_filtered.mat"),
        ("charge_d2_max_unpol.mat","D2", "charge_d2_max_unpol_filtered.mat"),
    ]

    gdir = os.path.join(REPO, "data", "generation")
    rows = []
    for src, design, out in plan:
        src_p = os.path.join(gdir, src)
        if not os.path.isfile(src_p):
            print(f"SKIP {src} (not found)")
            continue
        out_p = os.path.join(gdir, out)
        before, after = filter_charge_mat(src_p, out_p, design, fractions)
        rows.append((src, design, before, after))
        print(f"{src}: Jopt {before:.2f} -> {after:.2f} mA/cm2 (factor {after/before:.3f})")

    csv_out = os.path.join(REPO, "results", "fdtd", "spectral_filter_summary.csv")
    with open(csv_out, "w") as f:
        f.write("source_generation,design,Jopt_before,Jopt_after_filter,retained_fraction\n")
        for src, design, before, after in rows:
            f.write(f"{src},{design},{before:.3f},{after:.3f},{after/before:.4f}\n")
    print(f"\nWrote {csv_out}")
