"""Layer-resolved above-gap / below-gap current budget.

Every reported Jsc is compared against the ideal single-gap ceiling
(only photons with energy above the electrical bandgap of the absorbing
layer are counted). Uses the committed FDTD absorptance spectra and
ASTM G173 direct+circumsolar (AM1.5G) photon flux.

Output:
    results/fdtd/energy_budget.csv - per-design totals
    stdout - human-readable table
"""
import os
import numpy as np
import h5py

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q = 1.602176634e-19
HC = 1.98644586e-25       # J m
# ASTM G173 global tilt integrated over 300-1700 nm is 945.62 W/m^2 (of 1000.06 W/m^2 full band)
AM15G_300_1700_W = 945.62

# ASTM G173 direct+circumsolar reference (wavelength nm, spectral irradiance W/m2/nm)
# coarse ~50 nm resolution table sufficient for gap integrals
_G173 = np.array([
    [300, 0.00], [350, 0.30], [400, 0.65], [450, 1.65], [500, 1.85],
    [550, 1.80], [600, 1.60], [650, 1.55], [700, 1.42], [750, 1.20],
    [800, 1.08], [850, 0.94], [900, 0.79], [950, 0.60], [1000, 0.75],
    [1050, 0.72], [1100, 0.60], [1150, 0.36], [1200, 0.50], [1250, 0.42],
    [1300, 0.25], [1350, 0.05], [1400, 0.03], [1450, 0.35], [1500, 0.32],
    [1550, 0.27], [1600, 0.24], [1650, 0.20], [1700, 0.16],
])
# renormalize the coarse table so its integral matches the official 1000.06 W/m2
_lam = _G173[:, 0]
_S = _G173[:, 1]
_S = _S * (AM15G_300_1700_W / np.trapezoid(_S, _lam))


def photon_flux(lam_nm):
    """Photon flux density S*lam/hc, at requested wavelengths (nm)."""
    S = np.interp(lam_nm, _lam, _S)   # W/m2/nm at requested lam
    return S * (lam_nm * 1e-9) / HC   # photons/m2/s/nm


def load_absorptance(mat):
    with h5py.File(mat, "r") as h:
        lam = np.array(h["lam"]).ravel() * 1e9  # nm
        keys = {k: np.array(h[k]).ravel() for k in h.keys() if k != "lam"}
    return lam, keys


def j_from_A(lam_nm, A):
    flux = photon_flux(lam_nm)
    order = np.argsort(lam_nm)
    lam_s = lam_nm[order]
    A_s = np.clip(A[order], 0, 1)
    flux_s = flux[order]
    return Q * np.trapezoid(A_s * flux_s, lam_s) / 10.0   # mA/cm2


def split_above_below(lam_nm, A, gap_eV):
    lam_gap = 1239.842 / gap_eV
    m_above = lam_nm <= lam_gap
    lam_a = lam_nm[m_above]; A_a = A[m_above]
    lam_b = lam_nm[~m_above]; A_b = A[~m_above]
    return j_from_A(lam_a, A_a), j_from_A(lam_b, A_b), lam_gap


def ideal_ceiling(gap_eV, lam_max_nm=1700):
    lam = np.linspace(300, lam_max_nm, 4001)
    lam_gap = 1239.842 / gap_eV
    A = (lam <= lam_gap).astype(float)
    return j_from_A(lam, A), lam_gap


D1_pvk_gap = 1.55
D1_qd_gap = 1.03
D2_blend_gap = 1.44
D2_qd_first_exciton = 1.03  # underlying dot excitonic edge for reporting

d1_lam, d1 = load_absorptance(os.path.join(REPO, "results", "fdtd", "planar_results.mat"))
d2_lam, d2 = load_absorptance(os.path.join(REPO, "results", "fdtd", "embedded_5nm_results.mat"))

rows = []
print("=" * 88)
print("Layer-resolved above/below electrical-gap current  (mA/cm^2)")
print("=" * 88)

d1_pvk_above, d1_pvk_below, _ = split_above_below(d1_lam, d1["Apvk"], D1_pvk_gap)
d1_qd_above, d1_qd_below, _ = split_above_below(d1_lam, d1["Aqd"], D1_qd_gap)
d1_pvk_ceiling, _ = ideal_ceiling(D1_pvk_gap)
d1_qd_ceiling, _ = ideal_ceiling(D1_qd_gap)
print(f"D1 perovskite (Eg=1.55 eV, lam_gap 800 nm)")
print(f"  above: {d1_pvk_above:6.2f}   below: {d1_pvk_below:6.2f}   ideal-ceiling: {d1_pvk_ceiling:6.2f}")
print(f"D1 QD film  (Eg=1.03 eV, lam_gap 1204 nm)")
print(f"  above: {d1_qd_above:6.2f}   below: {d1_qd_below:6.2f}   ideal-ceiling: {d1_qd_ceiling:6.2f}")
d1_total_above = d1_pvk_above + d1_qd_above
d1_total_below = d1_pvk_below + d1_qd_below
print(f"D1 total active-layer: above {d1_total_above:.2f}  below {d1_total_below:.2f}  "
      f"below/total {100*d1_total_below/(d1_total_above+d1_total_below):.1f}%")
print()

# D2 blend: split against 1.44 eV effective gap; also report against 1.03 eV QD gap
d2_blend_above, d2_blend_below_144, _ = split_above_below(d2_lam, d2["Aactive"], D2_blend_gap)
_, d2_blend_below_103, _ = split_above_below(d2_lam, d2["Aactive"], D2_qd_first_exciton)
d2_ceiling_144, _ = ideal_ceiling(D2_blend_gap)
d2_ceiling_103, _ = ideal_ceiling(D2_qd_first_exciton)
d2_total = d2_blend_above + d2_blend_below_144
print(f"D2 blend  (Eg=1.44 eV, lam_gap 861 nm)")
print(f"  above: {d2_blend_above:6.2f}   below: {d2_blend_below_144:6.2f}   ideal-ceiling: {d2_ceiling_144:6.2f}")
print(f"D2 blend also split at 1.03 eV underlying dot excitonic edge")
print(f"  above 1.03 eV: {d2_total - d2_blend_below_103:6.2f}   below 1.03 eV: {d2_blend_below_103:6.2f}   ideal-ceiling: {d2_ceiling_103:6.2f}")
print(f"D2 subgap share of total: {100*d2_blend_below_144/d2_total:.1f}%")
print()

# Correction factors for reported Jsc values
D1_reported_Jsc_baseline = 39.07
D1_reported_Jsc_optimized = 42.65
D2_reported_Jsc_baseline = 42.23
D2_reported_Jsc_optimized = 50.17
D2_reported_Jsc_hiV = 36.92

def correct(jsc_reported, above, total):
    return jsc_reported * above / total

d1_frac = d1_total_above / (d1_total_above + d1_total_below)
d2_frac = d2_blend_above / d2_total

print("Single-gap-consistent (above-gap only) Jsc estimate for reported devices")
print(f"Rescaling factor: D1 {d1_frac:.3f}   D2 {d2_frac:.3f}")
print()
print(f"D1 planar baseline:   reported {D1_reported_Jsc_baseline:.2f}  -> single-gap {correct(D1_reported_Jsc_baseline, d1_total_above, d1_total_above+d1_total_below):.2f}")
print(f"D1 optimized:         reported {D1_reported_Jsc_optimized:.2f}  -> single-gap {correct(D1_reported_Jsc_optimized, d1_total_above, d1_total_above+d1_total_below):.2f}")
print(f"D2 planar baseline:   reported {D2_reported_Jsc_baseline:.2f}  -> single-gap {correct(D2_reported_Jsc_baseline, d2_blend_above, d2_total):.2f}")
print(f"D2 optimized:         reported {D2_reported_Jsc_optimized:.2f}  -> single-gap {correct(D2_reported_Jsc_optimized, d2_blend_above, d2_total):.2f}")
print(f"D2 high-voltage:      reported {D2_reported_Jsc_hiV:.2f}   -> single-gap {correct(D2_reported_Jsc_hiV, d2_blend_above, d2_total):.2f}")

out = os.path.join(REPO, "results", "fdtd", "energy_budget.csv")
with open(out, "w") as f:
    f.write("design,layer,Eg_eV,J_above_gap_mA_cm2,J_below_gap_mA_cm2,ideal_ceiling_mA_cm2\n")
    f.write(f"D1,MAPbI3,1.55,{d1_pvk_above:.3f},{d1_pvk_below:.3f},{d1_pvk_ceiling:.3f}\n")
    f.write(f"D1,PbS-QD,1.03,{d1_qd_above:.3f},{d1_qd_below:.3f},{d1_qd_ceiling:.3f}\n")
    f.write(f"D2,blend,1.44,{d2_blend_above:.3f},{d2_blend_below_144:.3f},{d2_ceiling_144:.3f}\n")
    f.write(f"D2,blend_at_1p03,1.03,{d2_total - d2_blend_below_103:.3f},{d2_blend_below_103:.3f},{d2_ceiling_103:.3f}\n")
print(f"\nWrote {out}")
