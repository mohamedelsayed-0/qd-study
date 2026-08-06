"""Matched-QD-loading planar (D1) controls.

Fairness check requested in review: the headline D1 carries a 300 nm continuous
PbS-QD film, whereas D2 at 20 / 25 vol.% carries only 170 / 212.5 nm of
QD-equivalent material inside the same 850 nm active layer. Placement is
therefore not the *only* changed variable -- loading differs too. This script
builds planar D1 devices whose QD volume is matched to D2 while holding the
850 nm active thickness fixed:

    t_qd = 170.0 nm  ->  matches D2 at 20 vol.%   (t_pvk = 680.0 nm)
    t_qd = 212.5 nm  ->  matches D2 at 25 vol.%   (t_pvk = 637.5 nm)

For each it produces, through the *committed* pipeline (same geometry builder,
same FDTD->CHARGE conversion as build_grating_fdtd.py, same above-gap estimator
as spectral_filter.py), both a raw and a spectrally filtered CHARGE-ready
generation file:

    data/generation/charge_d1_load{TAG}.mat            (raw)
    data/generation/charge_d1_load{TAG}_filtered.mat   (strict single-gap)

plus a summary row in results/fdtd/matched_loading.csv. The CHARGE solves are
run separately via the d1_load* presets added to scripts/charge/presets.py.

Flat stack, mesh accuracy 3 (matching final_device.lsf, i.e. the committed
baseline D1 charge_hedge). Headless; writes each fill before the next so a
killed run keeps completed fills.
"""
from __future__ import annotations

import os, sys, shutil, subprocess, tempfile, time
import numpy as np
import scipy.io as sio

c = 299792458.0
Q = 1.602176634e-19
REPO = r"C:\Users\alkin\OneDrive\Documents\qd-research\Quantom-dot-solar"
BASE = os.path.join(REPO, "models", "fdtd", "final_planar_double_junction.fsp")
COMMON = os.path.join(REPO, "scripts", "common").replace("\\", "/")
GENDIR = os.path.join(REPO, "data", "generation")
OUT = os.path.join(REPO, "results", "fdtd", "matched_loading.csv")

# QD gap 1.03 eV (D1 film), MAPbI3 gap 1.55 eV -- same edges as spectral_filter.
GAP_QD = 1.03
GAP_PVK = 1.55
# active layer is fixed at 850 nm = t_qd + t_pvk; QD film sits 160 nm above Au.
ACTIVE_NM = 850.0
QD_BOTTOM_NM = 160.0
CASES = [("170", 170.0), ("212p5", 212.5)]   # tag, t_qd (nm)

sys.path.append(r"C:\Program Files\Lumerical\v261\api\python")
sys.path.insert(0, os.path.join(REPO, "scripts"))
import lumapi
from spectral_filter import above_gap_fraction   # exact committed estimator


def cleanup():
    for exe in ("fdtd-engine-msmpi.exe", "fdtd-solutions.exe"):
        subprocess.run(["taskkill", "/F", "/IM", exe], capture_output=True)
    time.sleep(2)


def build_and_run(f, t_qd_nm):
    """Rebuild the planar stack at t_qd, run FDTD, return (z, prof, lam, Apvk, Aqd)."""
    t_pvk_nm = ACTIVE_NM - t_qd_nm
    f.eval(
        'addpath("%s"); '
        'span=50e-9; qd_mat="QD_1550_NH4I"; '
        't_au=100e-9; t_mox=10e-9; t_edt=50e-9; '
        't_qd=%ge-9; t_pvk=%ge-9; t_sno2=40e-9; t_fto=300e-9; '
        'build_planar_double_junction; '
        'setnamed("FDTD","mesh accuracy",3); '
        'setnamed("FDTD","simulation time",1200e-15); '
        'gen_x=0; gen_x_span=span; gen_y=0; gen_y_span=span; '
        'gen_z=(z_qd_bottom+z_pvk_top)/2; gen_z_span=z_pvk_top-z_qd_bottom; '
        'prepare_generation_group;' % (COMMON, t_qd_nm, t_pvk_nm)
    )
    try:
        f.setresource("FDTD", 1, "processes", 1); f.setresource("FDTD", 1, "threads", 4)
    except Exception:
        pass
    print("  running FDTD acc 3, t_qd=%g nm, t_pvk=%g nm" % (t_qd_nm, t_pvk_nm))
    f.run()
    f.runanalysis("solar_generation")

    G = f.getresult("solar_generation", "G")
    x = np.asarray(G["x"]).ravel(); y = np.asarray(G["y"]).ravel(); z = np.asarray(G["z"]).ravel()
    raw = np.asarray(G["G"], dtype=float)
    g = raw if raw.shape == (len(x), len(y), len(z)) else raw.transpose(2, 1, 0)
    prof = np.trapezoid(np.trapezoid(g, x=x, axis=0), x=y, axis=0) / ((x[-1] - x[0]) * (y[-1] - y[0]))

    # layer-resolved absorptance for the above-gap fractions (final_device.lsf)
    f.eval(
        'fmon=getdata("P_pvk_qd","f"); lam=c/fmon; '
        'Apvk=(-transmission("P_pvk_top"))-(-transmission("P_pvk_qd")); '
        'Aqd=(-transmission("P_pvk_qd"))-(-transmission("P_qd_edt"));'
    )
    lam = np.asarray(f.getv("lam"), dtype=float).ravel()
    Apvk = np.asarray(f.getv("Apvk"), dtype=float).ravel()
    Aqd = np.asarray(f.getv("Aqd"), dtype=float).ravel()
    return z, prof, lam, Apvk, Aqd


def save_charge(path, z, prof):
    Gc = np.broadcast_to(prof[None, :, None], (2, len(z), 2)).copy()
    sio.savemat(path, {"x": np.array([[-25e-9], [25e-9]]), "y": z[:, None],
                       "z": np.array([[-0.5e-6], [0.5e-6]]), "G": Gc},
                do_compression=True)


def solve_case(tag, t_qd_nm):
    tmp = tempfile.mkdtemp(prefix="qd_match_")
    work = os.path.join(tmp, "work.fsp")
    shutil.copy2(BASE, work); cleanup()
    f = lumapi.FDTD(hide=True)
    try:
        f.load(work); f.switchtolayout()
        z, prof, lam, Apvk, Aqd = build_and_run(f, t_qd_nm)
    finally:
        try: f.close()
        except Exception: pass
        cleanup(); shutil.rmtree(tmp, ignore_errors=True)

    jopt_raw = float(Q * np.trapezoid(prof, x=z) / 10.0)
    agf_qd = above_gap_fraction(lam * 1e9, Aqd, GAP_QD)
    agf_pvk = above_gap_fraction(lam * 1e9, Apvk, GAP_PVK)

    # spectral filter with this geometry's boundaries + fresh fractions
    z_nm = z * 1e9
    qd_top_nm = QD_BOTTOM_NM + t_qd_nm
    scale = np.ones_like(prof)
    scale[(z_nm >= QD_BOTTOM_NM) & (z_nm < qd_top_nm)] = agf_qd
    scale[(z_nm >= qd_top_nm) & (z_nm <= QD_BOTTOM_NM + ACTIVE_NM)] = agf_pvk
    prof_filt = prof * scale
    jopt_filt = float(Q * np.trapezoid(prof_filt, x=z) / 10.0)

    save_charge(os.path.join(GENDIR, "charge_d1_load%s.mat" % tag), z, prof)
    save_charge(os.path.join(GENDIR, "charge_d1_load%s_filtered.mat" % tag), z, prof_filt)
    print("  t_qd=%g nm  Jopt_raw=%.3f  agf_qd=%.4f  agf_pvk=%.4f  Jopt_filt=%.3f (%.3f)"
          % (t_qd_nm, jopt_raw, agf_qd, agf_pvk, jopt_filt, jopt_filt / jopt_raw))
    return jopt_raw, agf_qd, agf_pvk, jopt_filt


if __name__ == "__main__":
    tags = sys.argv[1:] or [t for t, _ in CASES]
    if not os.path.isfile(OUT):
        with open(OUT, "w") as fh:
            fh.write("tag,t_qd_nm,t_pvk_nm,Jopt_raw_mA_cm2,agf_qd_1p03,"
                     "agf_mapbi3_1p55,Jopt_filtered_mA_cm2,retained_fraction\n")
    for tag, t_qd_nm in CASES:
        if tag not in tags:
            continue
        print("== case %s (t_qd=%g nm) ==" % (tag, t_qd_nm))
        jr, aq, ap, jf = solve_case(tag, t_qd_nm)
        with open(OUT, "a") as fh:
            fh.write("%s,%g,%g,%.3f,%.4f,%.4f,%.3f,%.4f\n"
                     % (tag, t_qd_nm, ACTIVE_NM - t_qd_nm, jr, aq, ap, jf, jf / jr))
    print("done ->", OUT)
