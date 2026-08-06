"""Han 2018 absorber-level validation.

Han et al. (Small 14, 1801016, 2018) fabricated a PbS-in-MAPbI3 hybrid absorber
(4 nm dots, 0.5 mg/mL ~ 0.033 vol%) and measured (i) UV-vis absorption extending
modestly toward 900 nm and (ii) an IPCE showing no significant contribution
beyond ~820 nm -- i.e. the near-infrared absorption the dots add is not
collected as current.

This script reproduces the *optical* half of that observation with the present
FDTD workflow, using the in-repo measured 4 nm PbS n/k (matching Han's dot size)
and MAPbI3 host, at Han's loadings. For each fill it reports:
  - the active-layer optical current Jopt,
  - the above-gap fraction of that current (photons above the MAPbI3 gap,
    800 nm), computed with the same estimator as scripts/spectral_filter.py.

The prediction under test: at Han's loading the dots add negligible optical
current, and what little they add sits below the 800 nm collection edge -- so
the strict single-gap reading discards it, matching Han's IPCE. Pure MAPbI3 and
our 20 vol% blend are included as anchors.

Headless. Writes results/fdtd/han_validation.csv incrementally so a killed run
keeps completed rows. One FDTD solve per fill at mesh accuracy 2.
"""
from __future__ import annotations

import os, sys, shutil, subprocess, tempfile, time
import numpy as np

c = 299792458.0
Q = 1.602176634e-19
H = 6.62607015e-34
REPO = r"C:\Users\alkin\OneDrive\Documents\qd-research\Quantom-dot-solar"
BASE = os.path.join(REPO, "models", "fdtd", "final_planar_double_junction.fsp")
QD4 = os.path.join(REPO, "data", "optical", "qd", "by_size", "PbS-QD_4nm.txt")
OUT = os.path.join(REPO, "results", "fdtd", "han_validation.csv")
GAP_MAPBI3 = 1.55            # eV -> 800 nm collection edge
# All cases go through the same Maxwell-Garnett sampled-medium path so that the
# pure reference and the blends differ only by dot content, not by native-vs-
# resampled material (which alone costs ~2 mA/cm2 and would bias the trend).
# fill = 1e-6 is the f->0 pure reference; 3.3e-4 / 6.6e-4 are Han's 0.5 / 1
# mg/mL loadings; 0.20 is this work's anchor and the method check against the
# committed D2 above-gap fraction (0.591).
FILLS = [1e-6, 3.3e-4, 6.6e-4, 0.20]

sys.path.append(r"C:\Program Files\Lumerical\v261\api\python")
import lumapi

OBJS = ["Glass","FTO","SnO2","Perovskite","QD-Pbs","p-PbS-EDT","MOo3","Au","FDTD",
        "PlaneWave","R","T","P_pvk_top","P_pvk_qd","P_qd_edt","solar_generation"]


def photon_flux(lam_nm):
    """AM1.5G photon flux on the coarse table used by spectral_filter.py."""
    tab = np.array([
        [300,0.02],[350,0.20],[400,0.62],[450,1.05],[500,1.20],[550,1.28],
        [600,1.30],[650,1.28],[700,1.22],[750,1.15],[800,1.08],[850,0.94],
        [900,0.79],[950,0.60],[1000,0.75],[1050,0.99],[1100,0.75],[1150,0.49],
        [1200,0.65],[1250,0.60],[1300,0.55],[1350,0.13],[1400,0.05],[1450,0.20],
        [1500,0.39],[1550,0.44],[1600,0.42],[1650,0.39],[1700,0.34]])
    S = np.interp(lam_nm, tab[:,0], tab[:,1])          # W m-2 nm-1 (coarse)
    return S * (lam_nm*1e-9) / (H*c)                    # photons per unit


def above_gap_fraction(lam_nm, A, gap_eV):
    order = np.argsort(lam_nm)
    lam = lam_nm[order]; A = np.clip(A[order], 0, 1)
    integrand = A * photon_flux(lam)
    tot = np.trapezoid(integrand, lam)
    above = lam <= (1239.842/gap_eV)
    top = np.trapezoid(integrand[above], lam[above])
    return float(top/tot) if tot > 0 else 0.0


def mk_blend(f, fill):
    hostmat = f.getnamed("Perovskite", "material")
    Nw = 400
    lam = np.linspace(300e-9, 1700e-9, Nw); fg = c/lam
    # One vector call: getfdtdindex returns an (Nw,1) complex index array in the
    # same order as fg. (A per-point scalar call returns a (1,1) array, which
    # newer NumPy refuses to cast with complex().)
    nh = np.asarray(f.getfdtdindex(hostmat, fg, float(fg.min()), float(fg.max()))).ravel()
    rows = []
    for ln in open(QD4):
        p = ln.split()
        if len(p) == 3:
            try: rows.append([float(v) for v in p])
            except ValueError: pass
    dq = np.array(rows)
    nq = np.interp(lam, dq[:,0]*1e-9, dq[:,1]) + 1j*np.interp(lam, dq[:,0]*1e-9, dq[:,2])
    eps_h = nh**2; eps_i = nq**2
    beta = fill*(eps_i-eps_h)/(eps_i+2*eps_h)
    eps_eff = eps_h*(1+2*beta)/(1-beta)
    sd = np.zeros((Nw,2), dtype=complex); sd[:,0]=fg; sd[:,1]=eps_eff
    f.eval("if(materialexists('QD_han')==0){ m=addmaterial('Sampled 3D data'); setmaterial(m,'name','QD_han'); }")
    f.putv("sdh", sd); f.eval("setmaterial('QD_han','sampled data',sdh);")


def cleanup():
    for exe in ("fdtd-engine-msmpi.exe","fdtd-solutions.exe"):
        subprocess.run(["taskkill","/F","/IM",exe], capture_output=True)
    time.sleep(3)


def active_absorptance(f):
    """Active-layer absorptance and its optical current, using the exact
    power-box difference and solar weighting of the committed pipeline
    (scripts/final_embedded_5nm_effective.lsf):

        Aactive = -transmission("P_pvk_top") + transmission("P_qd_edt")

    Returns (lam_nm, A_active, Iam, Jactive_opt).
    """
    f.eval(
        'fmon=getdata("P_pvk_top","f"); lam=c/fmon; Nlam=length(lam); '
        'Aactive=-transmission("P_pvk_top")+transmission("P_qd_edt"); '
        'sw=solar(0); si=solar(1); Iam=matrix(Nlam,1); '
        'for(w_i=1:Nlam){ Iam(w_i)=interp(si,sw,lam(w_i)); } '
        'Jint=0; '
        'for(w_i=1:Nlam-1){ dl=lam(w_i+1)-lam(w_i); '
        '  Jint=Jint+0.5*((Aactive(w_i)*Iam(w_i)*lam(w_i)'
        '    +Aactive(w_i+1)*Iam(w_i+1)*lam(w_i+1))/(h*c))*dl; } '
        'Jact=abs(1.602176634e-19*Jint*0.1);'
    )
    lam = np.asarray(f.getv("lam"), dtype=float).ravel() * 1e9
    A = np.asarray(f.getv("Aactive"), dtype=float).ravel()
    Iam = np.asarray(f.getv("Iam"), dtype=float).ravel()
    Jact = float(np.asarray(f.getv("Jact")).ravel()[0])
    return lam, A, Iam, Jact


def above_gap_fraction_weighted(lam_nm, A, Iam, gap_eV):
    """Above-gap fraction using the model's own solar irradiance Iam."""
    order = np.argsort(lam_nm)
    lam = lam_nm[order]; A = np.clip(A[order], 0, 1); Iam = Iam[order]
    integrand = A * Iam * lam          # photon-current weighting ~ A*I*lam
    tot = np.trapezoid(integrand, lam)
    above = lam <= (1239.842/gap_eV)
    top = np.trapezoid(integrand[above], lam[above])
    return float(top/tot) if tot > 0 else 0.0


def solve(fill):
    tmp = tempfile.mkdtemp(prefix="qd_han_")
    work = os.path.join(tmp, "work.fsp")
    shutil.copy2(BASE, work); cleanup()
    f = lumapi.FDTD(hide=True)
    try:
        f.load(work); f.switchtolayout()
        f.setnamed("FDTD", "mesh accuracy", 2)
        f.setnamed("FDTD", "simulation time", 1200e-15)
        try:
            f.setresource("FDTD",1,"processes",1); f.setresource("FDTD",1,"threads",4)
        except Exception: pass
        # Single blend absorber matching the committed D2 pipeline geometry:
        # drop the discrete QD film and extend the perovskite object across the
        # full 850 nm active region that P_pvk_top / P_qd_edt bracket, so the
        # power-box extraction is valid. (Thickness is the 850 nm modelled
        # active layer, not Han's 490 nm; the above-gap fraction that this test
        # validates is thickness-independent.)
        if int(f.getnamednumber("QD-Pbs")):
            f.select("QD-Pbs"); f.delete()
        f.setnamed("Perovskite", "z min", 160e-9)
        if fill > 0:
            mk_blend(f, fill)
            f.setnamed("Perovskite", "material", "QD_han")
        print(f"running fill={fill:.4g} at acc 2, 850 nm active")
        f.run()
        lam, A, Iam, jopt = active_absorptance(f)
        agf = above_gap_fraction_weighted(lam, A, Iam, GAP_MAPBI3)
        np.savetxt(os.path.join(REPO,"results","fdtd",f"han_absorptance_f{fill:.0e}.csv"),
                   np.column_stack([lam, A]), delimiter=",",
                   header="lam_nm,A_active", comments="")
        print(f"  fill={fill:.4g}  Jopt_active={jopt:.3f}  above-gap frac={agf:.4f}")
        return jopt, agf
    except Exception as e:
        print(f"FAIL fill={fill}: {str(e)[:120]}")
        return None, None
    finally:
        try: f.close()
        except Exception: pass
        cleanup(); shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    fills = [float(a) for a in sys.argv[1:]] or FILLS
    if not os.path.isfile(OUT):
        with open(OUT,"w") as fh:
            fh.write("fill_volfrac,han_mg_per_mL,Jopt_active_mA_cm2,above_gap_fraction_800nm\n")
    def mg_label(fill):
        for k, v in ((1e-6,0.0),(3.3e-4,0.5),(6.6e-4,1.0)):
            if abs(fill-k) <= 0.01*k:
                return v
        return float("nan")           # 20 vol% anchor has no mg/mL analogue
    for fill in fills:
        jopt, agf = solve(fill)
        if jopt is not None:
            with open(OUT,"a") as fh:
                fh.write("%.6g,%.3g,%.3f,%.4f\n" % (fill, mg_label(fill), jopt, agf))
    print("done ->", OUT)
