"""Build and run the optimized (textured-grating) FDTD models.

Regenerates the optical generation profiles behind the d1_max / d2_max presets
from the committed baseline FSP models:

    python scripts/build_grating_fdtd.py d1_max   # grating 600/200 + front glass texture, 1.03 eV film
    python scripts/build_grating_fdtd.py d2_max   # grating 600/350, 25%-fill Maxwell-Garnett blend

Output: results/fdtd/G_<config>.mat (raw 3D generation) and the CHARGE-ready
laterally averaged profile data/generation/charge_<config>.mat. Runs on a
scratch copy of the FSP; the committed models are never modified. Mesh
accuracy: d2_max at 3 (converged); d1_max at 2 (accuracy 3 exceeds memory on
the reference machine — see results/fdtd/mesh_convergence.csv).
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import scipy.io as sio

c = 299792458.0
Q = 1.602176634e-19
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(REPO, "models", "fdtd", "final_planar_double_junction.fsp")
QDTXT = os.path.join(REPO, "data", "optical", "qd", "by_ligand", "PbS_1200_NH4I.txt")
QD5 = os.path.join(REPO, "data", "optical", "qd", "by_size", "PbS-QD_5nm.txt")

OBJS = ["Glass", "FTO", "SnO2", "Perovskite", "QD-Pbs", "p-PbS-EDT", "MOo3", "Au",
        "FDTD", "PlaneWave", "R", "T", "P_pvk_top", "P_pvk_qd", "P_qd_edt",
        "solar_generation"]

CONFIGS = {
    "d1_max": dict(acc=2, period=600e-9, tooth=200e-9, fill=None, front_tex=True,
                   override=12e-9),
    "d2_max": dict(acc=3, period=600e-9, tooth=350e-9, fill=0.25, front_tex=False,
                   override=16e-9),
    "d1_flat_ref": dict(acc=2, period=600e-9, tooth=0.0, fill=None, front_tex=False,
                        override=12e-9),
}


def load_lumapi():
    import importlib
    roots = sorted(__import__("pathlib").Path("C:/Program Files/Lumerical").glob("v*/api/python"),
                   reverse=True)
    if not roots:
        raise RuntimeError("No Lumerical Python API found")
    sys.path.insert(0, str(roots[0]))
    return importlib.import_module("lumapi")


def mk_blend(f, fill):
    hostmat = f.getnamed("Perovskite", "material")
    Nw = 200
    lam = np.linspace(300e-9, 1700e-9, Nw)
    fg = c / lam
    nh = np.array([complex(f.getfdtdindex(hostmat, fi, fi, fi)) for fi in fg])
    rows = []
    for ln in open(QD5):
        p = ln.split()
        if len(p) == 3:
            try:
                rows.append([float(v) for v in p])
            except ValueError:
                pass
    dq = np.array(rows)
    nq = np.interp(lam, dq[:, 0] * 1e-9, dq[:, 1]) + 1j * np.interp(lam, dq[:, 0] * 1e-9, dq[:, 2])
    eps_h = nh ** 2
    eps_i = nq ** 2
    beta = fill * (eps_i - eps_h) / (eps_i + 2 * eps_h)
    eps_eff = eps_h * (1 + 2 * beta) / (1 - beta)
    sd = np.zeros((Nw, 2), dtype=complex)
    sd[:, 0] = fg
    sd[:, 1] = eps_eff
    f.eval("if(materialexists('QD_blendX')==0){ m=addmaterial('Sampled 3D data'); "
           "setmaterial(m,'name','QD_blendX'); }")
    f.putv("sdb", sd)
    f.eval("setmaterial('QD_blendX','sampled data',sdb);")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", choices=list(CONFIGS))
    args = ap.parse_args()
    cfg = CONFIGS[args.config]

    lumapi = load_lumapi()
    tmp = tempfile.mkdtemp(prefix="qd_grating_")
    work = os.path.join(tmp, "work.fsp")
    shutil.copy2(BASE, work)
    f = lumapi.FDTD(hide=True)
    try:
        f.load(work)
        f.switchtolayout()
        f.setnamed("FDTD", "mesh accuracy", cfg["acc"])
        f.setnamed("FDTD", "simulation time", 1200e-15)
        try:
            f.setresource("FDTD", 1, "processes", 1)
            f.setresource("FDTD", 1, "threads", 4)
        except Exception:
            pass
        if cfg["fill"] is not None:
            mk_blend(f, cfg["fill"])
            f.select("QD-Pbs")
            f.delete()
            f.setnamed("Perovskite", "z min", 160e-9)
            f.setnamed("Perovskite", "material", "QD_blendX")
        else:
            d = np.loadtxt(QDTXT)
            lamq = d[:, 0] * 1e-9
            sd = np.zeros((len(lamq), 2), dtype=complex)
            sd[:, 0] = c / lamq
            sd[:, 1] = (d[:, 1] + 1j * d[:, 2]) ** 2
            f.eval("if(materialexists('QD_1p03')==0){ m=addmaterial('Sampled 3D data'); "
                   "setmaterial(m,'name','QD_1p03'); }")
            f.putv("sdq", sd)
            f.eval("setmaterial('QD_1p03','sampled data',sdq);")
            f.setnamed("QD-Pbs", "material", "QD_1p03")
        P = cfg["period"]
        ys = f.getnamed("FDTD", "y span")
        for o in OBJS:
            if int(f.getnamednumber(o)):
                f.setnamed(o, "x", 0)
                f.setnamed(o, "x span", P)
        if cfg["tooth"] > 0:
            f.addrect(); f.set("name", "Au_tooth")
            f.set("x", 0); f.set("x span", P / 2)
            f.set("y", 0); f.set("y span", ys)
            f.set("z min", 100e-9); f.set("z max", 100e-9 + cfg["tooth"])
            f.set("material", "Au (Gold) - Palik Copy 1")
            f.set("override mesh order from material database", 1); f.set("mesh order", 1)
            f.addmesh(); f.set("name", "tooth_mesh")
            f.set("x", 0); f.set("x span", P / 2 + 60e-9)
            f.set("y", 0); f.set("y span", ys)
            f.set("z min", 80e-9); f.set("z max", 100e-9 + cfg["tooth"] + 40e-9)
            f.set("override x mesh", 1); f.set("dx", cfg["override"])
            f.set("override z mesh", 1); f.set("dz", cfg["override"])
        if cfg["front_tex"]:
            f.setnamed("FDTD", "z max", 2100e-9)
            f.setnamed("PlaneWave", "z", 1900e-9)
            try:
                f.setnamed("R", "z", 2000e-9)
            except Exception:
                pass
            gm = f.getnamed("Glass", "material")
            f.addrect(); f.set("name", "gtex")
            f.set("x", 0); f.set("x span", P / 2)
            f.set("y", 0); f.set("y span", ys)
            f.set("z min", 1650e-9); f.set("z max", 1850e-9)
            f.set("material", gm)
        print("running", args.config)
        f.run()
        f.runanalysis("solar_generation")
        G = f.getresult("solar_generation", "G")
        x = np.asarray(G["x"]).ravel(); y = np.asarray(G["y"]).ravel(); z = np.asarray(G["z"]).ravel()
        raw = np.asarray(G["G"], dtype=float)
        g = raw if raw.shape == (len(x), len(y), len(z)) else raw.transpose(2, 1, 0)
        prof = np.trapezoid(np.trapezoid(g, x=x, axis=0), x=y, axis=0) / ((x[-1] - x[0]) * (y[-1] - y[0]))
        jopt = Q * np.trapezoid(prof, x=z) / 10.0
        print("Jopt = %.2f mA/cm2" % jopt)
        raw_out = os.path.join(REPO, "results", "fdtd", "G_%s.mat" % args.config)
        sio.savemat(raw_out, {"G_x": x, "G_y": y, "G_z": z, "G_3d": g})
        Gc = np.broadcast_to(prof[None, :, None], (2, len(z), 2)).copy()
        charge_out = os.path.join(REPO, "data", "generation", "charge_%s.mat" % args.config)
        sio.savemat(charge_out, {"x": np.array([[-25e-9], [25e-9]]), "y": z[:, None],
                                 "z": np.array([[-0.5e-6], [0.5e-6]]), "G": Gc},
                    do_compression=True)
        print("wrote", raw_out)
        print("wrote", charge_out)
    finally:
        try:
            f.close()
        except Exception:
            pass
        for exe in ("fdtd-engine-msmpi.exe", "fdtd-solutions.exe"):
            subprocess.run(["taskkill", "/F", "/IM", exe], capture_output=True)
        time.sleep(2)
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
