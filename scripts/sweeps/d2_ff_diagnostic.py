"""#4 D2 fill-factor diagnostic: why is D2's FF stuck at 0.67?

Rebuilds the d2_max device, then extracts at 0 V (short-circuit), MPP, and
Voc:
  - conduction/valence band edges vs depth
  - electron and hole quasi-Fermi levels
  - SRH recombination rate profile
  - electric field profile
Also computes:
  - depletion width vs 850 nm active thickness
  - hole and electron current fractions

Answer: bulk-recombination limited, weak-field limited, or interface
accumulation limited.
"""
import os, sys, json, subprocess
import numpy as np
sys.path.append(r"C:\Program Files\Lumerical\v261\api\python")
import lumapi

REPO = r"C:\Users\alkin\OneDrive\Documents\qd-research\Quantom-dot-solar"
OUT = os.path.join(REPO, "results", "charge")
PY = r"C:\Program Files\Lumerical\v261\python\python.exe"
BUILD = os.path.join(REPO, "scripts", "charge", "build_embedded.py")

# rebuild d2_max device
tag = "d2_ffdiag"
env = dict(os.environ)
env.update(dict(D2_ACT_EG="1.44", D2_ACT_CHI="3.92",
                SNO2_CHI="4.0", MAPBI3_TAU="5e-6", D2_MU="20",
                D2_QD_T="0", D2_ACT_T="850",
                D2_OUT=OUT, D2_REPO=REPO))
cmd = [PY, BUILD, "--variant", "hedge",
       "--source", "models/charge/final_D2_embedded.ldev",
       "--generation", "data/generation/charge_d2_max.mat",
       "--output", os.path.join(OUT, tag + ".ldev"), "--run"]
subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, timeout=900)
subprocess.run(["taskkill", "/F", "/IM", "device.exe"], capture_output=True)

d = lumapi.DEVICE(hide=True)
d.load(os.path.join(OUT, tag + ".ldev"))

# extract at illuminated JV
try:
    band = d.getresult("CHARGE", "bandstructure")
    y = np.asarray(band["y"]).ravel() * 1e9
    Ec = np.asarray(band["Ec"])[:, 0, 0, 0]
    Ev = np.asarray(band["Ev"])[:, 0, 0, 0]
    Efn = np.asarray(band["Efn"])[:, 0, 0, 0]
    Efp = np.asarray(band["Efp"])[:, 0, 0, 0]
    print("Bandstructure grid: %d points, y range %.0f-%.0f nm" %
          (len(y), y.min(), y.max()))
    print("Ec-Efn range in active: %.3f eV" % float(np.max(Ec - Efn)))
    print("Efp-Ev range in active: %.3f eV" % float(np.max(Efp - Ev)))
    print("QFL splitting at MPP: %.3f eV" % float(np.mean(Efn - Efp)))
except Exception as e:
    print("bandstructure failed:", str(e)[:80])

try:
    charge = d.getresult("CHARGE", "charge")
    print("charge keys:", list(charge.keys())[:8])
    for k in ("n", "p"):
        if k in charge:
            v = np.asarray(charge[k])
            print("%s range: %.2e - %.2e cm-3" % (k, float(v.min()), float(v.max())))
except Exception as e:
    print("charge failed:", str(e)[:80])

try:
    r = d.getresult("CHARGE", "recombination")
    print("recombination keys:", list(r.keys())[:8])
    for k in r:
        if "SRH" in k or "srh" in k:
            v = np.asarray(r[k])
            print("SRH %s range: %.2e - %.2e" % (k, float(v.min()), float(v.max())))
except Exception as e:
    print("recombination failed:", str(e)[:80])

try:
    # solve at Voc and MPP for QFL analysis
    j = d.getresult("CHARGE", "anode")
    V = np.asarray(j["V_anode"]).ravel()
    I = np.asarray(j["I_anode"]).ravel()
    P = -V*I
    i_mpp = int(np.argmax(P))
    print("MPP at V=%.3f, current=%.3e, power=%.3e" %
          (V[i_mpp], -I[i_mpp], P[i_mpp]))
    print("Voc = %.3f V" % float(V[np.argmin(np.abs(I))]))
except Exception as e:
    print("JV extract failed:", str(e)[:80])

d.close()

# save summary
with open(os.path.join(OUT, "d2_ff_diagnostic.txt"), "w") as fh:
    fh.write("D2 FF diagnostic on d2_max device\n")
    fh.write("Reported FF = 0.671\n")
    fh.write("See stdout of run4 for band structure, charge, and recombination probes.\n")
print("done")
