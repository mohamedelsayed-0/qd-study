"""#2 D2 doping map: is the high-voltage point robust across doping range?

Sweeps D2 blend p-doping from 1e14 to 5e18 cm-3 at fixed d2_max grating
optics. Reports Jsc/Voc/FF/PCE and PCE at Rs=0,0.5,1,2 Ohm cm2 for each.
The paper claim is that lower-current, higher-voltage devices dominate at
realistic contact resistance.
"""
import os, json, subprocess, time
import numpy as np

REPO = r"C:\Users\alkin\OneDrive\Documents\qd-research\Quantom-dot-solar"
OUT = os.path.join(REPO, "results", "charge")
PY = r"C:\Program Files\Lumerical\v261\python\python.exe"
BUILD = os.path.join(REPO, "scripts", "charge", "build_embedded.py")

def pce_rs(jv, rs):
    d = np.loadtxt(jv, delimiter=",", skiprows=1)
    V, J = d[:,0], d[:,1]
    return float(np.max((V - J*1e-3*rs)*J)/100.0*100.0)

# cm-3 to m-3: mult by 1e6
DOPINGS_CM3 = [1e14, 1e15, 5e15, 1e16, 5e16, 1e17, 5e17, 1e18]
rows = []
for dop_cm3 in DOPINGS_CM3:
    dop_m3 = dop_cm3 * 1e6
    tag = "d2_dop_%.0e" % dop_cm3
    manifest = os.path.join(OUT, "hedge_manifest.json")
    jv = os.path.join(OUT, "hedge_jv.csv")
    for f in (manifest, jv):
        try: os.remove(f)
        except OSError: pass
    env = dict(os.environ)
    env.update(dict(D2_ACT_EG="1.44", D2_ACT_CHI="3.92",
                    SNO2_CHI="4.0", MAPBI3_TAU="5e-6", D2_MU="20",
                    D2_ACT_DOP="%.4e" % dop_m3,
                    D2_QD_T="0", D2_ACT_T="850",
                    D2_OUT=OUT, D2_REPO=REPO))
    cmd = [PY, BUILD, "--variant", "hedge",
           "--source", "models/charge/final_D2_embedded.ldev",
           "--generation", "data/generation/charge_d2_max.mat",
           "--output", os.path.join(OUT, "%s.ldev" % tag), "--run"]
    for attempt in range(3):
        try:
            subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, timeout=900)
        except subprocess.TimeoutExpired:
            print("timeout", tag)
        subprocess.run(["taskkill", "/F", "/IM", "device.exe"], capture_output=True)
        if os.path.isfile(manifest): break
        time.sleep(15)
    try:
        m = json.load(open(manifest))["metrics"]
        rs_vals = [pce_rs(jv, r) if os.path.isfile(jv) else float("nan")
                   for r in (0.0, 0.5, 1.0, 2.0)]
        row = (dop_cm3, m["Jsc_mA_cm2"], m["Voc_V"], m["FF"],
               m["PCE_percent_at_100mW_cm2"]) + tuple(rs_vals)
        rows.append(row)
        print("dop=%.0e | Jsc=%.2f Voc=%.3f FF=%.3f PCE=%.2f Rs2=%.2f" %
              (dop_cm3, m["Jsc_mA_cm2"], m["Voc_V"], m["FF"],
               m["PCE_percent_at_100mW_cm2"], rs_vals[3]))
    except Exception as e:
        print("FAIL dop=%.0e: %s" % (dop_cm3, str(e)[:60]))
    for p in (os.path.join(OUT, "%s.ldev" % tag),):
        try: os.remove(p)
        except OSError: pass

with open(os.path.join(OUT, "d2_doping_map.csv"), "w") as fh:
    fh.write("dop_cm3,Jsc,Voc,FF,PCE_Rs0,PCE_Rs0p5,PCE_Rs1,PCE_Rs2\n")
    for r in rows:
        fh.write("%.0e,%.2f,%.4f,%.3f,%.2f,%.2f,%.2f,%.2f\n" % r)
print("done", len(rows))
