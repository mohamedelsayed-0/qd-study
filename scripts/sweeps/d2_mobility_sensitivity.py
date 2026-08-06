"""#1 D2 mobility sensitivity: is the 24.2% headline carried by mu=20 cm2/Vs?

Sweeps blend mobility from realistic pure-perovskite (~50) down to the raw
QD film value (0.02) at the optimized d2_max grating optics, holding all
else fixed. Reports Jsc/Voc/FF/PCE and PCE at Rs=1.
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

MUS = ["0.02", "0.1", "1", "5", "10", "20", "50"]
rows = []
for mu in MUS:
    tag = "d2_mu%s" % mu.replace(".", "p")
    manifest = os.path.join(OUT, "hedge_manifest.json")
    jv = os.path.join(OUT, "hedge_jv.csv")
    for f in (manifest, jv):
        try: os.remove(f)
        except OSError: pass
    env = dict(os.environ)
    env.update(dict(D2_ACT_EG="1.44", D2_ACT_CHI="3.92",
                    SNO2_CHI="4.0", MAPBI3_TAU="5e-6", D2_MU=mu,
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
            print("timeout", tag, attempt+1)
        subprocess.run(["taskkill", "/F", "/IM", "device.exe"], capture_output=True)
        if os.path.isfile(manifest): break
        time.sleep(15)
    try:
        m = json.load(open(manifest))["metrics"]
        prs = pce_rs(jv, 1.0) if os.path.isfile(jv) else float("nan")
        row = (mu, m["Jsc_mA_cm2"], m["Voc_V"], m["FF"],
               m["PCE_percent_at_100mW_cm2"], prs)
        rows.append(row)
        print("mu=%s | Jsc=%.2f Voc=%.3f FF=%.3f PCE=%.2f Rs1=%.2f" % row)
    except Exception as e:
        print("FAIL mu=%s: %s" % (mu, str(e)[:60]))
    for p in (os.path.join(OUT, "%s.ldev" % tag),):
        try: os.remove(p)
        except OSError: pass

with open(os.path.join(OUT, "d2_mobility_sensitivity.csv"), "w") as fh:
    fh.write("mu_cm2Vs,Jsc,Voc,FF,PCE,PCE_Rs1\n")
    for r in rows:
        fh.write("%s,%.2f,%.4f,%.3f,%.2f,%.2f\n" % r)
print("done", len(rows))
