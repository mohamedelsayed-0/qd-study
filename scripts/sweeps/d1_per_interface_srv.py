"""#3 Per-interface SRV: which interface actually controls each design?

D1 has three interfaces (SnO2/pvk, pvk/QD internal, QD/EDT). The common SRV
sweep folded them; this script sweeps each individually at S in {1e2, 1e4,
1e6} cm/s with e=h. Uses the new SRV_qedt_n/p, SRV_pvqd_n/p, SRV_snpv_n/p
knobs added in the D1 builder.
"""
import os, json, subprocess, time

REPO = r"C:\Users\alkin\OneDrive\Documents\qd-research\Quantom-dot-solar"
OUT = os.path.join(REPO, "results", "charge")
PY = r"C:\Program Files\Lumerical\v261\python\python.exe"
BUILD = os.path.join(REPO, "scripts", "charge", "build_planar.py")

INTERFACES = ["qedt", "pvqd", "snpv"]
S_VALS = ["1e2", "1e4", "1e6"]
rows = []

def run_one(env_extra, tag):
    manifest = os.path.join(OUT, "hedge_manifest.json")
    jv = os.path.join(OUT, "hedge_jv.csv")
    for f in (manifest, jv):
        try: os.remove(f)
        except OSError: pass
    env = dict(os.environ)
    env.update(dict(QD_CHI="3.85", QD_EG="1.03", EDT_CHI="3.9",
                    MAPBI3_TAU="1e-7", SNO2_CHI="4.0", QD_TAU="1e-5",
                    FF_OUT=OUT, FF_REPO=REPO))
    env.update(env_extra)
    cmd = [PY, BUILD, "--variant", "hedge",
           "--source", "models/charge/hedge_1p03eV.ldev",
           "--generation", "data/generation/charge_hedge.mat",
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
        return (m["Jsc_mA_cm2"], m["Voc_V"], m["FF"], m["PCE_percent_at_100mW_cm2"])
    except Exception:
        return None
    finally:
        try: os.remove(os.path.join(OUT, "%s.ldev" % tag))
        except OSError: pass

# baseline (all SRV zero)
r = run_one({}, "d1_srv_base")
if r:
    rows.append(("baseline", 0.0, "e+h") + r)
    print("baseline | Jsc=%.2f Voc=%.3f FF=%.3f PCE=%.2f" % r)

for iface in INTERFACES:
    for s in S_VALS:
        # symmetric e=h
        tag = "d1_srv_%s_%s" % (iface, s.replace("+", "p"))
        env = {"SRV_%s_n" % iface: s, "SRV_%s_p" % iface: s}
        r = run_one(env, tag)
        if r:
            rows.append((iface, float(s), "e+h") + r)
            print("%s S=%s e+h | Jsc=%.2f Voc=%.3f FF=%.3f PCE=%.2f" %
                  (iface, s, r[0], r[1], r[2], r[3]))
        # electron-only high, hole zero
        tag = "d1_srv_%s_e_%s" % (iface, s.replace("+", "p"))
        env = {"SRV_%s_n" % iface: s, "SRV_%s_p" % iface: "1e2"}
        r = run_one(env, tag)
        if r:
            rows.append((iface, float(s), "e_only") + r)
            print("%s S=%s e_only | Jsc=%.2f Voc=%.3f FF=%.3f PCE=%.2f" %
                  (iface, s, r[0], r[1], r[2], r[3]))

with open(os.path.join(OUT, "d1_interface_map.csv"), "w") as fh:
    fh.write("interface,S_cm_s,mode,Jsc,Voc,FF,PCE\n")
    for r in rows:
        fh.write("%s,%.0e,%s,%.2f,%.4f,%.3f,%.2f\n" % r)
print("done", len(rows))
