"""Reproducible CHARGE device presets.

Each preset is the exact builder + generation + material parameter set behind a
headline result in the study. Examples:

    python scripts/charge/presets.py --list
    python scripts/charge/presets.py --preset final_d1 --run
    python scripts/charge/presets.py --preset final_d2 --run

Sweep points reuse a base preset with extra overrides (repeatable --set):

    python scripts/charge/presets.py --preset final_d1 --set SRV=1e3 --run      # interface recomb
    python scripts/charge/presets.py --preset final_d1 --set QD_EG=0.80 --run   # QD bandgap
    python scripts/charge/presets.py --preset final_d2 --set MAPBI3_TAU=1e-6 --run

Parameter sources are documented in logs/methods.md and scripts/charge/PRESETS.md.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CHARGE = REPO / "scripts" / "charge"


def lumerical_python() -> str:
    """Prefer the Lumerical python (ships lumapi); fall back to current."""
    found = sorted(Path("C:/Program Files/Lumerical").glob("v*/python/python.exe"), reverse=True)
    return str(found[0]) if found else sys.executable


PRESETS = {
    "final_d1": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_hedge.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "1e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5",
        },
        "expect": "Jsc ~39.1, Voc ~0.736, FF ~0.72, PCE ~20.6%",
    },
    "final_d2": {
        "builder": "build_embedded.py",
        "source": "models/charge/final_D2_embedded.ldev",
        "generation": "data/generation/charge_embedded.mat",
        "env": {
            "D2_ACT_EG": "1.44", "D2_ACT_CHI": "3.92",
            "SNO2_CHI": "4.0", "MAPBI3_TAU": "5e-6", "D2_QD_T": "0", "D2_ACT_T": "850",
        },
        "expect": "Jsc ~42.2, Voc ~0.715, FF ~0.67, PCE ~20.3%",
    },
    "d1_max": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_d1_max.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "5e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5", "QD_MU": "0.1",
        },
        "expect": "Jsc ~42.7, Voc ~0.741, FF ~0.76, PCE ~24.0%",
    },
    "d2_max": {
        "builder": "build_embedded.py",
        "source": "models/charge/final_D2_embedded.ldev",
        "generation": "data/generation/charge_d2_max.mat",
        "env": {
            "D2_ACT_EG": "1.44", "D2_ACT_CHI": "3.92",
            "SNO2_CHI": "4.0", "MAPBI3_TAU": "5e-6", "D2_MU": "20",
            "D2_QD_T": "0", "D2_ACT_T": "850",
        },
        "expect": "Jsc ~50.2, Voc ~0.719, FF ~0.67, PCE ~24.2%",
    },
    "d2_high_voltage": {
        "builder": "build_embedded.py",
        "source": "models/charge/final_D2_embedded.ldev",
        "generation": "data/generation/charge_d2_max.mat",
        "env": {
            "D2_ACT_EG": "1.44", "D2_ACT_CHI": "3.92",
            "SNO2_CHI": "4.0", "MAPBI3_TAU": "5e-6", "D2_MU": "20",
            "D2_ACT_DOP": "1e22", "D2_QD_T": "0", "D2_ACT_T": "850",
        },
        "expect": "Jsc ~36.9, Voc ~0.915, FF ~0.74, PCE ~24.9%",
    },
    # ---- strict single-gap variants (spectrally filtered generation) ----
    "final_d1_filtered": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_hedge_filtered.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "1e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5",
        },
        "expect": "single-gap D1 planar, Jsc ~31, PCE lower than 20.6%",
    },
    "final_d2_filtered": {
        "builder": "build_embedded.py",
        "source": "models/charge/final_D2_embedded.ldev",
        "generation": "data/generation/charge_embedded_filtered.mat",
        "env": {
            "D2_ACT_EG": "1.44", "D2_ACT_CHI": "3.92",
            "SNO2_CHI": "4.0", "MAPBI3_TAU": "5e-6", "D2_QD_T": "0", "D2_ACT_T": "850",
        },
        "expect": "single-gap D2 planar, Jsc ~25, PCE lower than 20.3%",
    },
    "d1_max_filtered": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_d1_max_filtered.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "5e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5", "QD_MU": "0.1",
        },
        "expect": "single-gap D1 optimized, Jsc ~34, PCE ~19%",
    },
    "d2_max_filtered": {
        "builder": "build_embedded.py",
        "source": "models/charge/final_D2_embedded.ldev",
        "generation": "data/generation/charge_d2_max_filtered.mat",
        "env": {
            "D2_ACT_EG": "1.44", "D2_ACT_CHI": "3.92",
            "SNO2_CHI": "4.0", "MAPBI3_TAU": "5e-6", "D2_MU": "20",
            "D2_QD_T": "0", "D2_ACT_T": "850",
        },
        "expect": "single-gap D2 optimized, Jsc ~30, PCE ~14-15%",
    },
    "d2_high_voltage_filtered": {
        "builder": "build_embedded.py",
        "source": "models/charge/final_D2_embedded.ldev",
        "generation": "data/generation/charge_d2_max_filtered.mat",
        "env": {
            "D2_ACT_EG": "1.44", "D2_ACT_CHI": "3.92",
            "SNO2_CHI": "4.0", "MAPBI3_TAU": "5e-6", "D2_MU": "20",
            "D2_ACT_DOP": "1e22", "D2_QD_T": "0", "D2_ACT_T": "850",
        },
        "expect": "single-gap doped D2, Jsc ~22, PCE ~14-15%",
    },
    # ---- matched-QD-loading D1 controls (review fairness check) ----
    # Same conservative transport as final_d1; only the QD film thickness and
    # its generation profile change, so the QD *volume* matches D2 at 20 / 25
    # vol.% while the 850 nm active thickness is held fixed.
    "d1_load170": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_d1_load170.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "1e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5",
            "QD_T": "170", "MAPBI3_T": "680",
        },
        "expect": "matched-loading D1, 170 nm QD (=D2 20 vol.%), raw regime",
    },
    "d1_load170_filtered": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_d1_load170_filtered.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "1e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5",
            "QD_T": "170", "MAPBI3_T": "680",
        },
        "expect": "matched-loading D1, 170 nm QD, strict single-gap",
    },
    "d1_load212p5": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_d1_load212p5.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "1e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5",
            "QD_T": "212.5", "MAPBI3_T": "637.5",
        },
        "expect": "matched-loading D1, 212.5 nm QD (=D2 25 vol.%), raw regime",
    },
    "d1_load212p5_filtered": {
        "builder": "build_planar.py",
        "source": "models/charge/hedge_1p03eV.ldev",
        "generation": "data/generation/charge_d1_load212p5_filtered.mat",
        "env": {
            "QD_CHI": "3.85", "QD_EG": "1.03", "EDT_CHI": "3.9",
            "MAPBI3_TAU": "1e-7", "SNO2_CHI": "4.0", "QD_TAU": "1e-5",
            "QD_T": "212.5", "MAPBI3_T": "637.5",
        },
        "expect": "matched-loading D1, 212.5 nm QD, strict single-gap",
    },
    # Conformal grating: infrastructure only. Setting CONFORMAL_TOOTH_H_NM
    # >0 (e.g. via --set on any D1 preset) installs a rectangular Au tooth
    # with stepped EDT/QD shells wrapping it. Jsc is currently under-reported
    # by ~12x because DEVICE's imported-generation dataset does not broadcast
    # across the widened sim x-domain — DEVICE support ticket needed. See
    # logs/manuscript_methods.md §Limitations for detail.
}


def run_preset(name: str, extra_env: dict, do_run: bool) -> int:
    p = PRESETS[name]
    out = REPO / "results" / "charge" / f"{name}.ldev"
    env = dict(os.environ)
    env.update(p["env"])
    env.update(extra_env)
    env.setdefault("FF_OUT", str(REPO / "results" / "charge"))
    env.setdefault("D2_OUT", str(REPO / "results" / "charge"))
    env.setdefault("FF_REPO", str(REPO))
    env.setdefault("D2_REPO", str(REPO))
    cmd = [
        lumerical_python(), str(CHARGE / p["builder"]),
        "--variant", "hedge", "--source", p["source"],
        "--generation", p["generation"], "--output", str(out),
    ]
    if do_run:
        cmd.append("--run")
    print(f"preset: {name}   expect: {p['expect']}")
    if extra_env:
        print("override:", extra_env)
    result = subprocess.run(cmd, cwd=str(REPO), env=env)
    subprocess.run(["taskkill", "/F", "/IM", "device.exe"], capture_output=True)
    manifest = REPO / "results" / "charge" / "hedge_manifest.json"
    if do_run and manifest.is_file():
        m = json.loads(manifest.read_text())["metrics"]
        print("RESULT  Jsc=%.2f  Voc=%.3f  FF=%.3f  PCE=%.2f%%" % (
            m["Jsc_mA_cm2"], m["Voc_V"], m["FF"], m["PCE_percent_at_100mW_cm2"]))
    return result.returncode


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="list presets and exit")
    ap.add_argument("--preset", choices=list(PRESETS))
    ap.add_argument("--set", action="append", default=[], metavar="ENV=VAL",
                    help="extra parameter override (repeatable), e.g. --set SRV=1e3")
    ap.add_argument("--run", action="store_true", help="solve CHARGE (omit to only build)")
    args = ap.parse_args()

    if args.list or not args.preset:
        for name, spec in PRESETS.items():
            print("%-10s -> %-18s | %s" % (name, spec["builder"], spec["expect"]))
        return

    try:
        extra = dict(kv.split("=", 1) for kv in args.set)
    except ValueError:
        ap.error("--set expects ENV=VAL pairs")
    sys.exit(run_preset(args.preset, extra, args.run))


if __name__ == "__main__":
    main()
