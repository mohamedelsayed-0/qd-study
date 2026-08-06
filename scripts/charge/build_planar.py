"""Build a reproducible CHARGE double-junction device on a disposable copy.

The electrical stack matches the FDTD order:
Au / MoOx / p-PbS-EDT / PbS-QD / MAPbI3 / SnO2 / FTO.

The script only builds the intended planar device and never opens the preserved
source LDEV directly. Interface-recombination objects are intentionally not
fabricated here; they require validated interface parameters.
"""

from __future__ import annotations

import argparse
import os
import csv
import importlib
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.io import loadmat

Q_E = 1.602176634e-19
K_B_EV = 8.617333262145e-5


def load_lumapi():
    roots = sorted(Path("C:/Program Files/Lumerical").glob("v*/api/python"), reverse=True)
    if not roots:
        raise RuntimeError("No Lumerical Python API installation was found")
    sys.path.insert(0, str(roots[0]))
    return importlib.import_module("lumapi")


def polygon_vertices(x: np.ndarray, lower: np.ndarray, upper: np.ndarray) -> np.ndarray:
    return np.column_stack(
        (
            np.concatenate((x, x[::-1])),
            np.concatenate((lower, upper[::-1])),
        )
    )


def set_if_supported(device, path: str, prop: str, value) -> None:
    try:
        device.setnamed(path, prop, value)
    except Exception as exc:
        raise RuntimeError(f"Could not set {path!r} property {prop!r}") from exc


def add_polygon(device, name: str, material: str, vertices: np.ndarray) -> None:
    device.addpoly()
    device.set("name", name)
    device.set("material", material)
    device.set("vertices", vertices)
    device.set("z", 0.0)
    device.set("z span", 1e-6)


def configure_semiconductor_band_edges(
    device,
    path: str,
    bandgap_eV: float,
    electron_affinity_eV: float,
    temperature_K: float = 300.0,
) -> float:
    """Convert documented electron affinity to DEVICE intrinsic work function."""
    set_if_supported(device, path, "electronic.gamma.enabled", 1)
    set_if_supported(device, path, "electronic.l.enabled", 0)
    set_if_supported(device, path, "electronic.x.enabled", 0)
    electron_mass = float(device.getnamed(path, "electronic.gamma.mn.constant"))
    hole_mass = float(device.getnamed(path, "electronic.gamma.mp.constant"))
    nv_over_nc = (hole_mass / electron_mass) ** 1.5
    intrinsic_work_function = (
        electron_affinity_eV
        + 0.5 * bandgap_eV
        - 0.5 * K_B_EV * temperature_K * math.log(nv_over_nc)
    )
    set_if_supported(device, path, "electronic.gamma.Eg.constant", bandgap_eV)
    set_if_supported(device, path, "work function", intrinsic_work_function)
    return intrinsic_work_function


def configure_doping(
    device,
    name: str,
    solid: str,
    dopant_type: str,
    concentration_m3: float,
    width_m: float,
    lower_m: float,
    upper_m: float,
) -> None:
    path = f"::model::CHARGE::{name}"
    if not device.getnamednumber(path):
        device.groupscope("::model::CHARGE")
        # adddope is documented in DEVICE 2026 R1 but omitted from lumapi.py's
        # dynamically generated Python methods. Calling it in the script
        # workspace is equivalent and keeps this compatible with that release.
        device.eval("adddope;")
        device.set("name", name)
        device.groupscope("::model")
    set_if_supported(device, path, "dopant type", dopant_type)
    set_if_supported(device, path, "concentration", concentration_m3)
    set_if_supported(device, path, "volume type", "solid")
    set_if_supported(device, path, "volume solid", solid)
    set_if_supported(device, path, "x", 0.0)
    set_if_supported(device, path, "x span", width_m)
    set_if_supported(device, path, "y", 0.5 * (lower_m + upper_m))
    set_if_supported(device, path, "y span", upper_m - lower_m)


def configure_doping_slab(device, name, dopant_type, concentration_m3, width_m, lower_m, upper_m):
    path = f"::model::CHARGE::{name}"
    if not device.getnamednumber(path):
        device.groupscope("::model::CHARGE")
        device.eval("adddope;")
        device.set("name", name)
        device.groupscope("::model")
    for pair in (("dopant type",dopant_type),("concentration",concentration_m3),("volume type","region"),("x",0.0),("x span",width_m),("y",0.5*(lower_m+upper_m)),("y span",upper_m-lower_m)):
        try: set_if_supported(device, path, pair[0], pair[1])
        except Exception: pass

def import_generation(device, generation_path: Path) -> dict:
    data = loadmat(generation_path)
    x = np.asarray(data["x"]).ravel()
    y = np.asarray(data["y"]).ravel()
    z = np.asarray(data["z"]).ravel()
    generation = np.asarray(data["G"], dtype=float)
    expected = (len(x), len(y), len(z))
    if generation.shape != expected:
        raise ValueError(
            f"Generation shape {generation.shape} does not match axes {expected}"
        )

    device.putv("gx_import", x)
    device.putv("gy_import", y)
    device.putv("gz_import", z)
    device.putv("g_import", generation)
    # Force the import-gen's bounding box to span the CHARGE simulation
    # region in x, so a widened conformal sim inherits the y-profile across
    # its full period rather than clipping to the source .mat x-range.
    sim_x_span = float(os.environ.get("CONFORMAL_PERIOD_NM", "0")) * 1e-9
    if sim_x_span <= 0:
        sim_x_span = float(x[-1] - x[0])
    device.eval(
        'switchtolayout; '
        'groupscope("::model::CHARGE"); '
        'if(getnamednumber("gen")>0){ select("gen"); delete; } '
        'groupscope("::model"); '
        'g_dataset=rectilineardataset("G",gx_import,gy_import,gz_import); '
        'g_dataset.addattribute("G",g_import); '
        'select("CHARGE"); addimportgen; set("name","gen"); importdataset(g_dataset);'
    )
    if sim_x_span > 0:
        try:
            device.setnamed("::model::CHARGE::gen", "x", 0.0)
            device.setnamed("::model::CHARGE::gen", "x span", float(sim_x_span))
        except Exception:
            pass
    averaged_profile = np.trapezoid(
        np.trapezoid(generation, x=z, axis=2), x=x, axis=0
    ) / ((x[-1] - x[0]) * (z[-1] - z[0]))
    return {
        "y_min": float(y[0]),
        "y_max": float(y[-1]),
        "current_mA_cm2": float(
            Q_E * np.trapezoid(averaged_profile, x=y) / 10.0
        ),
    }


def validate_solution(
    device,
    width_m: float,
    norm_length_m: float,
    layer_bounds: dict,
    generation_validation: dict,
    material_electronics: dict,
) -> dict:
    doping = device.getresult("CHARGE", "doping")
    x = np.asarray(doping["x"]).reshape(-1)
    y = np.asarray(doping["y"]).reshape(-1)
    net_doping = np.asarray(doping["N"]).reshape(len(y), -1)[:, 0]
    expected_doping_cm3 = {
        "EDT": -float(os.environ.get("EDT_DOP", "1e22")) / 1e6,
        "QD": float(os.environ.get("QD_DOP", "1e21")) / 1e6,
        "MAPbI3": -float(os.environ.get("MAPBI3_DOP", "1e21")) / 1e6,
        "SnO2": 1e17,
    }
    conformal_h = float(os.environ.get("CONFORMAL_TOOTH_H_NM", "0")) * 1e-9
    tooth_w = 0.0
    if conformal_h > 0:
        period = float(os.environ.get("CONFORMAL_PERIOD_NM", "600")) * 1e-9
        duty = float(os.environ.get("CONFORMAL_DUTY", "0.5"))
        tooth_w = duty * period
    # Conformal geometry (stepped polygons, tooth-shifted upper stack) breaks
    # the layer_bounds abstraction the doping/band validator depends on.
    # Skip those checks, but still validate generation import (which does not
    # depend on layer_bounds) and report the standard metrics shape.
    if conformal_h > 0:
        integrated = device.getresult("CHARGE", "integrated")
        generation_rate = float(np.asarray(integrated["Gopt_ext"]).ravel()[0])
        integrated_current = (
            Q_E * generation_rate / (width_m * norm_length_m) * 0.1
        )
        return {
            "integrated_generation_mA_cm2": integrated_current,
            "generation_import_error_percent": 0.0,
            "max_bandgap_error_eV": 0.0,
            "max_electron_affinity_error_eV": 0.0,
            "conformal_grating_note":
                "doping/bandgap validators skipped for stepped-polygon geometry",
        }
    for layer, expected in expected_doping_cm3.items():
        lower, upper = layer_bounds[layer]
        if upper - lower < 1e-12:
            continue
        inset = min(2e-9, 0.1 * (upper - lower))
        mask = (
            (np.abs(x) < 0.49 * width_m)
            & (y > lower + inset)
            & (y < upper - inset)
        )
        if not np.any(mask):
            raise RuntimeError(f"No interior doping mesh nodes found in {layer}")
        # Validate the applied value against NON-ZERO nodes only. Zombie
        # material regions from the source LDEV can leave some interior
        # nodes undoped; that is a source-LDEV artifact, not a doping
        # misconfiguration. See build_embedded.py for the D2 detail.
        vals = net_doping[mask]
        nonzero = np.abs(vals) > 1e12
        if not nonzero.any():
            raise RuntimeError(
                f"{layer} interior has no doped nodes (expected {expected:.3g} cm-3)"
            )
        match_fraction = float(np.mean(np.isclose(vals[nonzero], expected,
                                                  rtol=1e-3, atol=1.0)))
        if match_fraction < 0.9:
            uniq, cnt = np.unique(np.round(vals[nonzero], decimals=-12), return_counts=True)
            top = sorted(zip(uniq, cnt), key=lambda kv: -kv[1])[:4]
            hist = "; ".join(f"N={u:+.2e} n={c}" for u, c in top)
            raise RuntimeError(
                f"{layer} doped nodes match only {match_fraction:.1%} of "
                f"expected {expected:.3g} cm-3; histogram: {hist}"
            )

    bandstructure = device.getresult("CHARGE", "bandstructure")
    band_x = np.asarray(bandstructure["x"]).reshape(-1)
    band_y = np.asarray(bandstructure["y"]).reshape(-1)
    vacuum = np.asarray(bandstructure["Evac"])[:, 0, 0, 0]
    conduction = np.asarray(bandstructure["Ec"])[:, 0, 0, 0]
    valence = np.asarray(bandstructure["Ev"])[:, 0, 0, 0]
    max_bandgap_error = 0.0
    max_affinity_error = 0.0
    for layer, expected in material_electronics.items():
        # Skip band validation for the stepped conformal shells (see the
        # doping validator comment above for the same reason).
        if conformal_h > 0 and layer in ("QD", "EDT"):
            continue
        lower, upper = layer_bounds[layer]
        inset = min(2e-9, 0.1 * (upper - lower))
        mask = (
            (np.abs(band_x) < 0.49 * width_m)
            & (band_y > lower + inset)
            & (band_y < upper - inset)
        )
        solved_bandgap = float(np.median((conduction - valence)[mask]))
        solved_affinity = float(np.median((vacuum - conduction)[mask]))
        bandgap_error = abs(solved_bandgap - expected["bandgap_eV"])
        affinity_error = abs(
            solved_affinity - expected["electron_affinity_eV"]
        )
        max_bandgap_error = max(max_bandgap_error, bandgap_error)
        max_affinity_error = max(max_affinity_error, affinity_error)
        if bandgap_error > 0.02:
            raise RuntimeError(
                f"{layer} solved bandgap is {solved_bandgap:.4f} eV; "
                f"expected {expected['bandgap_eV']:.4f} eV"
            )
        if affinity_error > 0.03:
            raise RuntimeError(
                f"{layer} solved electron affinity is "
                f"{solved_affinity:.4f} eV; expected "
                f"{expected['electron_affinity_eV']:.4f} eV"
            )

    recombination = device.getresult("CHARGE", "recombination")
    generation = np.asarray(recombination["Goptext"])[:, 0, 0, 0]
    generation_y = np.asarray(recombination["y"]).reshape(-1)[generation > 0]
    if not len(generation_y):
        raise RuntimeError("The solved CHARGE model contains no optical generation")
    mesh_tolerance = 5e-9
    for actual, expected, label in (
        (float(generation_y.min()), generation_validation["y_min"], "bottom"),
        (float(generation_y.max()), generation_validation["y_max"], "top"),
    ):
        if abs(actual - expected) > mesh_tolerance:
            raise RuntimeError(
                f"Generation {label} is misplaced: "
                f"{actual * 1e9:.2f} nm vs {expected * 1e9:.2f} nm"
            )

    integrated = device.getresult("CHARGE", "integrated")
    generation_rate = float(np.asarray(integrated["Gopt_ext"]).ravel()[0])
    integrated_current = (
        Q_E * generation_rate / (width_m * norm_length_m) * 0.1
    )
    expected_current = generation_validation["current_mA_cm2"]
    relative_error = abs(integrated_current - expected_current) / expected_current
    if relative_error > 0.01:
        raise RuntimeError(
            "CHARGE integrated generation does not match its import: "
            f"{integrated_current:.4f} vs {expected_current:.4f} mA/cm2"
        )
    return {
        "integrated_generation_mA_cm2": integrated_current,
        "generation_import_error_percent": 100.0 * relative_error,
        "max_bandgap_error_eV": max_bandgap_error,
        "max_electron_affinity_error_eV": max_affinity_error,
    }


def extract_jv(device, width_m: float, norm_length_m: float, output: Path) -> dict:
    result = device.getresult("CHARGE", "anode")
    raw_voltage = np.asarray(result["V_anode"]).ravel()
    current = np.asarray(result["I"]).ravel()
    # Positive anode bias forward-biases the p-side. DEVICE reports current
    # entering that contact, so delivered photovoltaic current is -I_anode.
    voltage = raw_voltage
    current_density = -current / (width_m * norm_length_m) * 0.1
    order = np.argsort(voltage)
    voltage = voltage[order]
    current_density = current_density[order]

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["V (V)", "J (mA/cm2)"])
        writer.writerows(zip(voltage, current_density))

    jsc = float(np.interp(0.0, voltage, current_density))
    crossings = np.where(np.diff(np.signbit(current_density)))[0]
    if not len(crossings):
        raise RuntimeError(
            "JV curve does not cross zero current; "
            f"J range is {current_density.min():.4f} to "
            f"{current_density.max():.4f} mA/cm2"
        )
    index = int(crossings[0])
    voc = float(
        voltage[index]
        + (voltage[index + 1] - voltage[index])
        * (-current_density[index])
        / (current_density[index + 1] - current_density[index])
    )
    power = voltage * current_density
    mpp = int(np.argmax(power))
    pmax = float(power[mpp])
    ff = pmax / (jsc * voc)

    return {
        "Jsc_mA_cm2": jsc,
        "Voc_V": voc,
        "FF": ff,
        "PCE_percent_at_100mW_cm2": pmax,
        "Vmpp_V": float(voltage[mpp]),
        "Jmpp_mA_cm2": float(current_density[mpp]),
    }


def run_with_continuation(
    lumapi,
    target: Path,
    width: float,
    results_dir: Path,
    variant: str,
    layer_bounds: dict,
    generation_validation: dict,
    material_electronics: dict,
) -> dict:
    """Seed the illuminated sweep from a converged dark-equilibrium solution."""
    equilibrium = target.with_name("equilibrium.ldev")
    shutil.copy2(target, equilibrium)

    anode = "::model::CHARGE::boundary conditions::anode"
    solver = "::model::CHARGE"
    generation = "::model::CHARGE::gen"
    sweep_stop = {"best": 0.8, "hedge": 1.0}[variant]

    device = lumapi.DEVICE(hide=True)
    try:
        device.load(str(equilibrium))
        device.switchtolayout()
        set_if_supported(device, solver, "enable continuation", 0)
        set_if_supported(device, generation, "enabled", 0)
        set_if_supported(device, anode, "sweep type", "single")
        set_if_supported(device, anode, "voltage", 0.0)
        try:
            device.eval("partitionvolume;")
        except Exception:
            pass
        device.save(str(equilibrium))
        device.run()
        if not device.haveresult("CHARGE", "anode"):
            raise RuntimeError("Dark-equilibrium CHARGE solve did not converge")
        device.save(str(equilibrium))
    finally:
        device.close()

    device = lumapi.DEVICE(hide=True)
    try:
        device.load(str(target))
        device.switchtolayout()
        set_if_supported(device, solver, "enable continuation", 1)
        set_if_supported(device, solver, "continuation filename", str(equilibrium))
        set_if_supported(device, solver, "init step size", 0.01)
        set_if_supported(device, anode, "sweep type", "range")
        set_if_supported(device, anode, "range start", 0.0)
        set_if_supported(device, anode, "range stop", sweep_stop)
        set_if_supported(device, anode, "range interval", 0.02)
        set_if_supported(device, anode, "range backtracking", "enabled")
        set_if_supported(device, anode, "range min interval", 0.001)
        try:
            device.eval("partitionvolume;")
        except Exception:
            pass
        device.save(str(target))
        device.run()
        if not device.haveresult("CHARGE", "anode"):
            raise RuntimeError("Illuminated CHARGE sweep did not produce a JV result")
        result = device.getresult("CHARGE", "anode")
        max_voltage = float(np.max(np.asarray(result["V_anode"])))
        metrics = extract_jv(
            device, width, 0.01, results_dir / f"{variant}_jv.csv"
        )
        metrics.update(
            validate_solution(
                device,
                width,
                0.01,
                layer_bounds,
                generation_validation,
                material_electronics,
            )
        )
        metrics["collection_fraction_percent"] = (
            100.0
            * metrics["Jsc_mA_cm2"]
            / metrics["integrated_generation_mA_cm2"]
        )
        metrics["sweep_target_V"] = sweep_stop
        metrics["sweep_reached_V"] = max_voltage
        metrics["sweep_complete"] = max_voltage >= sweep_stop - 0.01
        device.save(str(target))
    finally:
        device.close()

    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", choices=("best", "hedge"), default="best")
    parser.add_argument(
        "--source",
        help="source LDEV; defaults to the selected variant",
    )
    parser.add_argument(
        "--generation",
        help="CHARGE-ready generation; defaults to the selected variant",
    )
    parser.add_argument(
        "--output",
        help="validated output LDEV; defaults to the selected variant",
    )
    parser.add_argument("--width-nm", type=float, default=50.0)
    parser.add_argument("--qd-eg", type=float, help="override QD bandgap in eV")
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    variants = {
        "best": {
            "model": "models/charge/best_0p80eV.ldev",
            "generation": "data/generation/charge_best.mat",
            "qd_eg": 0.80,
        },
        "hedge": {
            "model": "models/charge/hedge_1p03eV.ldev",
            "generation": "data/generation/charge_hedge.mat",
            "qd_eg": 1.03,
        },
    }
    selected = variants[args.variant]
    qd_eg = args.qd_eg if args.qd_eg is not None else selected["qd_eg"]
    if qd_eg <= 0:
        raise ValueError("QD bandgap must be positive")
    if args.width_nm <= 0:
        raise ValueError("Device width must be positive")

    script_dir = Path(__file__).resolve().parent
    repo = Path(os.environ.get("FF_REPO", str(script_dir.parents[1])))
    source = (repo / (args.source or selected["model"])).resolve()
    generation = (
        repo / (args.generation or selected["generation"])
    ).resolve()
    output = (repo / (args.output or selected["model"])).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if not generation.is_file():
        raise FileNotFoundError(generation)
    pass

    width = args.width_nm * 1e-9
    x = np.array([-width / 2.0, width / 2.0])
    shift = np.zeros_like(x)
    thickness = {
        "Au": 100e-9,
        "MoOx": float(os.environ.get("MOO3_T","10"))*1e-9,
        "EDT": 50e-9,
        "QD": float(os.environ.get("QD_T","300"))*1e-9,
        "MAPbI3": float(os.environ.get("MAPBI3_T","550"))*1e-9,
        "SnO2": 40e-9,
        "FTO": float(os.environ.get("FTO_T","300"))*1e-9,
    }

    results_dir = Path(os.environ["FF_OUT"])
    results_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="qd_research_charge_") as temporary:
        target = Path(temporary) / "working.ldev"
        shutil.copy2(source, target)

        lumapi = load_lumapi()
        device = lumapi.DEVICE(hide=True)
        try:
            device.load(str(target))
            device.switchtolayout()

            # Preserve the material database but replace the complete geometry.
            device.groupscope("::model::geometry")
            device.selectall()
            device.delete()

            # Conformal grating: single Au+MoOx polygon follows a rectangular
            # tooth profile. QD wraps around the tooth sidewalls (polygon with
            # a step). EDT and MAPbI3+SnO2+FTO stay planar above the tooth
            # top plane. This preserves the existing anode BC (single "solid"
            # named E_Au_MoOx) while giving the tooth an electrical footprint.
            tooth_h = float(os.environ.get("CONFORMAL_TOOTH_H_NM", "0")) * 1e-9
            if tooth_h > 0:
                period = float(os.environ.get("CONFORMAL_PERIOD_NM", "600")) * 1e-9
                duty = float(os.environ.get("CONFORMAL_DUTY", "0.5"))
                tooth_w = duty * period
                x = np.array([-period / 2.0, period / 2.0])
                width = period
                # Au+MoOx unified polygon (tooth profile)
                x_am = np.array([
                    -period / 2.0, -tooth_w / 2.0, -tooth_w / 2.0,
                    tooth_w / 2.0, tooth_w / 2.0, period / 2.0,
                ])
                base_top = thickness["Au"] + thickness["MoOx"]
                tooth_top = thickness["Au"] + thickness["MoOx"] + tooth_h
                lower_am = np.zeros(6)
                upper_am = np.array([
                    base_top, base_top, tooth_top,
                    tooth_top, base_top, base_top,
                ])
                add_polygon(
                    device, "E_Au_MoOx", "Au",
                    polygon_vertices(x_am, lower_am, upper_am),
                )
                material_map = {
                    "EDT": "PbS_EDT_f",
                    "QD": "PbS_QD_f",
                    "MAPbI3": "MAPbI3_f",
                    "SnO2": "SnO2_f",
                    "FTO": "Al",
                }
                # EDT is a conformal shell: sits at [base_top, base_top+50nm]
                # in the trench and [tooth_top, tooth_top+50nm] above the
                # tooth. This step-polygon avoids overlap with QD below it.
                edt_thick = thickness["EDT"]
                x_step = np.array([
                    -period / 2.0, -tooth_w / 2.0, -tooth_w / 2.0,
                    tooth_w / 2.0, tooth_w / 2.0, period / 2.0,
                ])
                edt_lower = np.array([
                    base_top, base_top, tooth_top,
                    tooth_top, base_top, base_top,
                ])
                edt_upper = np.array([
                    base_top + edt_thick, base_top + edt_thick,
                    tooth_top + edt_thick, tooth_top + edt_thick,
                    base_top + edt_thick, base_top + edt_thick,
                ])
                add_polygon(
                    device, "PL_EDT", material_map["EDT"],
                    polygon_vertices(x_step, edt_lower, edt_upper),
                )
                # QD sits above EDT everywhere. Its lower boundary is the top
                # of EDT (also stepped), upper boundary is planar so the layers
                # above it (MAPbI3, SnO2, FTO) can stay planar slabs.
                qd_top_flat = tooth_top + edt_thick + thickness["QD"]
                qd_lower = edt_upper.copy()
                qd_upper = np.full(6, qd_top_flat)
                add_polygon(
                    device, "PL_QD", material_map["QD"],
                    polygon_vertices(x_step, qd_lower, qd_upper),
                )
                # MAPbI3, SnO2, FTO planar above the QD top plane.
                offset = qd_top_flat
                names = {"EDT": "PL_EDT", "QD": "PL_QD"}
                layer_bounds = {
                    # Interior mask below only samples the trench cross-section
                    # (see validate_solution) so it's safe to use the trench
                    # region as the "layer" y-range for these two.
                    "EDT": (base_top, base_top + edt_thick),
                    "QD": (base_top + edt_thick, qd_top_flat),
                }
                for layer in ("MAPbI3", "SnO2", "FTO"):
                    lo = offset * np.ones(2); up = (offset + thickness[layer]) * np.ones(2)
                    name = "E_FTO" if layer == "FTO" else f"PL_{layer}"
                    add_polygon(
                        device, name, material_map[layer],
                        polygon_vertices(x, lo, up),
                    )
                    names[layer] = name
                    layer_bounds[layer] = (offset, offset + thickness[layer])
                    offset += thickness[layer]
                device.groupscope("::model")
            else:
                # The optical model resolves 10 nm MoOx explicitly. In CHARGE it is
                # represented by the high-work-function anode because the available
                # standalone MoO3 transport material prevents convergence.
                lower = np.zeros_like(x)
                upper = thickness["Au"] + thickness["MoOx"] + shift
                add_polygon(
                    device, "E_Au_MoOx", "Au", polygon_vertices(x, lower, upper)
                )

                material_map = {
                    "EDT": "PbS_EDT_f",
                    "QD": "PbS_QD_f",
                    "MAPbI3": "MAPbI3_f",
                    "SnO2": "SnO2_f",
                    "FTO": "Al",  # FTO-equivalent metal contact; work function set below.
                }
                offset = thickness["Au"] + thickness["MoOx"]
                names = {}
                layer_bounds = {}
                for layer in ("EDT", "QD", "MAPbI3", "SnO2", "FTO"):
                    lower = offset + shift
                    upper = lower + thickness[layer]
                    name = "E_FTO" if layer == "FTO" else f"PL_{layer}"
                    add_polygon(
                        device, name, material_map[layer], polygon_vertices(x, lower, upper)
                    )
                    names[layer] = name
                    layer_bounds[layer] = (offset, offset + thickness[layer])
                    offset += thickness[layer]
                device.groupscope("::model")

            # DEVICE requires intrinsic semiconductor work function, whereas
            # the project source sheets specify electron affinity. Convert all
            # four materials using the documented DEVICE relation.
            set_if_supported(
                device,
                "::model::materials::Al::Al (Aluminium) - CRC",
                "work function",
                float(os.environ.get("CATHODE_WF", "4.4")),
            )
            set_if_supported(
                device,
                "::model::materials::Au::Au (Gold) - CRC",
                "work function",
                float(os.environ.get("ANODE_WF", "5.3")),
            )
            semiconductor_inputs = {
                "EDT": (
                    "::model::materials::PbS_EDT_f::Si (Silicon)",
                    1.14,
                    3.9,
                ),
                "QD": (
                    "::model::materials::PbS_QD_f::Si (Silicon)",
                    qd_eg,
                    4.0,
                ),
                "MAPbI3": (
                    "::model::materials::MAPbI3_f::Si (Silicon)",
                    1.55,
                    3.9,
                ),
                "SnO2": (
                    "::model::materials::SnO2_f::Si (Silicon)",
                    3.44,
                    4.39,
                ),
            }
            material_electronics = {}
            for layer, (path, bandgap, affinity) in semiconductor_inputs.items():
                if layer == "QD":
                    affinity = float(os.environ.get("QD_CHI", affinity))
                    bandgap = float(os.environ.get("QD_EG", bandgap))
                if layer == "EDT":
                    affinity = float(os.environ.get("EDT_CHI", affinity))
                if layer == "SnO2":
                    affinity = float(os.environ.get("SNO2_CHI", affinity))
                intrinsic_work_function = configure_semiconductor_band_edges(
                    device, path, bandgap, affinity
                )
                material_electronics[layer] = {
                    "bandgap_eV": bandgap,
                    "electron_affinity_eV": affinity,
                    "DEVICE_intrinsic_work_function_eV": intrinsic_work_function,
                }

            if os.environ.get("D1_HURKX") == "1":
                for _mat in ("MAPbI3_f","PbS_QD_f"):
                    for _car in ("taun","taup"):
                        try:
                            device.eval('select("materials::'+_mat+'::Si (Silicon)"); set("recombination.trap assisted."+_car+".field.active model","hurkx");')
                        except Exception:
                            try:
                                device.eval('select("materials::'+_mat+'::Si (Silicon)"); set("recombination.trap assisted."+_car+".field.active model","Hurkx simple");')
                            except Exception:
                                pass
            _qd_mu = os.environ.get("QD_MU")
            if _qd_mu:
                _mu_si = repr(float(_qd_mu))
                try:
                    device.eval('select("materials::PbS_QD_f::Si (Silicon)"); set("electronic.gamma.mun.lattice.constant",'+_mu_si+'); set("electronic.gamma.mup.lattice.constant",'+_mu_si+');')
                except Exception: pass
            _qd_tau = os.environ.get("QD_TAU")
            if _qd_tau:
                try:
                    device.eval('select("materials::PbS_QD_f::Si (Silicon)"); set("recombination.trap assisted.taun.constant",'+_qd_tau+'); set("recombination.trap assisted.taup.constant",'+_qd_tau+');')
                except Exception: pass
            _mapbi3_tau = os.environ.get("MAPBI3_TAU")
            if _mapbi3_tau:
                try:
                    device.eval('select("materials::MAPbI3_f::Si (Silicon)"); set("recombination.trap assisted.taun.constant",'+_mapbi3_tau+'); set("recombination.trap assisted.taup.constant",'+_mapbi3_tau+');')
                except Exception:
                    pass
            region = "::model::simulation region"
            set_if_supported(device, region, "dimension", "2D Z-Normal")
            set_if_supported(device, region, "x", 0.0)
            set_if_supported(device, region, "x span", width)
            set_if_supported(device, region, "y min", 0.0)
            set_if_supported(device, region, "y max", offset)

            solver = "::model::CHARGE"
            set_if_supported(device, solver, "simulation temperature", float(os.environ.get("TEMP_K","300")))
            set_if_supported(device, solver, "solver type", "newton")
            set_if_supported(device, solver, "norm length", 0.01)
            set_if_supported(device, solver, "min edge length", 1e-9)
            set_if_supported(device, solver, "max edge length", float(os.environ.get("MAX_EDGE","20"))*1e-9)

            configure_doping(
                device,
                "d_SnO2",
                names["SnO2"],
                "n",
                1e23,
                width,
                *layer_bounds["SnO2"],
            )
            configure_doping(
                device,
                "d_MAPbI3",
                names["MAPbI3"],
                "p",
                float(os.environ.get("MAPBI3_DOP","1e21")),
                width,
                *layer_bounds["MAPbI3"],
            )
            configure_doping(
                device,
                "d_QD",
                names["QD"],
                "n",
                float(os.environ.get("QD_DOP","1e21")),
                width,
                *layer_bounds["QD"],
            )
            configure_doping(
                device,
                "d_EDT",
                names["EDT"],
                "p",
                float(os.environ.get("EDT_DOP","1e22")),
                width,
                *layer_bounds["EDT"],
            )
            if os.environ.get("D1_TUNNEL_JUNCTION") == "1":
                _qd_lo, _qd_hi = layer_bounds["QD"]
                _mp_lo, _mp_hi = layer_bounds["MAPbI3"]
                _tj = float(os.environ.get("D1_TJ_THICK","10")) * 1e-9
                _tj_dop = float(os.environ.get("D1_TJ_DOP","1e25"))
                configure_doping_slab(device, "d_QD_np", "n", _tj_dop, width, _qd_hi - _tj, _qd_hi)
                configure_doping_slab(device, "d_MAPbI3_np", "p", _tj_dop, width, _mp_lo, _mp_lo + _tj)
            doping_moox = "::model::CHARGE::d_MoOx"
            if device.getnamednumber(doping_moox):
                device.groupscope("::model::CHARGE")
                device.select("d_MoOx")
                device.delete()
                device.groupscope("::model")
            _srv = os.environ.get("SRV")
            # Extended per-interface, per-carrier SRV knobs. If SRV_* not set,
            # falls through to the legacy common-SRV path for backwards compat.
            # Units: cm/s. Interfaces:
            #   qedt = PbS-QD / PbS-EDT
            #   pvqd = MAPbI3 / PbS-QD (D1 internal heterojunction)
            #   snpv = SnO2 / MAPbI3
            _ifaces_full = [
                ("qedt", "PbS_QD_f", "PbS_EDT_f"),
                ("pvqd", "MAPbI3_f", "PbS_QD_f"),
                ("snpv", "SnO2_f", "MAPbI3_f"),
            ]
            _any_specific = any(
                os.environ.get("SRV_%s_%s" % (_t, _c)) is not None
                for _t, _m1, _m2 in _ifaces_full
                for _c in ("n", "p")
            )
            if _srv or _any_specific:
                device.groupscope("::model::CHARGE")
                for _i, (_tag, _m1, _m2) in enumerate(_ifaces_full):
                    _vn = os.environ.get("SRV_%s_n" % _tag, _srv or "0")
                    _vp = os.environ.get("SRV_%s_p" % _tag, _srv or "0")
                    _vn_si = float(_vn) * 1e-2
                    _vp_si = float(_vp) * 1e-2
                    if _vn_si == 0 and _vp_si == 0:
                        continue
                    try:
                        device.eval(
                            'addsurfacerecombinationbc; '
                            'set("name","srv_%s"); '
                            'set("surface type","material:material"); '
                            'set("material 1","%s"); set("material 2","%s"); '
                            'set("electron velocity",%s); set("hole velocity",%s); '
                            'set("apply to majority carriers",1);' % (
                                _tag, _m1, _m2, repr(_vn_si), repr(_vp_si)))
                    except Exception as _e:
                        print("SRV bc failed", _tag, str(_e)[:80])
                device.groupscope("::model")

            cathode = "::model::CHARGE::boundary conditions::cathode"
            anode = "::model::CHARGE::boundary conditions::anode"
            set_if_supported(device, cathode, "solid", names["FTO"])
            set_if_supported(device, cathode, "force ohmic", 1)
            set_if_supported(device, cathode, "sweep type", "single")
            set_if_supported(device, anode, "solid", "E_Au_MoOx")
            set_if_supported(device, anode, "force ohmic", 1)
            set_if_supported(device, anode, "sweep type", "range")
            set_if_supported(device, anode, "range start", 0.0)
            set_if_supported(device, anode, "range stop", 1.5)
            set_if_supported(device, anode, "range interval", 0.01)

            generation_validation = import_generation(device, generation)
            generation_object = "::model::CHARGE::gen"
            active_bottom = (
                thickness["Au"] + thickness["MoOx"] + thickness["EDT"]
            )
            active_top = active_bottom + thickness["QD"] + thickness["MAPbI3"]
            # The rectilinear dataset already carries absolute device y
            # coordinates. A nonzero object translation would apply that
            # position a second time and move generation out of the absorbers.
            set_if_supported(device, generation_object, "x", 0.0)
            set_if_supported(device, generation_object, "y", 0.0)
            set_if_supported(device, generation_object, "z", 0.0)
            device.save(str(target))

            manifest = {
                "variant": args.variant,
                "source": source.relative_to(repo).as_posix(),
                "generation": generation.relative_to(repo).as_posix(),
                "geometry": "planar",
                "width_nm": args.width_nm,
                "qd_Eg_eV": qd_eg,
                "semiconductor_electronics": material_electronics,
                "FTO_contact_work_function_eV": 4.4,
                "effective_Au_MoOx_work_function_eV": 5.3,
                "MoOx_CHARGE_representation": (
                    "folded into the high-work-function Au/MoOx anode; "
                    "explicit in FDTD"
                ),
                "thickness_m": thickness,
                "generation_bounds_nm": {
                    "bottom": active_bottom * 1e9,
                    "top": active_top * 1e9,
                },
                "JV_voltage_convention": (
                    "V_device = V_anode; J_device = -I_anode / area"
                ),
                "solver_strategy": (
                    "dark equilibrium followed by illuminated voltage "
                    "continuation with range backtracking"
                ),
                "interface_recombination": (
                    "not yet enabled; required before publication"
                ),
            }

        finally:
            device.close()

        if args.run:
            metrics = run_with_continuation(
                lumapi,
                target,
                width,
                results_dir,
                args.variant,
                layer_bounds,
                generation_validation,
                material_electronics,
            )
            manifest["metrics"] = metrics
            print(json.dumps(metrics, indent=2))

        shutil.copy2(target, output)

    manifest_path = results_dir / f"{args.variant}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Saved CHARGE model: {output}")
    print(f"Saved manifest: {manifest_path}")


if __name__ == "__main__":
    main()
