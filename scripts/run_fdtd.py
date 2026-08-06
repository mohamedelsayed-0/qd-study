"""Run an FDTD script without creating persistent working copies in the repo."""

from __future__ import annotations

import argparse
import importlib
import re
import shutil
import sys
import tempfile
from pathlib import Path


def load_lumapi():
    roots = sorted(Path("C:/Program Files/Lumerical").glob("v*/api/python"), reverse=True)
    if not roots:
        raise RuntimeError("No Lumerical Python API installation was found")
    sys.path.insert(0, str(roots[0]))
    return importlib.import_module("lumapi")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("script", help="LSF filename without the .lsf extension")
    parser.add_argument(
        "--base",
        default="models/fdtd/final_planar_double_junction.fsp",
        help="base FSP copied into the Windows temporary directory",
    )
    parser.add_argument("--visible", action="store_true", help="show the FDTD window")
    parser.add_argument(
        "--pattern",
        choices=("square", "hex", "honeycomb"),
        help="run one explicit_5nm_patterns organization in a clean process",
    )
    parser.add_argument(
        "--mesh-nm",
        type=float,
        help="tagged lateral mesh refinement for one explicit pattern",
    )
    args = parser.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9_]+", args.script):
        raise ValueError("Script name may contain only letters, numbers, and underscores")
    if args.pattern and args.script != "explicit_5nm_patterns":
        raise ValueError("--pattern is only valid with explicit_5nm_patterns")
    if args.mesh_nm is not None:
        if not args.pattern:
            raise ValueError("--mesh-nm requires --pattern")
        if not 0.4 <= args.mesh_nm <= 2.0:
            raise ValueError("--mesh-nm must be between 0.4 and 2.0 nm")

    script_dir = Path(__file__).resolve().parent
    repo = script_dir.parent
    source = (repo / args.base).resolve()
    if repo not in source.parents or not source.is_file():
        raise FileNotFoundError(f"Base model is not inside the repository: {source}")

    search_dirs = (
        script_dir,
        script_dir / "qd_patterns",
        script_dir / "common",
    )
    matches = [
        directory / f"{args.script}.lsf"
        for directory in search_dirs
        if (directory / f"{args.script}.lsf").is_file()
    ]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one script named {args.script}.lsf; found {len(matches)}"
        )

    with tempfile.TemporaryDirectory(prefix="qd_research_fdtd_") as temporary:
        work_dir = Path(temporary)
        working_model = work_dir / "working_model.fsp"
        shutil.copy2(source, working_model)
        print(f"Temporary model: {working_model}")

        lumapi = load_lumapi()
        session = lumapi.FDTD(hide=not args.visible)
        try:
            session.load(str(working_model))
            path_commands = " ".join(
                f'addpath("{str(path).replace(chr(92), "/")}");'
                for path in search_dirs
            )
            output_dir = str(work_dir).replace("\\", "/")
            setup = ""
            run_tag = args.pattern
            if args.pattern:
                setup = f'requested_pattern="{args.pattern}"; '
            if args.mesh_nm is not None:
                mesh_tag = f"{args.mesh_nm:g}".replace(".", "p")
                run_tag = f"{args.pattern}_mesh{mesh_tag}nm"
                setup += (
                    f"pattern_mesh_xy_nm={args.mesh_nm:g}; "
                    f'pattern_result_tag="{run_tag}"; '
                )
            session.eval(
                f'{path_commands} cd("{output_dir}"); {setup}{args.script};'
            )
        finally:
            session.close()

        if args.script in {"final_device", "final_embedded_5nm_effective"}:
            output_sets = {
                "final_device": {
                    "final_planar_double_junction.fsp": (
                        repo
                        / "models"
                        / "fdtd"
                        / "final_planar_double_junction.fsp"
                    ),
                    "G_planar_double_junction_raw.mat": (
                        repo / "data" / "generation" / "planar_fdtd_raw.mat"
                    ),
                    "final_device_out.mat": (
                        repo / "results" / "fdtd" / "planar_results.mat"
                    ),
                },
                "final_embedded_5nm_effective": {
                    "final_embedded_5nm_effective.fsp": (
                        repo
                        / "models"
                        / "fdtd"
                        / "final_embedded_5nm_effective.fsp"
                    ),
                    "G_embedded_5nm_raw.mat": (
                        repo
                        / "data"
                        / "generation"
                        / "embedded_5nm_fdtd_raw.mat"
                    ),
                    "embedded_5nm_results.mat": (
                        repo / "results" / "fdtd" / "embedded_5nm_results.mat"
                    ),
                },
            }
            outputs = output_sets[args.script]
            for source_name, destination in outputs.items():
                generated = work_dir / source_name
                if not generated.is_file():
                    raise FileNotFoundError(f"Expected FDTD output was not created: {generated}")
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(generated, destination)
                print(f"Updated: {destination}")
        else:
            if args.script == "explicit_5nm_patterns":
                if args.pattern:
                    expected = (
                        f"qd_5nm_{run_tag}_out.mat",
                        f"qd_5nm_{run_tag}_results.csv",
                        f"qd_5nm_{run_tag}_log.txt",
                    )
                else:
                    expected = (
                        "qd_5nm_pattern_out.mat",
                        "qd_5nm_pattern_results.csv",
                        "qd_5nm_pattern_log.txt",
                    )
                missing = [name for name in expected if not (work_dir / name).is_file()]
                if missing:
                    raise FileNotFoundError(
                        "Explicit pattern run was incomplete; missing "
                        + ", ".join(missing)
                    )
            output_dir = repo / "results" / "fdtd" / args.script
            output_dir.mkdir(parents=True, exist_ok=True)
            for generated in work_dir.iterdir():
                if (
                    generated.suffix.lower()
                    in {".mat", ".txt", ".csv", ".npz", ".log"}
                    and not generated.name.endswith("_p0.log")
                ):
                    destination = output_dir / generated.name
                    shutil.copy2(generated, destination)
                    print(f"Saved: {destination}")

    print(f"Completed {args.script}; temporary working files removed")


if __name__ == "__main__":
    main()
