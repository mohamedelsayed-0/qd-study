"""Consolidate explicit-QD pattern runs and assess ranking resolution."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parents[1]
RESULTS_DIR = REPO / "results" / "fdtd" / "explicit_5nm_patterns"
NAME = re.compile(
    r"qd_5nm_(square|hex|honeycomb)"
    r"(?:_mesh(?P<mesh>[0-9]+p[0-9]+)nm)?_results\.csv"
)


def main() -> None:
    records = []
    for path in sorted(RESULTS_DIR.glob("qd_5nm_*_results.csv")):
        match = NAME.fullmatch(path.name)
        if not match:
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != 1:
            raise ValueError(f"Expected one result row in {path}")
        row = rows[0]
        mesh_text = row.get("mesh_xy_nm") or "0.75"
        record = {
            "pattern": match.group(1),
            "mesh_xy_nm": float(mesh_text),
            "J0_mA_cm2": float(row["J0_mA_cm2"]),
            "J90_mA_cm2": float(row["J90_mA_cm2"]),
            "Junpolarized_mA_cm2": float(row["Junpolarized_mA_cm2"]),
            "achieved_fill": float(row["achieved_fill"]),
            "nearest_spacing_nm": float(row["nearest_spacing_nm"]),
            "max_energy_bound_error": float(row["max_energy_bound_error"]),
        }
        records.append(record)

    expected = {
        (pattern, mesh)
        for pattern in ("square", "hex", "honeycomb")
        for mesh in (0.75, 0.6)
    }
    actual = {(row["pattern"], row["mesh_xy_nm"]) for row in records}
    if actual != expected:
        raise RuntimeError(f"Pattern/mesh coverage mismatch: {actual}")

    records.sort(key=lambda row: (-row["mesh_xy_nm"], row["pattern"]))
    summary_csv = RESULTS_DIR / "summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    finest_mesh = min(row["mesh_xy_nm"] for row in records)
    refined = [row for row in records if row["mesh_xy_nm"] == finest_mesh]
    ranking = sorted(
        refined, key=lambda row: row["Junpolarized_mA_cm2"], reverse=True
    )
    pattern_spread = (
        ranking[0]["Junpolarized_mA_cm2"]
        - ranking[-1]["Junpolarized_mA_cm2"]
    )
    mesh_shifts = {}
    for pattern in ("square", "hex", "honeycomb"):
        values = [
            row["Junpolarized_mA_cm2"]
            for row in records
            if row["pattern"] == pattern
        ]
        mesh_shifts[pattern] = max(values) - min(values)
    max_mesh_shift = max(mesh_shifts.values())
    max_polarization_spread = max(
        abs(row["J0_mA_cm2"] - row["J90_mA_cm2"]) for row in records
    )
    uncertainty = max(max_mesh_shift, max_polarization_spread)
    ranking_resolved = pattern_spread > uncertainty

    summary = {
        "dot_diameter_nm": 5.0,
        "target_fill_fraction": 0.20,
        "achieved_fill_fraction": records[0]["achieved_fill"],
        "finest_mesh_nm": finest_mesh,
        "finest_mesh_ranking": [
            {
                "pattern": row["pattern"],
                "Junpolarized_mA_cm2": row["Junpolarized_mA_cm2"],
            }
            for row in ranking
        ],
        "finest_mesh_pattern_spread_mA_cm2": pattern_spread,
        "mesh_shift_by_pattern_mA_cm2": mesh_shifts,
        "max_mesh_shift_mA_cm2": max_mesh_shift,
        "max_polarization_spread_mA_cm2": max_polarization_spread,
        "ranking_resolved": ranking_resolved,
        "conclusion": (
            "Pattern ranking is not resolved at the current numerical "
            "uncertainty; treat square, hex, and honeycomb as optically "
            "equivalent until a finer isotropic mesh converges."
            if not ranking_resolved
            else "Pattern ranking exceeds the current numerical uncertainty."
        ),
    }
    summary_json = RESULTS_DIR / "summary.json"
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(summary["conclusion"])
    print(f"Saved: {summary_csv}")
    print(f"Saved: {summary_json}")


if __name__ == "__main__":
    main()
