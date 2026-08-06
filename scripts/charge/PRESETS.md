# CHARGE presets

Reproducible parameter sets behind every headline electrical result. Run through
`presets.py`, which selects the builder, generation profile, and material
parameters, then solves CHARGE.

```powershell
python scripts/charge/presets.py --list
python scripts/charge/presets.py --preset final_d1 --run
python scripts/charge/presets.py --preset final_d2 --run
```

Output: `results/charge/<preset>.ldev` + `results/charge/hedge_manifest.json`
(Jsc, Voc, FF, PCE). Uses the Lumerical python automatically (ships `lumapi`).

## Headline devices

| Preset | Builder | Generation | Jsc | Voc | FF | PCE |
|---|---|---|---:|---:|---:|---:|
| `final_d1` | build_planar.py | charge_hedge.mat | 39.07 | 0.736 | 0.716 | 20.59% |
| `final_d2` | build_embedded.py | charge_embedded.mat | 42.23 | 0.715 | 0.671 | 20.27% |
| `d1_max` | build_planar.py | charge_d1_max.mat | 42.65 | 0.741 | 0.759 | 24.00% |
| `d2_max` | build_embedded.py | charge_d2_max.mat | 50.17 | 0.719 | 0.671 | 24.21% |
| `d2_high_voltage` | build_embedded.py | charge_d2_max.mat | 36.92 | 0.915 | 0.737 | 24.91% |

All five verified to reproduce the values above from a clean checkout. The
`*_max` generation profiles come from the textured-grating FDTD runs
(D1 grating 600/200 nm + front glass texture; D2 grating 600/350 nm, 25% fill).
`d2_high_voltage` uses the same optics as `d2_max` and raises blend doping
to 1e16 cm^-3 (env D2_ACT_DOP=1e22, in m^-3), trading current for voltage.

**Reported Jsc caveat.** Every preset assumes each absorbed photon produces
one collectable carrier — this is CHARGE's default with no wavelength cutoff
at the electrical gap. FDTD absorptance extends below every absorber's gap,
so the reported Jsc is an intermediate-band upper bound. Under strict
single-gap physics the corrected Jsc values drop by ~20% for D1 and ~40% for
D2 (see `scripts/energy_budget.py` and `results/fdtd/energy_budget.csv`).

## Parameter sets

`final_d1` — QD_CHI 3.85 eV, QD_EG 1.03 eV, EDT_CHI 3.9 eV, MAPBI3_TAU 1e-7 s,
SNO2_CHI 4.0 eV, QD_TAU 1e-5 s.

`final_d2` — D2_ACT_EG 1.44 eV, D2_ACT_CHI 3.92 eV, SNO2_CHI 4.0 eV,
MAPBI3_TAU 5e-6 s, blend collapsed (D2_QD_T 0, D2_ACT_T 850 nm).

Physical justification for each value is in `logs/methods.md` and `charge_params/`.

## Sweeps

Sweep points reuse a base preset with repeatable `--set ENV=VAL` overrides, so
every figure is reproducible without a separate script:

```powershell
# interface recombination (fig4)
python scripts/charge/presets.py --preset final_d1 --set SRV=1e3 --run
python scripts/charge/presets.py --preset final_d2 --set SRV=1e3 --run
# QD bandgap (fig5)
python scripts/charge/presets.py --preset final_d1 --set QD_EG=0.80 --run
# D2 blend lifetime (fig10 tornado)
python scripts/charge/presets.py --preset final_d2 --set MAPBI3_TAU=1e-6 --run
```

Override keys: `QD_CHI QD_EG EDT_CHI SNO2_CHI QD_TAU MAPBI3_TAU SRV QD_DOP
EDT_DOP MAPBI3_DOP QD_T MAPBI3_T FTO_T MOO3_T TEMP_K MAX_EDGE` (D1);
`D2_ACT_EG D2_ACT_CHI SNO2_CHI MAPBI3_TAU SRV D2_ACT_DOP D2_ACT_T TEMP_K` (D2).
Contacts are ohmic; `ANODE_WF`/`CATHODE_WF` exist but are no-ops by construction.
Note: thickness-changing knobs shift the mesh against the saved source partition;
electrical solves are only validated at the default stack height (see log).
