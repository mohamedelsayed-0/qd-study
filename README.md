# QD / perovskite solar cell

3D FDTD optics + 2D CHARGE drift-diffusion. Stack: FTO / SnO2 / MAPbI3 / PbS-QD / p-PbS-EDT / MoOx / Au.
D1 = planar QD film. D2 = 5 nm PbS QDs embedded in MAPbI3 (Maxwell-Garnett).

## Results

Two interpretations: **intermediate-band** (raw generation, every absorbed photon collects) and **strict single-gap** (spectrally filtered generation).

### Intermediate-band (upper bound)

| Device | Jsc | Voc | FF | PCE | Rs 0.5 / 1 / 2 Ω·cm² |
|---|---:|---:|---:|---:|---|
| D1 planar | 39.07 | 0.736 | 0.716 | 20.59% | — |
| D2 planar | 42.23 | 0.715 | 0.671 | 20.27% | — |
| **D1 optimized** | 42.7 | 0.741 | 0.759 | **24.00%** | 23.20 / 22.43 / 20.94 |
| **D2 optimized** | 50.2 | 0.719 | 0.671 | **24.21%** | 23.27 / 22.34 / 20.57 |
| **D2 high-voltage** | 36.9 | 0.915 | 0.737 | **24.91%** | 24.40 / 23.90 / 22.93 |

Rs crossover at 0.76 Ω·cm². Unpolarized: D1 23.89%, D2 24.09%. D2 doping sweep peaks higher: PCE 27.50% at 1e17 cm⁻³ (see `d2_doping_map_dense.csv`, fig13).

### Strict single-gap

| Device | Jsc | Voc | FF | PCE | Rs 0.5 / 1 / 2 Ω·cm² |
|---|---:|---:|---:|---:|---|
| D1 planar | 31.03 | 0.724 | 0.721 | 16.20% | 15.71 / 15.22 / 14.70 |
| D2 planar | 24.91 | 0.687 | 0.690 | 11.81% | 11.51 / 11.20 / 10.85 |
| **D1 optimized** | 33.86 | 0.729 | 0.762 | **18.80%** | 18.31 / 17.79 / 16.84 |
| **D2 optimized** | 29.69 | 0.691 | 0.694 | 14.23% | 13.79 / 13.35 / 12.84 |
| D2 high-voltage | 21.81 | 0.909 | 0.745 | 14.78% | 14.53 / 14.31 / 14.10 |

D1 wins at every Rs by 4–5 pt. Unpolarized: D1 18.77, D2 14.15.

## Reproduce

```powershell
python scripts/charge/presets.py --list
python scripts/charge/presets.py --preset d1_max --run
python scripts/charge/presets.py --preset d2_max --run
python scripts/spectral_filter.py
python scripts/charge/presets.py --preset d1_max_filtered --run
python scripts/build_grating_fdtd.py d1_max
python scripts/make_paper_figures.py
```
