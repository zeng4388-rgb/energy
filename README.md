# FRB Isotropic Equivalent Energy (E_iso) Function

Bayesian framework to measure the isotropic equivalent energy function of Fast Radio Bursts (FRBs), adapted from [Luo et al. 2018, MNRAS, 481, 2320](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.2320L/abstract) and [Luo et al. 2020, MNRAS, 494, 665](https://ui.adsabs.harvard.edu/abs/2020MNRAS..494..665L/abstract).

## Modifications from Original Code

1. **Python 3** - Converted from Python 2.7 to Python 3
2. **E_iso calculation** - Uses isotropic equivalent energy instead of luminosity
3. **CHIME/FRB data** - Uses CHIME Catalog 1, Catalog 2, Repeater and Baseband data
4. **Simplified structure** - All code and data in one folder with relative paths

## Requirements

- Python 3.7+
- NumPy >= 1.14
- SciPy >= 1.0
- PyMultiNest (see [GitHub](https://github.com/JohannesBuchner/PyMultiNest))
- Matplotlib
- MPI (for parallel MultiNest runs)

## Data Files

| File | Description |
|------|-------------|
| `chime_cat1.json` | CHIME/FRB Catalog 1 (600 FRBs, 2018-2019) |
| `chime_cat2.json` | CHIME/FRB Catalog 2 (5045 FRBs) |
| `chime_repeaters.json` | CHIME/FRB 2023 Repeater Catalog (151 bursts) |
| `chime_pre_repeaters.json` | Previously Published Repeaters (345 bursts) |

## Code Files

| File | Description |
|------|-------------|
| `frb_util.py` | Core utility classes (Cosmology, Telescope, AstroDistribution, EventRate, Loadfiles) |
| `nest_samp.py` | Bayesian inference for CHIME/FRB real data |
| `nest_simu.py` | Bayesian inference for simulated data |
| `pltpost.py` | Posterior distribution plotting |
| `simufrb.py` | FRB simulation using E_iso Schechter function |

## Quick Start

### 1. Test with CHIME Catalog 1 (single survey, no MPI)

```bash
python3 nest_samp.py -fc chime_cat1.json -o test_cat1 -g ALG_YMW16 -upper 1
```

### 2. Run full CHIME Catalog 1 analysis

```bash
bash run_chime_cat1.sh
```

### 3. Run CHIME Catalog 2 analysis (one-off and repeaters)

```bash
bash run_chime_cat2.sh
```

### 4. Run repeating FRB analysis

```bash
bash run_chime_repeat.sh
```

### 5. Simulate and recover E_iso function

```bash
bash run_simu.sh
```

### 6. Plot posterior distributions

```bash
bash draw_results.sh
```

## Parameters

The E_iso Schechter function has 6 parameters:

| Parameter | Description | Prior Range (flat) |
|-----------|-------------|-------------------|
| `phis` | Event rate density (Gpc^-3 yr^-1) | [1e-3, 3e4] |
| `alpha` | Power-law index | [-3, 1.1] |
| `logEiso_max` | log10 characteristic energy (erg) | [48, 53] |
| `logEiso_0` | log10 minimum energy (erg) | [43, 48] |
| `mu_w` | Mean log pulse width | [-1, 2] |
| `sigma_w` | Std log pulse width | [0.1, 1.0] |

## E_iso Calculation

The isotropic equivalent energy is calculated as:

```
E_iso = fluence / (1+z) * 4*pi*D_L^2 * dnu
```

where:
- `fluence` is the observed fluence in Jy ms
- `z` is the redshift (estimated from DM)
- `D_L` is the luminosity distance
- `dnu` is the intrinsic bandwidth (default 1000 MHz)

## CHIME Survey Parameters

| Parameter | Catalog 1 | Catalog 2 |
|-----------|-----------|-----------|
| Frequency | 400-800 MHz | 400-800 MHz |
| Bandwidth | 400 MHz | 400 MHz |
| FoV | ~200 deg^2 | ~200 deg^2 |
| Gain | ~1.4 K/Jy | ~1.4 K/Jy |
| Tsys | ~25 K | ~25 K |
| Npol | 2 | 2 |
| SNR threshold | ~10 | ~10 |
| Obs time | ~213 days | ~1 year |

## Citation

If you use this code, please cite:
- Luo et al. 2018, MNRAS, 481, 2320
- Luo et al. 2020, MNRAS, 494, 665
- CHIME/FRB Collaboration 2021, ApJS, 257, 59 (Catalog 1)
- CHIME/FRB Collaboration 2023 (Repeater Catalog)
