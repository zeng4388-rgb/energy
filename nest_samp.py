#!/usr/bin/env python3
"""
nest_samp.py - Bayesian inference for CHIME/FRB sample (E_iso version)
Uses PyMultiNest to fit the isotropic equivalent energy function
Supports CHIME Catalog 1, Catalog 2, Repeater data

Usage:
    python nest_samp.py -fc chime_cat1.json -o cat1_ALG_YMW16 -g ALG_YMW16 -upper 1
    python nest_samp.py -fc chime_cat2.json -o cat2_ALG_YMW16 -g ALG_YMW16 -upper 1
    python nest_samp.py -fc chime_cat1.json -o cat1_repeaters -g ALG_YMW16 -upper 1 -repeaters
"""

import numpy as np
from collections import Counter
import time
import pymultinest
import warnings
import sys
import argparse

from frb_util import *

# Initialize classes
dis = AstroDistribution()
cos = Cosmology()
er = EventRate()
tel = Telescope()
lf = Loadfiles()

dnu = 400.  # Reference bandwidth in MHz (CHIME bandwidth)


def lnlik(vpar):
    """Log-likelihood function (E_iso version)"""
    try:
        norm = np.zeros(vFOV.shape)
        for i in range(len(norm)):
            norm[i] = dis.Norm1D(vSN0[i], vBW[i], vNpol[i], vG[i], vTs[i],
                                 dnu, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        loglik_fdm = np.zeros(vN.shape)
        for i in range(len(vN)):
            loglik_fdm[i] = np.sum(
                dis.log_distr_fdmw(dnu, vLOGF_2d[i], vDME_2d[i], vLOGW_2d[i],
                                   vpar[1], vpar[2], vpar[3], vpar[4], vpar[5],
                                   gtype=fgt) - np.log(norm[i]))
        loglik_norm = np.sum(loglik_fdm)
        rho = np.zeros(vFOV.shape)
        for i in range(len(rho)):
            rho[i] = er.rate_2d(vSN0[i], vBW[i], vNpol[i], vG[i], vTs[i],
                                dnu, vpar[0], vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        loglik_poi = np.sum(er.log_dis_poi(rho, vN, vFOV, vTime))
        res = loglik_norm + loglik_poi
        return res
    except Exception as e:
        print(f'Numerical error: @ {vpar}, {e}')
        return -1e99


def myprior(cube, ndim, nparams):
    """Prior transformation for MultiNest"""
    for i in range(ndim):
        cube[i] = vpar_range[i, 0] + cube[i] * (vpar_range[i, 1] - vpar_range[i, 0])
    if bolupper:
        cube[0] = np.power(10., cube[0])
        cube[3] = np.log10(cube[3])


def myloglike(cube, ndim, nparams):
    """Log-likelihood wrapper for MultiNest"""
    cube2 = np.zeros(ndim)
    for i in range(0, ndim):
        cube2[i] = cube[i]
    res = lnlik(cube2)
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='E_iso Function Measurements for CHIME/FRB sample')
    parser.add_argument('-fc', action='store', dest='fcat', type=str,
                        help='Input CHIME FRB catalog file (JSON)')
    parser.add_argument('-o', action='store', dest='fout', type=str,
                        help='Output file prefix')
    parser.add_argument('-g', action='store', dest='fgt', type=str,
                        help='Host galaxy DM type: ALG_NE2001, ALG_YMW16, ETG, LTG_NE2001, LTG_YMW16')
    parser.add_argument('-upper', dest='bolupper', action='store_true', default=False,
                        help='Use flat prior for phis and logEiso0')
    parser.add_argument('-halo', dest='bolhalo', action='store_true', default=False,
                        help='Subtract DM contribution from dark matter halo (30 pc/cm3)')
    parser.add_argument('-repeaters', dest='repeaters', action='store_true', default=False,
                        help='Use only repeating FRBs')
    parser.add_argument('-oneoff', dest='oneoff', action='store_true', default=False,
                        help='Use only one-off FRBs')
    parser.add_argument('-bw', action='store', dest='default_bw', type=float, default=400.,
                        help='Default bandwidth in MHz for E_iso calculation (CHIME: 400)')

    args = parser.parse_args()
    fcat = args.fcat
    fout = args.fout
    fgt = args.fgt
    bolupper = args.bolupper
    bolhalo = args.bolhalo
    use_repeaters = args.repeaters
    use_oneoff = args.oneoff
    default_bw = args.default_bw

    # Load CHIME catalog
    if 'cat1' in fcat.lower() or 'cat1' in fcat:
        frb_cat = lf.LoadCHIMECat1(fcat)
    elif 'cat2' in fcat.lower() or 'cat2' in fcat:
        frb_cat = lf.LoadCHIMECat2(fcat)
    else:
        frb_cat = lf.LoadCHIMECat1(fcat)

    # Filter by repeater/one-off
    if use_repeaters:
        mask = frb_cat['is_repeater']
        print(f'Using {np.sum(mask)} repeating bursts')
    elif use_oneoff:
        mask = ~frb_cat['is_repeater']
        print(f'Using {np.sum(mask)} one-off bursts')
    else:
        mask = np.ones(len(frb_cat['fluence']), dtype=bool)
        print(f'Using all {np.sum(mask)} bursts')

    # Apply mask
    vFLUENCE = frb_cat['fluence'][mask]
    vDM = frb_cat['dm'][mask]
    vWIDTH = frb_cat['width'][mask]

    # Calculate excess DM
    if fgt.find('NE2001') >= 0:
        vDME = frb_cat['dm_exc_ne2001'][mask]
    else:
        vDME = frb_cat['dm_exc_ymw16'][mask]

    # Subtract halo DM contribution
    if bolhalo:
        vDME = vDME - 30.

    # Ensure positive values
    vDME = np.maximum(vDME, 0.1)

    # Calculate observed quantities
    vLOGF = np.log10(vFLUENCE)  # log10 fluence in Jy ms
    vLOGW = np.log10(np.maximum(vWIDTH, 0.001))  # log10 width in seconds

    # CHIME survey parameters (single survey)
    # CHIME/FRB: FoV ~ 200 deg^2, gain ~ 1.4 K/Jy, Tsys ~ 25K, BW ~ 400 MHz, Npol=2, SN0~10
    chime_fov = 200.0  # deg^2
    chime_gain = 1.4   # K/Jy
    chime_tsys = 25.0  # K
    chime_bw = 400.0   # MHz
    chime_npol = 2
    chime_sn0 = 10.0

    # Observation time (CHIME/FRB catalog 1: ~213 days, catalog 2: longer)
    if 'cat2' in fcat.lower():
        chime_time = 365.0 * 24.0  # ~1 year in hours
    else:
        chime_time = 213.0 * 24.0  # ~213 days in hours

    vN = np.array([len(vLOGF)])
    vSN0 = np.array([chime_sn0])
    vTs = np.array([chime_tsys])
    vG = np.array([chime_gain])
    vBW = np.array([chime_bw])
    vNpol = np.array([chime_npol])
    vFOV = np.array([chime_fov])
    vTime = np.array([chime_time])

    vLOGF_2d = [vLOGF]
    vDME_2d = [vDME]
    vLOGW_2d = [vLOGW]

    # Set prior ranges (E_iso version)
    # Parameters: [phis, alpha, logEiso_max, logEiso_min, mu_w, sigma_w]
    # CHIME data: E_iso ~ 10^37 to 10^42 erg, median ~ 10^39.9 erg
    # logEiso_max = characteristic/cut-off E_iso (where Schechter exponential kicks in)
    # logEiso_min = minimum E_iso (detection threshold)
    if bolupper:
        # Flat prior in log for phis, flat for logEiso0
        vpara = np.array([-3.0, -3.0, 39.0, 1e36, -1, 0.1])
        vparb = np.array([4.47, 1.1, 43.0, 1e40, 2, 1.0])
    else:
        vpara = np.array([1e-3, -3.0, 39.0, 37.0, -1, 0.1])
        vparb = np.array([3e4, 1.1, 43.0, 42.0, 2, 1.0])

    vpar_range = np.dstack((vpara.transpose(), vparb.transpose()))[0, :, :]

    print('------------ Parameter range -----------')
    print(vpar_range)
    a1 = time.perf_counter()
    print(f'Test likelihood: {myloglike(vpara, len(vpara), len(vpara))}')
    a2 = time.perf_counter()
    print(f'Time: {a1:.2f}, {a2:.2f}')
    print("Running Nest Sampling ...")

    # Create output directory
    import os
    os.makedirs('nest_out/samp', exist_ok=True)

    # Run MultiNest
    pymultinest.run(myloglike, myprior, len(vpara),
                    importance_nested_sampling=False,
                    resume=False,
                    verbose=True,
                    sampling_efficiency='model',
                    n_live_points=1000,
                    outputfiles_basename='nest_out/samp/' + fout)
