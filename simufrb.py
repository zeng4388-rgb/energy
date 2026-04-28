#!/usr/bin/env python3
"""
simufrb.py - FRB simulation using E_iso Schechter function
Simulates mock FRB observations for CHIME-like surveys

Usage:
    python simufrb.py -ns 100 -phis 1000 -alpha -1.5 -logeisomax 50.0 -logeiso0 44.0 \
        -dnu 1000 -mu 0.1 -sig 0.3 -fgt ALG_YMW16 -ga 1.4 -npol 2 -bw 400 -ts 25 -sn0 10 -fov 200 \
        -out simu/simdat_chime.txt
"""

import numpy as np
import random as rd
from scipy import integrate
from scipy.interpolate import interp1d
import time
import warnings
import sys
import argparse

from frb_util import *

dis = AstroDistribution()
cos = Cosmology()
er = EventRate()
tel = Telescope()


def Simu_FRBs(phis, alpha, logeisomax, logeiso0, mu_w, sigma_w, dnu, ns,
              fov, npol, g, tsys, bw, sn0):
    """
    Simulate FRB observations using E_iso Schechter function
    Returns array with columns: fluence, width, time, DM_excess, threshold,
                                 logE_iso, z, DM_igm, DM_host, DM_src
    """
    res = np.zeros((ns, 10))
    ns0 = ns
    nt = 0
    lamda = er.rate_2d(sn0, bw, npol, g, tsys, dnu, phis, alpha,
                       logeiso0, logeisomax, mu_w, sigma_w) * fov
    print(f'Expected rate: {lamda:.2f} events/year')
    print(f'Simulating {ns} FRBs...')

    while ns > 0:
        # Sampling E_iso from Schechter function
        vlogeiso = np.arange(logeiso0 - 1., 55., (55. - logeiso0 + 1.) / 10000)
        vlik = dis.Schechter_log(vlogeiso, 1, alpha, logeisomax)
        vlogeiso = Sampling1D(vlogeiso, vlik, logeiso0, 54, ns0)
        vEiso = np.power(10., vlogeiso)

        # Sampling beam position (epsilon)
        vlnEps = np.random.uniform(-np.log(2), 0, ns0)
        vEps = np.exp(vlnEps)

        # Sampling redshift from comoving volume
        vZg = np.arange(0, 3.1, 3.1 / 10000)
        vlik = dis.Distribution_volume(vZg)
        vZ = Sampling1D(vZg, vlik, 0, 3.0, ns0)

        # Sampling intrinsic pulse width (log-normal)
        vlogW0 = np.arange(-0.5, 1.5, 2. / 10000)
        vlik = dis.dis_logw(vlogW0, mu_w, sigma_w)
        vlogW0 = Sampling1D(vlogW0, vlik, -0.4, 1.4, ns0)
        vW = np.power(10., vlogW0) * (1 + vZ)

        # Sampling DM of host galaxies
        vDMH0 = np.arange(0, 5001., 5001. / 10000)
        vlik = dis.Distribution_HostGalaxyDM(vDMH0, fgalaxy_type=fgt)
        vDMH0 = Sampling1D(vDMH0, vlik, 0, 5000, ns0)
        vDMH = vDMH0 * np.sqrt(dis.SFR(vZ)) / np.sqrt(dis.SFR(0))

        # Sampling DM of local sources
        vDMS = np.random.uniform(0, 50, ns0)

        # Calculate DM of IGM with redshift
        vDMI = cos.DispersionMeasure_IGM(vZ)

        # Sum up to extragalactic DM
        vDME = (vDMH + vDMS) / (1 + vZ) + vDMI

        # Sensitivity threshold (fluence in Jy ms)
        vft = tel.RMEq(sn0, g, tsys, npol, bw, vW)

        # Calculate observed fluence from E_iso
        vFluence = cos.Eiso_to_Fluence(vZ, vEiso * vEps, dnu=dnu)

        # Selection: detect if fluence > threshold
        detected = vFluence > vft
        nlen = np.sum(detected)

        if nlen > 0:
            if nlen > ns:
                nlen = ns

            # Sample arrival times
            larray = np.repeat(lamda, nlen)
            vT = rd.expovariate(larray[0]) if len(larray) == 1 else np.array(
                [rd.expovariate(l) for l in larray[:nlen]])

            idx = np.where(detected)[0][:nlen]
            res[nt:(nt + nlen), 0] = vFluence[idx]       # Fluence (Jy ms)
            res[nt:(nt + nlen), 1] = vW[idx]              # Observed width (ms)
            res[nt:(nt + nlen), 2] = vT[:nlen]            # Arrival time
            res[nt:(nt + nlen), 3] = vDME[idx]            # Excess DM
            res[nt:(nt + nlen), 4] = vft[idx]             # Threshold fluence
            res[nt:(nt + nlen), 5] = vlogeiso[idx]        # log10 E_iso
            res[nt:(nt + nlen), 6] = vZ[idx]              # Redshift
            res[nt:(nt + nlen), 7] = vDMI[idx]            # IGM DM
            res[nt:(nt + nlen), 8] = vDMH[idx]            # Host DM
            res[nt:(nt + nlen), 9] = vDMS[idx]            # Source DM
            nt = nt + nlen

        ns = ns - nlen
        pct = float(nt) / ns0 * 100
        print(f"{pct:.1f}% mock FRBs have been simulated.")

    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FRB Simulator (E_iso version)')
    parser.add_argument('-ns', action='store', dest='Ns', type=int, default=100,
                        help='Number of FRBs to simulate')
    parser.add_argument('-phis', action='store', dest='phis', type=float, default=1000,
                        help='Event rate density (Gpc^-3 yr^-1)')
    parser.add_argument('-alpha', action='store', dest='alpha', type=float, default=-1.5,
                        help='Power-law index of E_iso function')
    parser.add_argument('-logeisomax', action='store', dest='logeisomax', type=float, default=41.0,
                        help='log10 characteristic E_iso (erg)')
    parser.add_argument('-logeiso0', action='store', dest='logeiso0', type=float, default=38.0,
                        help='log10 minimum E_iso (erg)')
    parser.add_argument('-dnu', action='store', dest='dnu', type=float, default=400,
                        help='Reference bandwidth (MHz)')
    parser.add_argument('-mu', action='store', dest='mu', type=float, default=-0.5,
                        help='Mean of log pulse width distribution')
    parser.add_argument('-sig', action='store', dest='sigma', type=float, default=0.3,
                        help='Std of log pulse width distribution')
    parser.add_argument('-fgt', action='store', dest='fgt', type=str, default='ALG_YMW16',
                        help='Host galaxy DM type')
    parser.add_argument('-ga', action='store', dest='gain', type=float, default=1.4,
                        help='Telescope gain (K/Jy)')
    parser.add_argument('-npol', action='store', dest='npol', type=int, default=2,
                        help='Number of polarization channels')
    parser.add_argument('-bw', action='store', dest='bw', type=float, default=400,
                        help='Bandwidth (MHz)')
    parser.add_argument('-ts', action='store', dest='tsys', type=float, default=25,
                        help='System temperature (K)')
    parser.add_argument('-sn0', action='store', dest='sn0', type=float, default=10,
                        help='SNR detection threshold')
    parser.add_argument('-fov', action='store', dest='fov', type=float, default=200,
                        help='Field of view (deg^2)')
    parser.add_argument('-out', action='store', dest='output', default='simu/simdat.txt',
                        help='Output file')

    args = parser.parse_args()
    Ns = args.Ns
    phis = args.phis
    alpha = args.alpha
    logeisomax = args.logeisomax
    logeiso0 = args.logeiso0
    dnu = args.dnu
    mu = args.mu
    sigma = args.sigma
    fgt = args.fgt
    gain = args.gain
    npol = args.npol
    bw = args.bw
    Ts = args.tsys
    sn0 = args.sn0
    fov = args.fov
    output = args.output

    res = Simu_FRBs(phis, alpha, logeisomax, logeiso0, mu, sigma, dnu,
                    Ns, fov, npol, gain, Ts, bw, sn0)

    import os
    os.makedirs(os.path.dirname(output) if os.path.dirname(output) else '.', exist_ok=True)
    np.savetxt(output, res, delimiter=' ',
               header="#S W T DMe thres logEiso Z DMi DMh DMs", comments='')
    print(f'Saved {len(res)} simulated FRBs to {output}')
