#!/usr/bin/env python3

import numpy as np
from scipy import interpolate
import time
import pymultinest
import warnings
import sys
import argparse

from frb_util import *
dis = AstroDistribution()
cos = Cosmology()
er = EventRate()
lf = Loadfiles()

dnu = 400.0  # 默认值，__main__ 中从 config 更新

sn0 = 10.
npol = 2.

# Two surveys
g1 = 0.7
bw1 = 300
tsys1 = 30
fov1 = 0.55

g2 = 0.05
bw2 = 300
tsys2 = 100
fov2 = 30

def lnlik(vpar):
    """能量版似然函数

    vpar: [phis, alpha, log_Es, log_E0, mu_w, sigma_w]
    """
    try:
        norm1 = dis.Norm1D_E(sn0, bw1, npol, g1, tsys1, dnu,
                             vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        norm2 = dis.Norm1D_E(sn0, bw2, npol, g2, tsys2, dnu,
                             vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        loglik_fdm1 = np.sum(
            dis.log_distr_efdmw(dnu, vLOGFLUX1, vDME1, vLOGW1,
                                vpar[1], vpar[2], vpar[3], vpar[4], vpar[5],
                                gtype=fgt) - np.log(norm1))
        loglik_fdm2 = np.sum(
            dis.log_distr_efdmw(dnu, vLOGFLUX2, vDME2, vLOGW2,
                                vpar[1], vpar[2], vpar[3], vpar[4], vpar[5],
                                gtype=fgt) - np.log(norm2))
        rho1 = er.rate_2d_E(sn0, bw1, npol, g1, tsys1, dnu,
                            vpar[0], vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        rho2 = er.rate_2d_E(sn0, bw2, npol, g2, tsys2, dnu,
                            vpar[0], vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        rho = np.array([rho1, rho2])
        loglik_poi = np.sum(er.log_dis_poi(rho, vN, vFOV, vT))
        res = loglik_fdm1 + loglik_fdm2 + loglik_poi
        return res
    except Exception:
        print('Numerical error: @', vpar)
        return -1e99

def myprior(cube, ndim, nparams):
    """从 config.json 读取先验范围（能量参数）"""
    cfg = config['prior']
    cube[0] = 10.0 ** (cfg['log_phis'][0] + cube[0] * (cfg['log_phis'][1] - cfg['log_phis'][0]))
    cube[1] = cfg['alpha'][0] + cube[1] * (cfg['alpha'][1] - cfg['alpha'][0])
    cube[2] = cfg['log_Es'][0] + cube[2] * (cfg['log_Es'][1] - cfg['log_Es'][0])
    cube[3] = cfg['log_E0'][0] + cube[3] * (cfg['log_E0'][1] - cfg['log_E0'][0])
    cube[4] = cfg['mu_w'][0] + cube[4] * (cfg['mu_w'][1] - cfg['mu_w'][0])
    cube[5] = cfg['sigma_w'][0] + cube[5] * (cfg['sigma_w'][1] - cfg['sigma_w'][0])

def myloglike(cube, ndim, nparams):
    cube2 = np.zeros(ndim)
    for i in range(0, ndim):
        cube2[i] = cube[i]
    res = lnlik(cube2)
    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='E_iso Measurements Program for Simulated FRBs')
    parser.add_argument('-f1', action='store', dest='simu1', type=str, help='Input Survey 1')
    parser.add_argument('-f2', action='store', dest='simu2', type=str, help='Input Survey 2')
    parser.add_argument('-o', action='store', dest='fout', type=str, help='Output file basename')
    parser.add_argument('-g', action='store', dest='fgt', type=str, help='Host galaxy type')
    parser.add_argument('--config', action='store', dest='config_path', type=str,
                        default='config.json', help='Path to config.json')
    args = parser.parse_args()
    simu1 = args.simu1
    simu2 = args.simu2
    fout = args.fout
    fgt = args.fgt

    # 加载全局配置
    config = load_config(args.config_path)
    dnu = config.get('analysis', {}).get('dnu', 400.0)
    cosmo_cfg = config.get('cosmology', {})
    cos = Cosmology(omegam=cosmo_cfg.get('Omega_m', 0.308),
                    omegal=cosmo_cfg.get('Omega_L', 0.692))

    cat1 = lf.LoadSimuData(simu1)
    vLOGFLUX1 = np.log10(cat1['S'])
    vLOGW1 = np.log10(cat1['W'])
    vDME1 = cat1['DMe']
    vdT1 = cat1['T']
    vN1 = len(vLOGFLUX1)
    vT1 = np.sum(vdT1)

    cat2 = lf.LoadSimuData(simu2)
    vLOGFLUX2 = np.log10(cat2['S'])
    vLOGW2 = np.log10(cat2['W'])
    vDME2 = cat2['DMe']
    vdT2 = cat2['T']
    vN2 = len(vLOGFLUX2)
    vT2 = np.sum(vdT2)

    vN = np.array([vN1, vN2])
    vFOV = np.array([fov1, fov2])
    vT = np.array([vT1, vT2])

    # 先验范围（从 config.json 读取）
    cfg = config['prior']
    vpara = np.array([cfg['log_phis'][0], cfg['alpha'][0], cfg['log_Es'][0],
                      cfg['log_E0'][0], cfg['mu_w'][0], cfg['sigma_w'][0]])
    vparb = np.array([cfg['log_phis'][1], cfg['alpha'][1], cfg['log_Es'][1],
                      cfg['log_E0'][1], cfg['mu_w'][1], cfg['sigma_w'][1]])

    vpar_range = np.dstack((vpara.transpose(), vparb.transpose()))[0, :, :]

    print('------------par range-----------')
    print(vpar_range)
    a1 = time.perf_counter()
    print(myloglike(vpara, len(vpara), len(vpara)))
    a2 = time.perf_counter()
    print(a1, a2)
    print("Running Nest Sampling ...")
    # run MultiNest
    import os
    output_dir = config['output']['nest_out_dir'] + 'simu/'
    os.makedirs(output_dir, exist_ok=True)
    pymultinest.run(myloglike, myprior, len(vpara),
                    importance_nested_sampling=False,
                    resume=False,
                    verbose=True,
                    sampling_efficiency='model',
                    n_live_points=1000,
                    outputfiles_basename=config['output']['nest_out_dir'] + 'simu/' + fout)
