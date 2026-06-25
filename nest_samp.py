#!/usr/bin/env python3

import numpy as np
from collections import Counter
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
tel = Telescope()
lf = Loadfiles()

dnu = 400.0  # 默认值，__main__ 中从 config 更新

def lnlik(vpar):
    """能量版似然函数

    vpar: [phis, alpha, log_Es, log_E0, mu_w, sigma_w]
    """
    try:
        # 归一化因子（能量版）
        norm = np.zeros(vFOV.shape)
        for i in range(len(norm)):
            norm[i] = dis.Norm1D_E(vSN0[i], vBW[i], vNpol[i], vG[i], vTs[i],
                                   dnu, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        # 逐事件似然（能量版，对 z 边际化）
        loglik_fdm = np.zeros(vN.shape)
        for i in range(len(vN)):
            loglik_fdm[i] = np.sum(
                dis.log_distr_efdmw(dnu, vLOGF_2d[i], vDME_2d[i], vLOGW_2d[i],
                                    vpar[1], vpar[2], vpar[3], vpar[4], vpar[5],
                                    gtype=fgt) - np.log(norm[i]))
        loglik_norm = np.sum(loglik_fdm)
        # 泊松似然（能量版事件率）
        rho = np.zeros(vFOV.shape)
        for i in range(len(rho)):
            rho[i] = er.rate_2d_E(vSN0[i], vBW[i], vNpol[i], vG[i], vTs[i],
                                  dnu, vpar[0], vpar[1], vpar[2], vpar[3],
                                  vpar[4], vpar[5])
        loglik_poi = np.sum(er.log_dis_poi(rho, vN, vFOV, vTime))
        res = loglik_norm + loglik_poi
        return res
    except Exception:
        print('Numerical error: @', vpar)
        return -1e99

def myprior(cube, ndim, nparams):
    """从 config.json 读取先验范围（能量参数）

    cube[0] 输出 linear phis（下游 rate_2d_E 需要线性值）
    cube[2], cube[3] 输出 log10(E*), log10(E0)
    """
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
    parser = argparse.ArgumentParser(description='E_iso Measurements Program for FRB sample')
    parser.add_argument('-fc', action='store', dest='fcat', type=str,
                        help='Input FRB catalog file (txt or fits)')
    parser.add_argument('-fs', action='store', dest='fsvy', type=str,
                        help='Input Survey Information file')
    parser.add_argument('-o', action='store', dest='fout', type=str,
                        help='Output file basename')
    parser.add_argument('-g', action='store', dest='fgt', type=str,
                        help='Host galaxy type (ETG, LTG_NE2001, LTG_YMW16, ALG_NE2001, ALG_YMW16)')
    parser.add_argument('-halo', action='store_true', dest='bolhalo',
                        help='Bool option: removing the DM from dark halo')
    parser.add_argument('--fits', action='store_true', dest='use_fits',
                        help='Load catalog from FITS file instead of txt')
    parser.add_argument('--config', action='store', dest='config_path', type=str,
                        default='config.json', help='Path to config.json')

    args = parser.parse_args()
    fcat = args.fcat
    fsvy = args.fsvy
    fout = args.fout
    fgt = args.fgt
    bolhalo = args.bolhalo
    use_fits = args.use_fits

    # 加载全局配置
    config = load_config(args.config_path)
    dnu = config.get('analysis', {}).get('dnu', 400.0)
    halo_dm = config.get('analysis', {}).get('halo_dm', 0.0)
    # 用 config 中的宇宙学参数更新全局 cos 实例
    cosmo_cfg = config.get('cosmology', {})
    cos = Cosmology(omegam=cosmo_cfg.get('Omega_m', 0.308),
                    omegal=cosmo_cfg.get('Omega_L', 0.692))

    # 如果未指定 fcat/fsvy，从 config 读取默认路径
    if fcat is None:
        fcat = config['data']['catalog_path']
        # 仅在用户未显式指定 --fits 时，根据文件扩展名自动检测
        if not use_fits:
            use_fits = fcat.endswith('.fits') or fcat.endswith('.fit')
    if fsvy is None:
        fsvy = config['data']['survey_info_path']

    # 加载数据
    if use_fits:
        vF, vW, vDM_obs, vDM_ne2001, vDM_ymw16, vSVY = lf.LoadFitsCatalog(fcat)
        # 根据银河系 DM 模型选择对应的河外 DM
        if fgt and fgt.find('NE2001') >= 0:
            vDME = vDM_ne2001  # 已经是 DM_obs - DM_MW(NE2001)
        else:
            vDME = vDM_ymw16   # 已经是 DM_obs - DM_MW(YMW16)
        if bolhalo:
            vDME = vDME - halo_dm
        vLOGF = np.log10(vF / vW)   # flux = fluence / width [Jy]
        vLOGW = np.log10(vW)
    else:
        frb_cat = lf.LoadCatalogue(fcat)
        if fgt and fgt.find('NE2001') >= 0:
            vDME = frb_cat['DM'] - frb_cat['DM_NE2001']
        else:
            vDME = frb_cat['DM'] - frb_cat['DM_YMW16']
        if bolhalo:
            vDME = vDME - halo_dm
        vLOGF = np.log10(frb_cat['S'])
        vLOGW = np.log10(frb_cat['W'])
        vSVY = frb_cat['SURVEY']

    if fgt and fgt.find('ETG') >= 0:
        fgt = 'ETG'

    # 加载巡天信息
    svy_info = lf.LoadSvyInfo(fsvy)

    # 按巡天分组
    cts = Counter(vSVY)

    # 自动检测巡天列表（从 tel_svy.txt 读取）
    survey_names = svy_info['SURVEY']
    vN = np.array([cts.get(s, 0) for s in survey_names])
    vSN0 = svy_info['SN0']
    vTs = svy_info['Tsys']
    vG = svy_info['Gain']
    vBW = svy_info['BW']
    vNpol = svy_info['Npol']
    vFOV = svy_info['FOV']
    vTime = svy_info['TIME']

    # 按巡天名称分组数据
    vLOGF_2d = [vLOGF[vSVY == s] for s in survey_names]
    vDME_2d = [vDME[vSVY == s] for s in survey_names]
    vLOGW_2d = [vLOGW[vSVY == s] for s in survey_names]

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
    output_dir = config['output']['nest_out_dir'] + 'samp/'
    os.makedirs(output_dir, exist_ok=True)
    pymultinest.run(myloglike, myprior, len(vpara),
                    importance_nested_sampling=False,
                    resume=False,
                    verbose=True,
                    sampling_efficiency='model',
                    n_live_points=1000,
                    outputfiles_basename=config['output']['nest_out_dir'] + 'samp/' + fout)
