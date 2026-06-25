import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import scipy.special as spf
import json
from pathlib import Path

def load_config(config_path='config.json'):
    """加载全局配置文件（使用相对路径，基于脚本所在目录）"""
    with open(Path(__file__).parent / config_path) as f:
        return json.load(f)

class Cosmology:
    _cache = {}  # (omegam, omegal, omegab) -> instance，避免重复积分

    def __init__(self, omegam=0.308, omegal=0.692, omegab=0.0484):
        key = (omegam, omegal, omegab)
        cached = Cosmology._cache.get(key)
        if cached is not None:
            self.__dict__.update(cached.__dict__)
            return
        self.Omega_m = omegam
        self.Omega_b = omegab
        self.Omega_L = omegal
        self.c = 2.9979245800e10
        self.pc2cm = 3.08567758e18
        self.km2cm = 1e5
        self.Mpc2cm = 3.08567758e24
        self.Gpc2cm = 3.08567758e27
        self.Jy2CGS = 1e-23
        self.MHz2Hz = 1e6
        self.Jyms2CGS = 1e-26
        self.h0 = 0.6781
        self.H0 = self.h0 * 100 * self.km2cm / self.Mpc2cm
        self.Rhoc = 1.88 * self.h0 * self.h0 * 1e-29    #The critical density of universe in gram/cm^3
        self.Nc = self.Rhoc / 1.6726e-24       #The number density of universe in Hydrogen atom, in units of 1/cm^3
        self.f_IGM = 0.83
        self.z0 = 0.8

        self.vz = np.arange(-7, 6, 0.03)
        self.vz = np.power(10., self.vz)
        self.vz[0] = 0
        func = lambda z: self.c / self.H(z)
        self.vd = [ integrate.quad(func, 0, zv)[0] for zv in self.vz ]
        self.cd_interp = interp1d(self.vz, self.vd)
                
        func2 = lambda z: (1+z) * self.c / self.H(z) * self.Omega_b * self.Nc
        self.vdm = [ integrate.quad(func2, 0, zv)[0] for zv in self.vz ]
        self.dmigm = interp1d(self.vz, self.vdm)
        
        self.vld = self.Luminosity_Distance(self.vz)
        self.Ld2z = interp1d(self.vld, self.vz)
        self.Ld1 = self.Luminosity_Distance(1.0)
        self.Cd1 = self.Comoving_Distance(1.0)
        Cosmology._cache[key] = self
        
    def E(self, z):
        """
        logarithmic time derivative of scale factor
        """
        return np.sqrt(self.Omega_m * np.power(1 + z, 3.0) + self.Omega_L)
    
    def H(self, z):
        """
        calculate the Hubble ratio
        z: redshifts
        """
        return self.H0 * self.E(z)
    
    def Comoving_Distance(self, z):
        """
        calculate the comoving distance
        z: redshifts
        """
        return self.cd_interp(z)

    def dVdOdz(self,z):
        """
        calculate the diffrential comoving volume dV/dz/dOmega in units of Gpc^3
        z: redshifts
        """
        drdz = self.c / self.H(z)
        cd2 = self.Comoving_Distance(z)
        cd2 = cd2 * cd2
        dcv = cd2 * drdz / (self.Gpc2cm * self.Gpc2cm * self.Gpc2cm)
        return dcv
    
    def Luminosity_Distance(self, z):
        """
        calculate the luminosity distance in Mpc
        z: redshifts
        """
        dl = (1 + z) * self.Comoving_Distance(z)
        return dl

    def Luminosity(self, z, f=1., dnu=1000.):
        """
        calculate the intrinsic luminosity from flux
        f: flux in unit of Jy
        dnu: intrinsic spectral width
        z: redshifts
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        lum = f * self.Jy2CGS * dnu * self.MHz2Hz * 4 * np.pi * ld2
        return lum

    def Energy(self, z, flu=1.0, dnu=1000.):
        """
        calculate the intrinsic energy from fluence
        flu: fluence in units of Jy ms
        dnu: intrinsic spectral width, default at 1000 MHz
        z: redshift
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        ener = flu / (1+z) * self.Jyms2CGS * dnu * self.MHz2Hz * 4 * np.pi * ld2
        return ener
    
    def LuminosityDistance_to_z(self, ld):
        """
        converting luminosity distance to redshift
        """
        return self.Ld2z(ld)
    
    def Luminosity_Distance_dimless(self, z):
        """
        calculate the dimension less luminosity distance
        z: redshifts
        """
        dl = (1 + z) * self.Comoving_Distance(z) / (2. * self.Cd1)
        return dl
    
    def DispersionMeasure_IGM(self, z, chi=7./8):
        """
        calculate dispersion measure of intergalactic medium by integrating redshift
        """
        return self.dmigm(z) * self.f_IGM * chi / self.Mpc2cm * 1e6
        
    def Luminosity_to_Flux(self, z, lum, dnu=1000):
        """
        calculate the flux density observed at a given luminosity
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld*ld
        flux = lum/4/np.pi/ld2/dnu/self.MHz2Hz/self.Jy2CGS
        return flux
    
    def Energy_to_Flu(self, z, ener, dnu=1000):
        """
        calculate the observed fluence when knowing intrinsic energy
        """
        ld = self.Luminosity_Distance(z)
        ld2 = ld*ld
        flu = ener*(1+z)/4/np.pi/ld2/dnu/self.MHz2Hz/self.Jyms2CGS
        return flu

    def DMeq(self, z, dme, dmhost):
        """
        The DM equation which is solved to get redshift value
        """
        z = np.maximum(z, 1e-7)   # fsolve 可能试探负值，保护插值器
        dmi = self.DispersionMeasure_IGM(z)
        dmh = dmhost/(1+z)
        return dmi + dmh - dme

    def GetZ(self, dme, dmhost):
        """
        solving the differential equation to get redshift
        """
        z = fsolve(self.DMeq, self.z0, args=(dme, dmhost))
        z = float(z)
        return z

class Telescope:
    def __init__(self):
        self.MHz2Hz = 1e6
        self.ms2s = 1e-3

    def RMEq(self, snr, g, tsys, npol, bw, w):
        s = snr*tsys/g/np.sqrt(npol*bw*self.MHz2Hz*w*self.ms2s)
        return s

class AstroDistribution:
    def __init__(self):
        self.tel = Telescope()
        self.cos = Cosmology()
        self.vpar_etg = np.array([0.001713, 1.099, 0.2965, 0.01246, 1.055, 0.7262])
        self.vpar_ltg_ne2001 = np.array([0.01715, 1.062, 0.5202, 0.00416, 0.7227, 1.151])
        self.vpar_ltg_ymw16 = np.array([0.01561, 0.759, 0.3013, 0.01889, 1.042, 0.5791])
        self.vpar_alg_ne2001 = np.array([0.005485, 0.8665, 1.009, 0.01406, 1.069, 0.5069])
        self.vpar_alg_ymw16 = np.array([0.01199, 0.7597, 0.3082, 0.01735, 1.048, 0.6025])
        self.Zmax = 5
        self.Zmin = 2e-6
        self.DMsmax = 50.
        self.Wmax = 20
        self.Wmin = 0.05

    def Schechter_log(self, logl, phis, alpha, logls):
        """
        The Schecheter luminosity function per logarithmic luminosity
        logl: the logarithmic luminosity
        phis: normalization constant
        logls: the cut-off luminosity
        alpha: power index
        """
        l = np.power(10., logl)
        ls = np.power(10., logls)
        phi = np.log(10) * phis * np.power(l / ls, (alpha + 1)) * np.exp(-l / ls)
        return phi

    def log_IntBeam(self, logl, alpha, logls, logl0):
        ratio = np.power(10., logl-logls)
        lik0 = gammainc(alpha+1, ratio) - gammainc(alpha+1, 2*ratio)
        lik = lik0/np.log(2)      # Beam efficiency from 50% to 100%
        lik = np.where(np.isfinite(lik) & (lik > 0), lik, 1e-199)
        loglik = np.log(lik)
        loglik[logl < logl0] = -1e99
        return loglik

    def IntLum(self, eps, alpha, logls, logl0):
        with np.errstate(invalid='ignore'):
             ratio = np.power(10., logl0-logls)
        return gammainc(alpha+1, ratio/eps)

    # --- E_iso 能量函数别名（数学形式与光度函数完全相同） ---
    def Schechter_E_log(self, loge, phis, alpha, logEs):
        """Schechter 能量函数 per dex（与 Schechter_log 数学形式相同）"""
        return self.Schechter_log(loge, phis, alpha, logEs)

    def IntE(self, eps, alpha, logEs, logE_min):
        """Schechter 能量累积积分（与 IntLum 数学形式相同）"""
        return self.IntLum(eps, alpha, logEs, logE_min)

    def log_IntBeam_E(self, loge, alpha, logEs, logE0):
        """波束卷积后的 Schechter 能量函数（与 log_IntBeam 数学形式相同）"""
        return self.log_IntBeam(loge, alpha, logEs, logE0)
    
    def Distribution_Local_galaxy_DM(self, dmv, vpar):
        """
        General form of DM distribution function of host galaxies, using double gaussian function in logarithmic DM
        vpar[i] is the i-th parameter of this distribution function
        """
        val = vpar[0] * np.exp(-np.power((dmv - vpar[1]) / vpar[2], 2.)) \
              + vpar[3] * np.exp(-np.power((dmv - vpar[4]) / vpar[5], 2.))
        return val

    def ThetaFunc(self, x):
        """
        Normalize the values positive
        """
        return 0.5 * (np.sign(x) + 1)
            
    def func_gaussian(self, dmv, vpar):
        """
        default gaussina distribution for non-galaxy case
        """
        dmoff = dmv - vpar[0]
        sig = vpar[1]
        sig = sig * sig
        return np.exp(-0.5 * dmoff * dmoff / sig) * self.ThetaFunc(dmv)
    
    def func_uniform(self, dmv, vpar):
        pdf = np.ones(dmv.shape)
        pdf[dmv >= vpar[1]] = 0
        pdf[dmv <= vpar[0]] = 0
        return pdf
    
    def Distribution_HostGalaxyDM(self, dmv0, fgalaxy_type=None, vpar=np.array([0,50])):
        """
        DM distribution functions of different type of host galaxies
        ETG: early-tyep galaxies
        LTG: late-type galaxies
        ALG: all the galaxies
        NE2001: the referenced galaxy electron density using NE2001 model
        YMW16: the referenced galaxyt electron density using YMW16 model
        """
        ind = dmv0 < 0
        if not fgalaxy_type: 
            fgalaxy_type = self.func_gaussian
        elif fgalaxy_type == 'ETG':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_etg)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'LTG_NE2001':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), 
                    self.vpar_ltg_ne2001)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'LTG_YMW16':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), 
                    self.vpar_ltg_ymw16)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'ALG_NE2001':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_alg_ne2001)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        elif fgalaxy_type == 'ALG_YMW16':
            dmv = dmv0.copy()
            dmv[dmv0<1e-9] = np.ones(dmv[dmv0<1e-9].shape)*1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_alg_ymw16)
            res[dmv0<=0] = np.zeros(res[dmv0<=0].shape)
            return res
        else:
            return fgalaxy_type(dmv0,vpar)
    
    
    def log_Distribution_HostGalaxyDM(self, dmv, fgalaxy_type=None, vpar=np.array([0,50])):
        if not fgalaxy_type:
            dmoff = dmv[dmv > 0] - vpar[0]
            sig = vpar[1]
            sig = sig * sig
            res = dmv.copy()
            res[dmv < 0] = -1e99
            res[dmv > 0] = -0.5 * dmoff * dmoff / sig
            return res
        else:
            val = self.Distribution_HostGalaxyDM(dmv, fgalaxy_type=fgalaxy_type)
            ind=val <= 0
            indv=val>0
            val[indv] = np.log(val[indv])
            val[ind] = np.ones(val[ind].shape) * -1e99
            return val
    
    def SFR(self, z):
        """
        Star-forming history, the values taken from Hopkins & Beacom (2016)
        """
        sfr = (0.017 + 0.13 * z)/(1 + np.power(z/3.3, 5.3))
        return sfr 

    def kappa(self, z):
        """
        Normalized SFH from redshift of z to redshift of 0 (nearby universe)
        """
        return np.sqrt(self.SFR(0)/self.SFR(z))

    def Distribution_volume(self, z):
        """
        Differential comoving volume
        """
        r = self.cos.Comoving_Distance(z)/self.cos.Comoving_Distance(1)
        pv = r * r / self.cos.E(z)
        return pv
    
    def log_Distribution_volume(self, z):
        """
        Logarithimc differential comoving volume
        """
        if isinstance(z, np.ndarray):
            ind = z<0
            pv = np.log(self.Distribution_volume(z))
            pv[ind] = -1e99
            return pv
        else:
            if z<0:
                return -1e99
            else:
                return np.log(self.Distribution_volume(z))

    def IntDMsrc(self, u1, u2, vpar):
        """
        Analytic integral when marginalizing distribution of DMsrc
        """
        a1 = vpar[0]
        b1 = vpar[1]
        c1 = vpar[2]
        a2 = vpar[3]
        b2 = vpar[4]
        c2 = vpar[5]
        k1 = np.power(10.,b1)*a1*c1*np.exp(c1*c1*np.log(10)*np.log(10)/4)
        k2 = np.power(10.,b2)*a2*c2*np.exp(c2*c2*np.log(10)*np.log(10)/4)
        q1 = (c1*c1*np.log(10)*np.log(10)+b1*np.log(100)-2*np.log(u1))/c1/np.log(100)
        q2 = (c1*c1*np.log(10)*np.log(10)+b1*np.log(100)-2*np.log(u2))/c1/np.log(100)
        q3 = (c2*c2*np.log(10)*np.log(10)+b2*np.log(100)-2*np.log(u1))/c2/np.log(100)
        q4 = (c2*c2*np.log(10)*np.log(10)+b2*np.log(100)-2*np.log(u2))/c2/np.log(100)
        int_h = np.log(10)*np.sqrt(np.pi)/2 * (k1*(spf.erf(q2)-spf.erf(q1))+k2*(spf.erf(q4)-spf.erf(q3)))
        int_hs = int_h/self.DMsmax   #   integral including uniform DMsrc
        return int_hs

    def log_IntDMsrc(self, u1, u2, gtype=None):
        """
        Logarithmic integrals of above marginalization in different galaxy type cases
        """
        if not gtype:
            gtype = self.func_gaussian
            res = np.log(np.maximum(self.IntDMsrc(u1, u2, np.array([0, 50])), 1e-300))
            return res
        elif gtype == 'ETG':
            res = np.log(np.maximum(self.IntDMsrc(u1, u2, self.vpar_etg), 1e-300))
            return res
        elif gtype == 'LTG_NE2001':
            res = np.log(np.maximum(self.IntDMsrc(u1, u2, self.vpar_ltg_ne2001), 1e-300))
            return res
        elif gtype == 'LTG_YMW16':
            res = np.log(np.maximum(self.IntDMsrc(u1, u2, self.vpar_ltg_ymw16), 1e-300))
            return res
        elif gtype == 'ALG_NE2001':
            res = np.log(np.maximum(self.IntDMsrc(u1, u2, self.vpar_alg_ne2001), 1e-300))
            return res
        elif gtype == 'ALG_YMW16':
            res = np.log(np.maximum(self.IntDMsrc(u1, u2, self.vpar_alg_ymw16), 1e-300))
            return res
        else:
            return -1e99
           
    def dis_logw(self, logw0, mu, sigma):
        a = 1./np.sqrt(2.*np.pi*sigma*sigma)
        b = np.exp(-(logw0-mu)*(logw0-mu)/2/sigma/sigma)
        return a*b

    def log_dis_logw(self, logw0, mu, sigma):
        a = 2*np.pi*sigma*sigma
        b = -(logw0-mu)*(logw0-mu)/2/sigma/sigma
        return b-1/2.*np.log(a)

    def log_distr_fdmwz(self, dnu, logflux, dme, logw, z, alpha, logls, logl0, mu, sigma, gtype=None):
        """
        Logarithmic joint distribution function of flux, DM and redshift
        """
        flux = np.power(10., logflux)
        logl = np.log10(self.cos.Luminosity(z, f=flux, dnu=dnu))
        logint1 = self.log_IntBeam(logl, alpha, logls, logl0)
        #print logint1
        logw0 = logw - np.log10(1+z)
        logfw = self.log_dis_logw(logw0, mu, sigma)
        logfz = self.log_Distribution_volume(z)
        dmi = self.cos.DispersionMeasure_IGM(z)
        u1 = (dme-dmi)*(1+z)*self.kappa(z)
        u2 = ((dme-dmi)*(1+z)-self.DMsmax)*self.kappa(z)
        #print u2
        logint2 = np.ones(u2.shape) * (-1e99)
        ind = u2 > 0
        logint2[ind] = self.log_IntDMsrc(u1[ind],u2[ind],gtype=gtype)
        #print logint2
        loglikv = logint1 + logfz + logfw + logint2 + np.log(1+z)
        return loglikv

    def log_distr_efdmwz(self, dnu, logflux, dme, logw, z, alpha, logEs, logE0, mu, sigma, gtype=None):
        """能量版联合概率 p(E_iso, DM_host | z)

        将观测流量 logflux 和脉冲宽度 logw 耦合为 E_iso，再用 Schechter 能量函数评估。
        物理关系：E_iso = S * w * dnu * 4π * D_L² / (1+z)
        """
        flux = np.power(10., logflux)
        w_ms = np.power(10., logw)       # 脉冲宽度 [ms]
        fluence = flux * w_ms            # fluence = S * w [Jy·ms]
        loge = np.log10(self.cos.Energy(z, flu=fluence, dnu=dnu))
        logint1 = self.log_IntBeam_E(loge, alpha, logEs, logE0)
        logw0 = logw - np.log10(1+z)
        logfw = self.log_dis_logw(logw0, mu, sigma)
        logfz = self.log_Distribution_volume(z)
        dmi = self.cos.DispersionMeasure_IGM(z)
        u1 = (dme-dmi)*(1+z)*self.kappa(z)
        u2 = ((dme-dmi)*(1+z)-self.DMsmax)*self.kappa(z)
        logint2 = np.ones(u2.shape) * (-1e99)
        ind = u2 > 0
        logint2[ind] = self.log_IntDMsrc(u1[ind],u2[ind],gtype=gtype)
        loglikv = logint1 + logfz + logfw + logint2 + np.log(1+z)
        return loglikv

    def log_distr_fdmw(self, dnu, logflux, dme, logw, alpha, logls, logl0, mu, sigma, gtype=None):
        """
        Logarithmic joint distribution function of flux and DM after maginalization of redshift
        :param dnu: intrinsic bandwidth, i.e. 1GHz
        :param logflux: logarithmic flux
        :param dme: dme
        :param logls: logarithmic cut-off luminosity
        :param alpha: lf index
        :param logl0: minimum of lf L
        :param fgalaxy_type: galaxy type
        :return logarithmic likelihood:
        """
        #stepdms = 100/1000.
        #vdms = np.arange(0, 100, stepdm)
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 1000
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        lik = 0
        for z in vz:
            likv = np.exp(self.log_distr_fdmwz(dnu, logflux, dme, logw, z, alpha, logls, logl0, mu, sigma, gtype=gtype))
            lik += z * stepz * likv
        ind = lik > 0
        ind2 = lik <= 0
        loglik = lik.copy()
        loglik[ind] = np.log(lik[ind])
        loglik[ind2] = np.ones(loglik[ind2].shape) * -1e99
        return loglik

    def log_distr_efdmw(self, dnu, logflux, dme, logw, alpha, logEs, logE0, mu, sigma, gtype=None):
        """能量版联合概率，对红移 z 边际化

        对于有红移的事件，可直接使用 log_distr_efdmwz；本函数用于仅有 DM 的事件。
        """
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 200
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        lik = 0
        for z in vz:
            likv = np.exp(self.log_distr_efdmwz(dnu, logflux, dme, logw, z, alpha, logEs, logE0, mu, sigma, gtype=gtype))
            lik += z * stepz * likv
        ind = lik > 0
        ind2 = lik <= 0
        loglik = lik.copy()
        loglik[ind] = np.log(lik[ind])
        loglik[ind2] = np.ones(loglik[ind2].shape) * -1e99
        return loglik
    
    def Norm1D(self, sn0, bw, npol, g, tsys, dnu, alpha, logls, logl0, mu, sigma):
        """
        Normalization factor for dimensionless likelihood
        """
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        stepeps = (1-0.5) / 200.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.Wmax) - np.log10(self.Wmin)) / 1000.
        vlogw = np.arange(np.log10(self.Wmin), np.log10(self.Wmax), steplogw)
        nf = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            lt = self.cos.Luminosity(z, ft, dnu)
            loglt = np.log10(lt)
            ind = loglt < logl0
            loglt[ind] = logl0
            int_eps = np.zeros(loglt.shape)
            for i in np.arange(len(loglt)):
                int_eps[i] = np.sum(self.IntLum(veps, alpha, logls, loglt[i]*np.ones(veps.shape))/veps/np.log(2)*stepeps)
            int_w = np.sum(int_eps*self.dis_logw(vlogw, mu, sigma)*steplogw)
            fz = self.Distribution_volume(z)
            nf += z*stepz*fz*int_w
        if nf <= 0:
            nf = 1e-199
        return nf

    def Norm1D_E(self, sn0, bw, npol, g, tsys, dnu, alpha, logEs, logE0, mu, sigma):
        """能量版归一化因子

        检测阈值从 L_min 转换为 E_min：E_min = S_min * w_obs * dnu * 4π * D_L² / (1+z)
        """
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 200.
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        stepeps = (1-0.5) / 50.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.Wmax) - np.log10(self.Wmin)) / 100.
        vlogw = np.arange(np.log10(self.Wmin), np.log10(self.Wmax), steplogw)
        nf = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            fw = ft * vw   # F_min = S_min * w_obs [Jy·ms]
            et = self.cos.Energy(z, flu=fw, dnu=dnu)
            loget = np.log10(et)
            loget = np.where(np.isfinite(loget) & (loget >= logE0), loget, logE0)
            int_eps = np.zeros(loget.shape)
            for i in np.arange(len(loget)):
                int_eps[i] = np.sum(self.IntE(veps, alpha, logEs, loget[i]*np.ones(veps.shape))/veps/np.log(2)*stepeps)
            int_w = np.sum(int_eps*self.dis_logw(vlogw, mu, sigma)*steplogw)
            fz = self.Distribution_volume(z)
            nf += z*stepz*fz*int_w
        if nf <= 0:
            nf = 1e-199
        return nf

class EventRate:
    def __init__(self):
        self.cos = Cosmology()
        self.ad = AstroDistribution()
        self.tel = Telescope()
        self.s2h = 1./3600
        self.yr2hr = 365*24.
        self.rad2deg2 = 3282.806350011744
        self.Gpc2Mpc = 1e3
 
    def rate_2d(self, sn0, bw, npol, g, tsys, dnu, phis, alpha, logls, logl0, mu, sigma):
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        stepeps = (1-0.5)/ 200.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 1000.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        rho = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            lt = self.cos.Luminosity(z, ft, dnu)
            loglt = np.log10(lt)
            ind = loglt < logl0
            loglt[ind] = logl0
            int_eps = np.zeros(loglt.shape)
            for i in range(len(loglt)):
                int_eps[i] = np.sum(phis*self.ad.IntLum(veps, alpha, logls, loglt[i]*np.ones(veps.shape))/veps/np.log(2)*stepeps)
            int_w = np.sum(int_eps*self.ad.dis_logw(vlogw, mu, sigma)*steplogw)
            fz = self.cos.dVdOdz(z)/(1+z)
            rho += z*stepz*fz*int_w
        rho_deg = rho/self.rad2deg2/self.yr2hr
        return rho_deg

    def rate_2d_E(self, sn0, bw, npol, g, tsys, dnu, phis, alpha, logEs, logE0, mu, sigma):
        """能量版事件率密度 rho_deg [deg^-2 hr^-1]

        检测阈值从 L_min 转换为 E_min：E_min = S_min * w_obs * dnu * 4π * D_L² / (1+z)
        """
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 200.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        stepeps = (1-0.5)/ 50.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 100.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        rho = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            fw = ft * vw   # F_min = S_min * w_obs [Jy·ms]
            et = self.cos.Energy(z, flu=fw, dnu=dnu)
            loget = np.log10(et)
            loget = np.where(np.isfinite(loget) & (loget >= logE0), loget, logE0)
            int_eps = np.zeros(loget.shape)
            for i in range(len(loget)):
                int_eps[i] = np.sum(phis*self.ad.IntE(veps, alpha, logEs, loget[i]*np.ones(veps.shape))/veps/np.log(2)*stepeps)
            int_w = np.sum(int_eps*self.ad.dis_logw(vlogw, mu, sigma)*steplogw)
            fz = self.cos.dVdOdz(z)/(1+z)
            rho += z*stepz*fz*int_w
        rho_deg = rho/self.rad2deg2/self.yr2hr
        return rho_deg

    def log_dis_poi(self, rho, N, Omega, T):
        lamda = rho*Omega*T
        ind = lamda <= 0
        lamda[ind] = 1e-199
        loglik = N*np.log(lamda)-lamda-spf.gammaln(N+1)
        return loglik

    def Sens(self, sn0, g, tsys, npol, bw, mu, sigma):
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 1000.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        ints0 = 0
        intz = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            ints0 += z*stepz*np.sum(ft*self.ad.dis_logw(vlogw, mu, sigma)*steplogw)
            intz += z*stepz
        smin = ints0/intz
        return smin

    def Rate(self, logft, dnu, phis, alpha, logls, logl0, mu, sigma):
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        stepeps = (1-0.5) / 200.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 1000.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        rho = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = np.power(10, logft)
            lt = self.cos.Luminosity(z, ft, dnu)
            loglt = np.log10(lt)
            #ind = loglt < logl0
            #loglt[ind] = logl0
            if loglt < logl0:
                loglt = logl0
            int_eps = np.sum(phis*self.ad.IntLum(veps, alpha, logls, loglt*np.ones(veps.shape))/veps/np.log(2)*stepeps)
            int_w = np.sum(int_eps*self.ad.dis_logw(vlogw, mu, sigma)*steplogw)
            fz = self.cos.dVdOdz(z)/(1+z)
            rho += z*stepz*fz*int_w
        rho_deg = rho/self.rad2deg2/self.yr2hr
        return rho_deg

    def Rate_E(self, logft, dnu, phis, alpha, logEs, logE0, mu, sigma):
        """能量版事件率（对单一 flux 阈值）"""
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 200.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        stepeps = (1-0.5) / 50.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 100.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        rho = 0
        for z in vz:
            vw = np.power(10, vlogw)*(1+z)
            ft = np.power(10, logft)
            fw = ft * vw   # fluence = flux * width [Jy·ms]
            et = self.cos.Energy(z, flu=fw, dnu=dnu)
            loget = np.log10(et)
            loget = np.where(np.isfinite(loget) & (loget >= logE0), loget, logE0)
            int_eps = np.sum(phis*self.ad.IntE(veps, alpha, logEs, loget*np.ones(veps.shape))/veps/np.log(2)*stepeps)
            int_w = np.sum(int_eps*self.ad.dis_logw(vlogw, mu, sigma)*steplogw)
            fz = self.cos.dVdOdz(z)/(1+z)
            rho += z*stepz*fz*int_w
        rho_deg = rho/self.rad2deg2/self.yr2hr
        return rho_deg

    def Rfrb(self, phis, alpha, logls, loglmin):
        ratio = np.power(10., loglmin-logls)
        return phis*gammainc(alpha+1, ratio)

    def Rfrb_E(self, phis, alpha, logEs, logEmin):
        """能量版事件率（对单一能量阈值）"""
        ratio = np.power(10., logEmin-logEs)
        return phis*gammainc(alpha+1, ratio)   

class Loadfiles:
    #def __init__(self):

    def LoadCatalogue(self, fname):
        """遗留 TXT 加载器，仅兼容旧版手工整理的 TXT 格式
        注意：与 preprocess_catalog.py 输出的列名不兼容，主流程应使用 LoadFitsCatalog
        """
        cat = np.loadtxt(fname, dtype=str)
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0,i]] = cat[1:, i]

        cat2['S'] = np.array(cat2['S'], dtype=float)
        cat2['Seu'] = np.array(cat2['Seu'], dtype=float)
        cat2['Sel'] = np.array(cat2['Sel'], dtype=float)
        cat2['W'] = np.array(cat2['W'], dtype=float)
        cat2['Weu'] = np.array(cat2['Weu'], dtype=float)
        cat2['Wel'] = np.array(cat2['Wel'], dtype=float)
        cat2['F'] = np.array(cat2['F'], dtype=float)
        cat2['Feu'] = np.array(cat2['Feu'], dtype=float)
        cat2['Fel'] = np.array(cat2['Fel'], dtype=float)
        cat2['DM'] = np.array(cat2['DM'], dtype=float)
        cat2['DM_NE2001'] = np.array(cat2['DM_NE2001'], dtype=float)
        cat2['DM_YMW16'] = np.array(cat2['DM_YMW16'], dtype=float)
        cat2['SURVEY'] = np.array(cat2['SURVEY'])
        cat2['Gain'] = np.array(cat2['Gain'], dtype=float)
        cat2['Tsys'] = np.array(cat2['Tsys'], dtype=float)
        cat2['BW'] = np.array(cat2['BW'], dtype=float)
        cat2['Npol'] = np.array(cat2['Npol'], dtype=float)
        cat2['SN0'] = np.array(cat2['SN0'], dtype=float)
        return cat2

    def LoadSvyInfo(self, fname):
        cat = np.loadtxt(fname, dtype=str)
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0,i]] = cat[1:, i]
        cat2['SURVEY'] = np.array(cat2['SURVEY'])
        cat2['FOV'] = np.array(cat2['FOV'], dtype=float)
        cat2['TIME'] = np.array(cat2['TIME'], dtype=float)
        cat2['Gain'] = np.array(cat2['Gain'], dtype=float)
        cat2['Tsys'] = np.array(cat2['Tsys'], dtype=float)
        cat2['BW'] = np.array(cat2['BW'], dtype=float)
        cat2['Npol'] = np.array(cat2['Npol'], dtype=float)
        cat2['SN0'] = np.array(cat2['SN0'], dtype=float)
        return cat2

    def LoadSimuData(self, fname):
        """读取模拟 FRB 数据文件（# 开头为 header 行，空格分隔）"""
        with open(fname) as f:
            header_line = f.readline().strip()
            if header_line.startswith('#'):
                header_line = header_line[1:]
            col_names = header_line.split()
        cat = np.loadtxt(fname, dtype=float, comments='#', skiprows=1)
        cat2 = {}
        for i, name in enumerate(col_names):
            cat2[name] = cat[:, i]
        return cat2

    def LoadFitsCatalog(self, fname):
        """读取 CHIME/FRB Catalog 2 FITS 文件

        返回: vF(fluence Jy·ms), vW(width ms), vz(redshift), vDM_obs, vDM_ne2001, vDM_ymw16, vSVY

        CHIME Catalog 2 关键列名:
          fluence, fluence_err — fluence 及误差 [Jy ms]
          bc_width — 校准后脉冲宽度 [ms]
          dm_fitb — 拟合 DM [pc cm^-3]
          dm_exc_ne2001 — DM - DM_MW(NE2001)，即河外 DM
          dm_exc_ymw16 — DM - DM_MW(YMW16)，即河外 DM
          红移：Catalog 2 中无红移列，需后续用 P(z|DM) 处理
        """
        try:
            from astropy.io import fits
        except ImportError:
            raise ImportError("需要安装 astropy: pip install astropy")

        with fits.open(fname) as hdul:
            data = hdul[1].data
            col_names = data.columns.names

            # === fluence ===
            fluence_key = next((k for k in ['fluence', 'Fluence'] if k in col_names), None)
            if fluence_key is None:
                raise KeyError(f"未找到 fluence 列，可用列: {col_names}")
            vF = np.array(data[fluence_key], dtype=float)

            # === width ===
            width_key = next((k for k in ['bc_width', 'width_fitb', 'width', 'Width'] if k in col_names), None)
            if width_key is None:
                raise KeyError(f"未找到 width 列，可用列: {col_names}")
            vW = np.array(data[width_key], dtype=float)
            # 优先读取 header 标记的单位，否则用中位数启发式判断
            hdr = hdul[1].header
            wunit = hdr.get('WUNIT', '').lower()
            if wunit == 'ms':
                pass  # 已经是毫秒，无需转换
            else:
                med_w = np.nanmedian(vW)
                if np.isfinite(med_w) and med_w < 1.0:
                    vW = vW * 1000.0

            # === DM_obs（观测 DM） ===
            dm_key = next((k for k in ['dm_obs', 'dm_fitb', 'bonsai_dm', 'dm', 'DM'] if k in col_names), None)
            if dm_key is None:
                raise KeyError(f"未找到 DM 列，可用列: {col_names}")
            vDM_obs = np.array(data[dm_key], dtype=float)

            # === DM_exc_ne2001: DM - DM_MW,NE2001（河外 DM） ===
            ne2001_key = next((k for k in ['dm_exc_ne2001', 'DM_NE2001'] if k in col_names), None)
            if ne2001_key is None:
                raise KeyError(f"FITS 缺少 NE2001 河外 DM 列，可用列: {col_names}")
            vDM_ne2001 = np.array(data[ne2001_key], dtype=float)

            # === DM_exc_ymw16: DM - DM_MW,YMW16（河外 DM） ===
            ymw16_key = next((k for k in ['dm_exc_ymw16', 'DM_YMW16'] if k in col_names), None)
            if ymw16_key is None:
                raise KeyError(f"FITS 缺少 YMW16 河外 DM 列，可用列: {col_names}")
            vDM_ymw16 = np.array(data[ymw16_key], dtype=float)

            # === 红移（Catalog 2 通常无此列） ===
            z_key = next((k for k in ['redshift', 'z', 'Z'] if k in col_names), None)
            vz_raw = np.array(data[z_key], dtype=float) if z_key else np.full(len(vDM_obs), np.nan)

            # === 巡天（Catalog 2 为单一 CHIME 巡天，此列可选） ===
            svy_key = next((k for k in ['survey', 'SURVEY', 'telescope'] if k in col_names), None)
            if svy_key is not None:
                vSVY = np.array(data[svy_key])
            else:
                vSVY = np.array(['CHIME'] * len(vDM_obs))

        return vF, vW, vDM_obs, vDM_ne2001, vDM_ymw16, vSVY

    def resolve_redshift(self, vDM, vDM_gal, vZ_catalog):
        """统一红移来源：优先用 catalog 红移，否则从 DM 推断

        Args:
            vDM: 观测 DM 数组
            vDM_gal: 银河系 DM 贡献数组
            vZ_catalog: catalog 中的红移（NaN 表示缺失）

        Returns:
            vz: 统一的红移数组
        """
        cos = Cosmology()
        vz = np.zeros(len(vDM))
        for i in range(len(vDM)):
            if not np.isnan(vZ_catalog[i]):
                vz[i] = vZ_catalog[i]
            else:
                dm_cosmic = vDM[i] - vDM_gal[i]
                # 用典型 DM_host ≈ 50 pc cm⁻³ 估算（对应 μ_host=4.44, log-normal 中位数）
                dm_host_typical = 50.0
                if dm_cosmic > dm_host_typical:
                    vz[i] = cos.GetZ(dm_cosmic, dm_host_typical)
                else:
                    vz[i] = 0.01  # 默认最小红移
        return vz

def gammainc(alpha, x):
    if alpha==0:
        return -spf.expi(-x)

    elif (alpha<0):
        return (gammainc(alpha+1,x)-np.power(x, alpha)*np.exp(-x))/alpha

    else:
        return spf.gammaincc(alpha,x)*spf.gamma(alpha)

def getargv(argv, key):
    for i in range(0, len(argv)):
        arg = argv[i]
        if (arg == key):
            return argv[i + 1]

def chkargv(argv, key):
    for i in range(0, len(argv)):
        arg = argv[i]
        if (arg == key):
            return True
    return False

def Sampling1D(x, y, x1, x2, n):
    ymax = np.max(y)
    if ymax <= 0:
        raise ValueError("Sampling1D: 目标分布在采样范围内全为零，无法采样")
    nt = 0
    res = np.array([])
    while (nt<n):
        fuc = interp1d(x, y / ymax, bounds_error=False, fill_value=0.0)
        vx = np.random.uniform(x1, x2, n-nt)
        vy = np.random.uniform(0, 1, n-nt)
        res = np.append(res, vx[vy <= fuc(vx)])
        nt = len(res)
    return res

def SamplingND(fuc, par_range, maxv_ori, n):
    nt = 0
    maxv = maxv_ori
    res = np.array([])
    npar, m = par_range.shape
    res = res.reshape((0, npar))
    while (nt < n):
        vpar = np.random.uniform(0, 1, (n-nt, npar))
        for i in range(npar):
            lv = par_range[i,0]
            rv = par_range[i,1]
            vpar[:,i] = vpar[:,i] * (rv-lv) + lv
        
        vy = np.random.uniform(0, 1, n-nt)*maxv
        fv = fuc(vpar)
        if (np.max(fv) > maxv):
            maxv = np.max(fv)*1.5
            nt = 0
            res = np.array([])
            res = res.reshape((0, npar))
        else:
            res = np.vstack( (res, vpar[vy <= fv,:]))
            nt, m = res.shape
    return res
