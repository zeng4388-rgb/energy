"""
frb_util.py - FRB Isotropic Equivalent Energy (E_iso) Luminosity Function Utilities
Converted to Python 3 from RuiLuoAstro/frblf_erd (Luo et al. 2018, 2020)
Modified to use CHIME/FRB Catalog 1, Catalog 2, Repeater and Baseband data
Calculates FRB isotropic equivalent energy E_iso instead of luminosity
"""

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import scipy.special as spf
import json
import os


class Cosmology:
    def __init__(self, omegam=0.308, omegal=0.692, omegab=0.0484):
        self.Omega_m = omegam
        self.Omega_b = omegab
        self.Omega_L = omegal
        self.c = 2.9979245800e10          # cm/s
        self.pc2cm = 3.08567758e18
        self.km2cm = 1e5
        self.Mpc2cm = 3.08567758e24
        self.Gpc2cm = 3.08567758e27
        self.Jy2CGS = 1e-23
        self.MHz2Hz = 1e6
        self.Jyms2CGS = 1e-26
        self.h0 = 0.6781
        self.H0 = self.h0 * 100 * self.km2cm / self.Mpc2cm
        self.Rhoc = 1.88 * self.h0 * self.h0 * 1e-29    # g/cm^3
        self.Nc = self.Rhoc / 1.6726e-24                  # 1/cm^3
        self.f_IGM = 0.83
        self.z0 = 0.8

        # Pre-compute interpolation tables
        self.vz = np.arange(-7, 6, 0.03)
        self.vz = np.power(10., self.vz)
        self.vz[0] = 0

        func = lambda z: self.c / self.H(z)
        self.vd = [integrate.quad(func, 0, zv)[0] for zv in self.vz]
        self.cd_interp = interp1d(self.vz, self.vd)

        func2 = lambda z: (1 + z) * self.c / self.H(z) * self.Omega_b * self.Nc
        self.vdm = [integrate.quad(func2, 0, zv)[0] for zv in self.vz]
        self.dmigm = interp1d(self.vz, self.vdm)

        self.vld = self.Luminosity_Distance(self.vz)
        self.Ld2z = interp1d(self.vld, self.vz)
        self.Ld1 = self.Luminosity_Distance(1.0)
        self.Cd1 = self.Comoving_Distance(1.0)

    def E(self, z):
        """Logarithmic time derivative of scale factor"""
        return np.sqrt(self.Omega_m * np.power(1 + z, 3.0) + self.Omega_L)

    def H(self, z):
        """Hubble ratio"""
        return self.H0 * self.E(z)

    def Comoving_Distance(self, z):
        """Comoving distance in cm"""
        return self.cd_interp(z)

    def dVdOdz(self, z):
        """Differential comoving volume dV/dz/dOmega in Gpc^3"""
        drdz = self.c / self.H(z)
        cd2 = self.Comoving_Distance(z) ** 2
        dcv = cd2 * drdz / (self.Gpc2cm ** 3)
        return dcv

    def Luminosity_Distance(self, z):
        """Luminosity distance in Mpc"""
        dl = (1 + z) * self.Comoving_Distance(z)
        return dl

    def Luminosity(self, z, f=1., dnu=1000.):
        """Intrinsic luminosity from flux (Jy), spectral width dnu (MHz)"""
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        lum = f * self.Jy2CGS * dnu * self.MHz2Hz * 4 * np.pi * ld2
        return lum

    def Energy(self, z, flu=1.0, dnu=1000.):
        """
        Isotropic equivalent energy E_iso from fluence
        flu: fluence in Jy ms
        dnu: intrinsic spectral width in MHz (default 1000)
        z: redshift
        Returns E_iso in erg
        """
        ld = self.Luminosity_Distance(z)  # in Mpc
        ld2 = ld * ld
        # E_iso = fluence / (1+z) * 4*pi*DL^2 * dnu
        # flu in Jyms = 1e-26 erg/cm^2/Hz * 1e-3 s
        ener = flu / (1 + z) * self.Jyms2CGS * dnu * self.MHz2Hz * 4 * np.pi * ld2
        return ener

    def Eiso_to_Fluence(self, z, eiso, dnu=1000.):
        """Convert E_iso (erg) to observed fluence (Jy ms)"""
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        flu = eiso * (1 + z) / (4 * np.pi * ld2 * dnu * self.MHz2Hz * self.Jyms2CGS)
        return flu

    def Luminosity_to_Flux(self, z, lum, dnu=1000):
        """Flux density from luminosity"""
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        flux = lum / 4 / np.pi / ld2 / dnu / self.MHz2Hz / self.Jy2CGS
        return flux

    def Energy_to_Flu(self, z, ener, dnu=1000):
        """Observed fluence from intrinsic energy"""
        ld = self.Luminosity_Distance(z)
        ld2 = ld * ld
        flu = ener * (1 + z) / 4 / np.pi / ld2 / dnu / self.MHz2Hz / self.Jyms2CGS
        return flu

    def DispersionMeasure_IGM(self, z, chi=7. / 8):
        """Dispersion measure of IGM"""
        return self.dmigm(z) * self.f_IGM * chi / self.Mpc2cm * 1e6

    def DMeq(self, z, dme, dmhost):
        """DM equation solved to get redshift"""
        dmi = self.DispersionMeasure_IGM(z)
        dmh = dmhost / (1 + z)
        return dmi + dmh - dme

    def GetZ(self, dme, dmhost):
        """Solve DM equation for redshift"""
        z = fsolve(self.DMeq, self.z0, args=(dme, dmhost))
        z = float(z)
        return z


class Telescope:
    def __init__(self):
        self.MHz2Hz = 1e6
        self.ms2s = 1e-3

    def RMEq(self, snr, g, tsys, npol, bw, w):
        """Radiometer equation: sensitivity in Jy ms"""
        s = snr * tsys / g / np.sqrt(npol * bw * self.MHz2Hz * w * self.ms2s)
        return s


class AstroDistribution:
    def __init__(self):
        self.tel = Telescope()
        self.cos = Cosmology()
        # DM distribution parameters for different galaxy types
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
        Schechter luminosity/energy function per logarithmic unit
        logl: logarithmic luminosity or energy
        phis: normalization
        logls: cut-off luminosity/energy
        alpha: power index
        """
        l = np.power(10., logl)
        ls = np.power(10., logls)
        phi = np.log(10) * phis * np.power(l / ls, (alpha + 1)) * np.exp(-l / ls)
        return phi

    def log_Schechter_log(self, logl, alpha, logls, logl0):
        """Logarithmic Schechter function"""
        phi = (logl - logls) * (alpha + 1) * np.log(10.) - np.power(10., logl - logls)
        lik = phi.copy()
        lik[logl < logl0] = -1e99
        return lik

    def log_IntBeam(self, logl, alpha, logls, logl0):
        """Beam efficiency integration - uses custom gammainc for negative alpha support"""
        ratio = np.power(10., logl - logls)
        lik0 = gammainc(alpha + 1, ratio) - gammainc(alpha + 1, 2 * ratio)
        lik = lik0 / np.log(2)
        ind = lik <= 0
        lik[ind] = 1e-199
        loglik = np.log(lik)
        loglik[logl < logl0] = -1e99
        return loglik

    def IntLum(self, eps, alpha, logls, logl0):
        """Integrated luminosity/energy fraction - uses custom gammainc for negative alpha"""
        with np.errstate(invalid='ignore'):
            ratio = np.power(10., logl0 - logls)
        return gammainc(alpha + 1, ratio)

    def Distribution_Local_galaxy_DM(self, dmv, vpar):
        """Double Gaussian DM distribution for host galaxies"""
        val = vpar[0] * np.exp(-np.power((dmv - vpar[1]) / vpar[2], 2.)) \
              + vpar[3] * np.exp(-np.power((dmv - vpar[4]) / vpar[5], 2.))
        return val

    def ThetaFunc(self, x):
        """Heaviside step function"""
        return 0.5 * (np.sign(x) + 1)

    def func_gaussian(self, dmv, vpar):
        """Default Gaussian DM distribution"""
        dmoff = dmv - vpar[0]
        sig = vpar[1] * vpar[1]
        return np.exp(-0.5 * dmoff * dmoff / sig) * self.ThetaFunc(dmv)

    def func_uniform(self, dmv, vpar):
        """Uniform DM distribution"""
        pdf = np.ones(dmv.shape)
        pdf[dmv >= vpar[1]] = 0
        pdf[dmv <= vpar[0]] = 0
        return pdf

    def Distribution_HostGalaxyDM(self, dmv0, fgalaxy_type=None, vpar=np.array([0, 50])):
        """Host galaxy DM distribution for different galaxy types"""
        if not fgalaxy_type:
            fgalaxy_type = self.func_gaussian
        elif fgalaxy_type == 'ETG':
            dmv = dmv0.copy()
            dmv[dmv0 < 1e-9] = np.ones(dmv[dmv0 < 1e-9].shape) * 1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_etg)
            res[dmv0 <= 0] = np.zeros(dmv0[dmv0 <= 0].shape)
            return res
        elif fgalaxy_type == 'LTG_NE2001':
            dmv = dmv0.copy()
            dmv[dmv0 < 1e-9] = np.ones(dmv[dmv0 < 1e-9].shape) * 1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_ltg_ne2001)
            res[dmv0 <= 0] = np.zeros(dmv0[dmv0 <= 0].shape)
            return res
        elif fgalaxy_type == 'LTG_YMW16':
            dmv = dmv0.copy()
            dmv[dmv0 < 1e-9] = np.ones(dmv[dmv0 < 1e-9].shape) * 1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_ltg_ymw16)
            res[dmv0 <= 0] = np.zeros(dmv0[dmv0 <= 0].shape)
            return res
        elif fgalaxy_type == 'ALG_NE2001':
            dmv = dmv0.copy()
            dmv[dmv0 < 1e-9] = np.ones(dmv[dmv0 < 1e-9].shape) * 1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_alg_ne2001)
            res[dmv0 <= 0] = np.zeros(dmv0[dmv0 <= 0].shape)
            return res
        elif fgalaxy_type == 'ALG_YMW16':
            dmv = dmv0.copy()
            dmv[dmv0 < 1e-9] = np.ones(dmv[dmv0 < 1e-9].shape) * 1e-9
            res = self.Distribution_Local_galaxy_DM(np.log10(dmv), self.vpar_alg_ymw16)
            res[dmv0 <= 0] = np.zeros(dmv0[dmv0 <= 0].shape)
            return res
        else:
            return fgalaxy_type(dmv0, vpar)

    def log_Distribution_HostGalaxyDM(self, dmv, fgalaxy_type=None, vpar=np.array([0, 50])):
        """Log of host galaxy DM distribution"""
        if not fgalaxy_type:
            dmoff = dmv[dmv > 0] - vpar[0]
            sig = vpar[1] * vpar[1]
            res = dmv.copy()
            res[dmv < 0] = -1e99
            res[dmv > 0] = -0.5 * dmoff * dmoff / sig
            return res
        else:
            val = self.Distribution_HostGalaxyDM(dmv, fgalaxy_type=fgalaxy_type)
            ind = val <= 0
            indv = val > 0
            val[indv] = np.log(val[indv])
            val[ind] = np.ones(val[ind].shape) * -1e99
            return val

    def SFR(self, z):
        """Star formation rate history (Hopkins & Beacom 2006)"""
        sfr = (0.017 + 0.13 * z) / (1 + np.power(z / 3.3, 5.3))
        return sfr

    def kappa(self, z):
        """Normalized SFH"""
        return np.sqrt(self.SFR(0) / self.SFR(z))

    def Distribution_volume(self, z):
        """Differential comoving volume"""
        r = self.cos.Comoving_Distance(z) / self.cos.Comoving_Distance(1)
        pv = r * r / self.cos.E(z)
        return pv

    def log_Distribution_volume(self, z):
        """Log differential comoving volume"""
        if type(z) == np.ndarray:
            ind = z < 0
            pv = np.log(self.Distribution_volume(z))
            pv[ind] = -1e99
            return pv
        else:
            if z < 0:
                return -1e99
            else:
                return np.log(self.Distribution_volume(z))

    def IntDMsrc(self, u1, u2, vpar):
        """Analytic integral for DMsrc marginalization"""
        a1, b1, c1 = vpar[0], vpar[1], vpar[2]
        a2, b2, c2 = vpar[3], vpar[4], vpar[5]
        k1 = np.power(10., b1) * a1 * c1 * np.exp(c1 * c1 * np.log(10) * np.log(10) / 4)
        k2 = np.power(10., b2) * a2 * c2 * np.exp(c2 * c2 * np.log(10) * np.log(10) / 4)
        q1 = (c1 * c1 * np.log(10) * np.log(10) + b1 * np.log(100) - 2 * np.log(u1)) / c1 / np.log(100)
        q2 = (c1 * c1 * np.log(10) * np.log(10) + b1 * np.log(100) - 2 * np.log(u2)) / c1 / np.log(100)
        q3 = (c2 * c2 * np.log(10) * np.log(10) + b2 * np.log(100) - 2 * np.log(u1)) / c2 / np.log(100)
        q4 = (c2 * c2 * np.log(10) * np.log(10) + b2 * np.log(100) - 2 * np.log(u2)) / c2 / np.log(100)
        int_h = np.log(10) * np.sqrt(np.pi) / 2 * (
            k1 * (spf.erf(q2) - spf.erf(q1)) + k2 * (spf.erf(q4) - spf.erf(q3)))
        int_hs = int_h / self.DMsmax
        return int_hs

    def log_IntDMsrc(self, u1, u2, gtype=None):
        """Log of DMsrc marginalization integral"""
        if not gtype:
            return -1e99
        elif gtype == 'ETG':
            return np.log(self.IntDMsrc(u1, u2, self.vpar_etg))
        elif gtype == 'LTG_NE2001':
            return np.log(self.IntDMsrc(u1, u2, self.vpar_ltg_ne2001))
        elif gtype == 'LTG_YMW16':
            return np.log(self.IntDMsrc(u1, u2, self.vpar_ltg_ymw16))
        elif gtype == 'ALG_NE2001':
            return np.log(self.IntDMsrc(u1, u2, self.vpar_alg_ne2001))
        elif gtype == 'ALG_YMW16':
            return np.log(self.IntDMsrc(u1, u2, self.vpar_alg_ymw16))
        else:
            return -1e99

    def dis_logw(self, logw0, mu, sigma):
        """Log-normal pulse width distribution"""
        a = 1. / np.sqrt(2. * np.pi * sigma * sigma)
        b = np.exp(-(logw0 - mu) * (logw0 - mu) / 2 / sigma / sigma)
        return a * b

    def log_dis_logw(self, logw0, mu, sigma):
        """Log of log-normal pulse width distribution"""
        a = 2 * np.pi * sigma * sigma
        b = -(logw0 - mu) * (logw0 - mu) / 2 / sigma / sigma
        return b - 0.5 * np.log(a)

    def log_distr_fdmwz(self, dnu, logfluence, dme, logw, z, alpha, logeiso0, logeisomax, mu, sigma, gtype=None):
        """
        Log joint distribution of fluence, DM and redshift (E_iso version)
        dnu: intrinsic bandwidth (MHz)
        logfluence: log10 fluence (Jy ms)
        dme: excess DM
        logw: log10 observed pulse width (ms)
        z: redshift
        alpha: power-law index of E_iso function
        logeiso0: log10 minimum E_iso (erg)
        logeisomax: log10 characteristic E_iso (erg)
        mu, sigma: pulse width distribution params
        gtype: host galaxy DM type
        """
        fluence = np.power(10., logfluence)
        # Calculate E_iso from fluence
        logeiso = np.log10(self.cos.Energy(z, flu=fluence, dnu=dnu))
        logint1 = self.log_IntBeam(logeiso, alpha, logeisomax, logeiso0)
        logw0 = logw - np.log10(1 + z)
        logfw = self.log_dis_logw(logw0, mu, sigma)
        logfz = self.log_Distribution_volume(z)
        dmi = self.cos.DispersionMeasure_IGM(z)
        u1 = (dme - dmi) * (1 + z) * self.kappa(z)
        u2 = ((dme - dmi) * (1 + z) - self.DMsmax) * self.kappa(z)
        logint2 = np.ones(u2.shape) * (-1e99)
        ind = u2 > 0
        logint2[ind] = self.log_IntDMsrc(u1[ind], u2[ind], gtype=gtype)
        loglikv = logint1 + logfz + logfw + logint2 + np.log(1 + z)
        return loglikv

    def log_distr_fdmw(self, dnu, logfluence, dme, logw, alpha, logeiso0, logeisomax, mu, sigma, gtype=None):
        """
        Log joint distribution of fluence and DM after marginalizing redshift
        (E_iso version)
        """
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 1000
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        lik = 0
        for z in vz:
            likv = np.exp(self.log_distr_fdmwz(
                dnu, logfluence, dme, logw, z, alpha, logeiso0, logeisomax, mu, sigma, gtype=gtype))
            lik += z * stepz * likv
        ind = lik > 0
        ind2 = lik <= 0
        loglik = lik.copy()
        loglik[ind] = np.log(lik[ind])
        loglik[ind2] = np.ones(loglik[ind2].shape) * -1e99
        return loglik

    def Norm1D(self, sn0, bw, npol, g, tsys, dnu, alpha, logeiso0, logeisomax, mu, sigma):
        """
        Normalization factor for dimensionless likelihood (E_iso version)
        """
        stepz = (np.log(self.Zmax) - np.log(self.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.Zmin), np.log(self.Zmax), stepz))
        stepeps = (1 - 0.5) / 200.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.Wmax) - np.log10(self.Wmin)) / 100.
        vlogw = np.arange(np.log10(self.Wmin), np.log10(self.Wmax), steplogw)
        nf = 0
        for z in vz:
            vw = np.power(10, vlogw) * (1 + z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            # Convert fluence threshold to energy at this redshift
            et = self.cos.Energy(z, flu=ft, dnu=dnu)
            loget = np.log10(et)
            ind = loget < logeiso0
            loget[ind] = logeiso0
            int_eps = np.zeros(loget.shape)
            for i in np.arange(len(loget)):
                int_eps[i] = np.sum(
                    self.IntLum(veps, alpha, logeisomax, loget[i] * np.ones(veps.shape))
                    / veps / np.log(2) * stepeps)
            int_w = np.sum(int_eps * self.dis_logw(vlogw, mu, sigma) * steplogw)
            fz = self.Distribution_volume(z)
            nf += z * stepz * fz * int_w
        if nf <= 0:
            nf = 1e-199
        return nf


class EventRate:
    def __init__(self):
        self.cos = Cosmology()
        self.ad = AstroDistribution()
        self.tel = Telescope()
        self.s2h = 1. / 3600
        self.yr2hr = 365 * 24.
        self.rad2deg2 = 3282.806350011744
        self.Gpc2Mpc = 1e3

    def rate_2d(self, sn0, bw, npol, g, tsys, dnu, phis, alpha, logeiso0, logeisomax, mu, sigma):
        """Event rate density (E_iso version)"""
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        stepeps = (1 - 0.5) / 200.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 100.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        rho = 0
        for z in vz:
            vw = np.power(10, vlogw) * (1 + z)
            ft = self.tel.RMEq(sn0, g, tsys, npol, bw, vw)
            # Convert fluence threshold to energy
            et = self.cos.Energy(z, flu=ft, dnu=dnu)
            loget = np.log10(et)
            ind = loget < logeiso0
            loget[ind] = logeiso0
            int_eps = np.zeros(loget.shape)
            for i in range(len(loget)):
                int_eps[i] = np.sum(
                    phis * self.ad.IntLum(veps, alpha, logeisomax, loget[i] * np.ones(veps.shape))
                    / veps / np.log(2) * stepeps)
            int_w = np.sum(int_eps * self.ad.dis_logw(vlogw, mu, sigma) * steplogw)
            fz = self.cos.dVdOdz(z) / (1 + z)
            rho += z * stepz * fz * int_w
        rho_deg = rho / self.rad2deg2 / self.yr2hr
        return rho_deg

    def log_dis_poi(self, rho, N, Omega, T):
        """Poisson likelihood for event count"""
        lamda = rho * Omega * T
        ind = lamda <= 0
        lamda[ind] = 1e-199
        loglik = N * np.log(lamda) - lamda - spf.gammaln(N)
        return loglik

    def Rate(self, logft, dnu, phis, alpha, logeiso0, logeisomax, mu, sigma):
        """Event rate for given fluence threshold (E_iso version)"""
        stepz = (np.log(self.ad.Zmax) - np.log(self.ad.Zmin)) / 1000.
        vz = np.exp(np.arange(np.log(self.ad.Zmin), np.log(self.ad.Zmax), stepz))
        stepeps = (1 - 0.5) / 200.
        veps = np.arange(0.5, 1, stepeps)
        steplogw = (np.log10(self.ad.Wmax) - np.log10(self.ad.Wmin)) / 100.
        vlogw = np.arange(np.log10(self.ad.Wmin), np.log10(self.ad.Wmax), steplogw)
        rho = 0
        for z in vz:
            vw = np.power(10, vlogw) * (1 + z)
            ft = np.power(10, logft)
            et = self.cos.Energy(z, flu=ft, dnu=dnu)
            loget = np.log10(et)
            if loget < logeiso0:
                loget = logeiso0
            int_eps = np.sum(
                phis * self.ad.IntLum(veps, alpha, logeisomax, loget * np.ones(veps.shape))
                / veps / np.log(2) * stepeps)
            int_w = np.sum(int_eps * self.ad.dis_logw(vlogw, mu, sigma) * steplogw)
            fz = self.cos.dVdOdz(z) / (1 + z)
            rho += z * stepz * fz * int_w
        rho_deg = rho / self.rad2deg2 / self.yr2hr
        return rho_deg

    def Rfrb(self, phis, alpha, logeiso0, logeisomax):
        """Cumulative FRB rate above E_iso_min"""
        ratio = np.power(10., logeiso0 - logeisomax)
        return phis * spf.gammaincc(alpha + 1, ratio) * spf.gamma(alpha + 1)


class Loadfiles:
    """Load CHIME/FRB catalog data"""

    def LoadCatalogue(self, fname):
        """Load original format catalog (space-separated txt)"""
        cat = np.loadtxt(fname, dtype='str')
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0, i]] = cat[1:, i]
        for key in ['S', 'Seu', 'Sel', 'W', 'Weu', 'Wel', 'F', 'Feu', 'Fel',
                     'DM', 'DM_NE2001', 'DM_YMW16', 'Gain', 'Tsys', 'BW', 'Npol', 'SN0']:
            cat2[key] = np.array(cat2[key], dtype=float)
        cat2['SURVEY'] = np.array(cat2['SURVEY'])
        return cat2

    def LoadSvyInfo(self, fname):
        """Load survey info file"""
        cat = np.loadtxt(fname, dtype='str')
        row, col = cat.shape
        cat2 = {}
        for i in range(col):
            cat2[cat[0, i]] = cat[1:, i]
        for key in ['FOV', 'TIME', 'Gain', 'Tsys', 'BW', 'Npol', 'SN0']:
            cat2[key] = np.array(cat2[key], dtype=float)
        cat2['SURVEY'] = np.array(cat2['SURVEY'])
        return cat2

    def LoadCHIMECat1(self, fname='chime_cat1.json'):
        """
        Load CHIME/FRB Catalog 1 (JSON format)
        Returns dict with keys: fluence, fluence_err, dm, dm_exc_ne2001, dm_exc_ymw16,
        width, snr, flux, flux_err, ra, dec, repeater_name, tns_name, etc.
        """
        with open(fname, 'r') as f:
            data = json.load(f)

        cat = {'tns_name': [], 'repeater_name': [], 'fluence': [], 'fluence_err': [],
               'dm': [], 'dm_exc_ne2001': [], 'dm_exc_ymw16': [],
               'width': [], 'snr': [], 'flux': [], 'flux_err': [],
               'ra': [], 'dec': [], 'peak_freq': [], 'low_freq': [], 'high_freq': []}

        def safe_float(val, default=0.0):
            """Safely convert value to float, handling '<' prefixed strings"""
            if val is None or val == '' or val == -9999:
                return default
            if isinstance(val, str):
                val = val.strip().lstrip('<').lstrip('>')
                try:
                    return float(val)
                except ValueError:
                    return default
            try:
                if np.isnan(float(val)):
                    return default
                return float(val)
            except (ValueError, TypeError):
                return default

        for d in data:
            flu = d.get('fluence', -9999)
            flu = safe_float(flu, -9999)
            if flu <= 0:
                continue

            cat['tns_name'].append(d.get('tns_name', ''))
            rn = d.get('repeater_name', '')
            cat['repeater_name'].append(rn if rn not in ['', '-9999', None] else '')
            cat['fluence'].append(flu)
            cat['fluence_err'].append(safe_float(d.get('fluence_err', 0)))
            cat['dm'].append(safe_float(d.get('dm_fitb', 0)))
            cat['dm_exc_ne2001'].append(safe_float(d.get('dm_exc_ne2001', 0)))
            cat['dm_exc_ymw16'].append(safe_float(d.get('dm_exc_ymw16', 0)))
            cat['width'].append(safe_float(d.get('width_fitb', 0)))
            cat['snr'].append(safe_float(d.get('snr_fitb', 0)))
            cat['flux'].append(safe_float(d.get('flux', 0)))
            cat['flux_err'].append(safe_float(d.get('flux_err', 0)))
            cat['ra'].append(safe_float(d.get('ra', 0)))
            cat['dec'].append(safe_float(d.get('dec', 0)))
            cat['peak_freq'].append(safe_float(d.get('peak_freq', 600), 600))
            cat['low_freq'].append(safe_float(d.get('low_freq', 400), 400))
            cat['high_freq'].append(safe_float(d.get('high_freq', 800), 800))

        for key in cat:
            if key in ['tns_name', 'repeater_name']:
                cat[key] = np.array(cat[key])
            else:
                cat[key] = np.array(cat[key], dtype=float)

        # Separate repeaters and one-offs
        is_repeater = np.array([r != '' for r in cat['repeater_name']])
        cat['is_repeater'] = is_repeater
        return cat

    def LoadCHIMECat2(self, fname='chime_cat2.json'):
        """Load CHIME/FRB Catalog 2 (JSON format)"""
        with open(fname, 'r') as f:
            data = json.load(f)

        def safe_float(val, default=0.0):
            """Safely convert value to float, handling NaN and special strings"""
            if val is None or val == '' or val == -9999:
                return default
            if isinstance(val, str):
                val = val.strip().lstrip('<').lstrip('>')
                try:
                    return float(val)
                except ValueError:
                    return default
            try:
                if isinstance(val, float) and np.isnan(val):
                    return default
                return float(val)
            except (ValueError, TypeError):
                return default

        cat = {'tns_name': [], 'repeater_name': [], 'event_id': [],
               'fluence': [], 'fluence_err': [],
               'dm': [], 'dm_exc_ne2001': [], 'dm_exc_ymw16': [],
               'width': [], 'snr': [], 'flux': [], 'flux_err': [],
               'ra': [], 'dec': [], 'peak_freq': [], 'low_freq': [], 'high_freq': [],
               'excluded_flag': [], 'sidelobe_flag': []}

        for d in data:
            flu = safe_float(d.get('fluence', None), -9999)
            if flu <= 0:
                continue

            cat['tns_name'].append(d.get('tns_name', ''))
            rn = d.get('repeater_name', '')
            cat['repeater_name'].append(rn if rn not in ['', '-9999', None] else '')
            cat['event_id'].append(str(d.get('event_id', '')))
            cat['fluence'].append(flu)
            cat['fluence_err'].append(safe_float(d.get('fluence_err', 0)))
            cat['dm'].append(safe_float(d.get('dm_fitb', 0)))
            cat['dm_exc_ne2001'].append(safe_float(d.get('dm_exc_ne2001', 0)))
            cat['dm_exc_ymw16'].append(safe_float(d.get('dm_exc_ymw16', 0)))
            cat['width'].append(safe_float(d.get('width_fitb', 0)))
            cat['snr'].append(safe_float(d.get('snr_fitb', 0)))
            cat['flux'].append(safe_float(d.get('flux', 0)))
            cat['flux_err'].append(safe_float(d.get('flux_err', 0)))
            cat['ra'].append(safe_float(d.get('ra', 0)))
            cat['dec'].append(safe_float(d.get('dec', 0)))
            cat['peak_freq'].append(safe_float(d.get('peak_freq', 600), 600))
            cat['low_freq'].append(safe_float(d.get('low_freq', 400), 400))
            cat['high_freq'].append(safe_float(d.get('high_freq', 800), 800))
            cat['excluded_flag'].append(int(safe_float(d.get('excluded_flag', 0))))
            cat['sidelobe_flag'].append(int(safe_float(d.get('sidelobe_flag', 0))))

        for key in cat:
            if key in ['tns_name', 'repeater_name', 'event_id']:
                cat[key] = np.array(cat[key])
            else:
                cat[key] = np.array(cat[key], dtype=float)

        is_repeater = np.array([r != '' for r in cat['repeater_name']])
        cat['is_repeater'] = is_repeater
        return cat

    def LoadSimuData(self, fname):
        """Load simulated FRB data"""
        cat = np.loadtxt(fname, dtype='str', comments='XXX')
        row, col = cat.shape
        cat2 = {}
        cat[0, 0] = cat[0, 0][1:]  # Remove leading '#'
        for i in range(col):
            cat2[cat[0, i]] = cat[1:, i]
        for key in ['S', 'W', 'T', 'DMe', 'thres', 'logL', 'Z', 'DMi', 'DMh', 'DMs']:
            cat2[key] = np.array(cat2[key], dtype=float)
        return cat2


def gammainc(alpha, x):
    """Regularized lower incomplete gamma function"""
    if alpha == 0:
        return -spf.expi(-x)
    elif alpha < 0:
        return (gammainc(alpha + 1, x) - np.power(x, alpha) * np.exp(-x)) / alpha
    else:
        return spf.gammaincc(alpha, x) * spf.gamma(alpha)


def getargv(argv, key):
    """Get argument value from command line"""
    for i in range(0, len(argv)):
        if argv[i] == key:
            return argv[i + 1]
    return None


def chkargv(argv, key):
    """Check if argument exists in command line"""
    return key in argv


def Sampling1D(x, y, x1, x2, n):
    """1D rejection sampling"""
    nt = 0
    res = np.array([])
    while nt < n:
        fuc = interp1d(x, y / np.max(y))
        vx = np.random.uniform(x1, x2, n - nt)
        vy = np.random.uniform(0, 1, n - nt)
        res = np.append(res, vx[vy <= fuc(vx)])
        nt = len(res)
    return res


def SamplingND(fuc, par_range, maxv_ori, n):
    """N-dimensional rejection sampling"""
    nt = 0
    maxv = maxv_ori
    res = np.array([])
    npar, m = par_range.shape
    res = res.reshape((0, npar))
    while nt < n:
        vpar = np.random.uniform(0, 1, (n - nt, npar))
        for i in range(npar):
            lv = par_range[i, 0]
            rv = par_range[i, 1]
            vpar[:, i] = vpar[:, i] * (rv - lv) + lv
        vy = np.random.uniform(0, 1, n - nt) * maxv
        fv = fuc(vpar)
        if np.max(fv) > maxv:
            maxv = np.max(fv) * 1.5
            nt = 0
            res = np.array([])
            res = res.reshape((0, npar))
        else:
            res = np.vstack((res, vpar[vy <= fv, :]))
            nt, m = res.shape
        print(nt)
    return res
