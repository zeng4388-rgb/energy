#!/usr/bin/env python3
"""
pltpost.py - Plot posterior distributions from MultiNest output (E_iso version)
Supports both real CHIME data and simulated data

Usage:
    python pltpost.py -f nest_out/samp/cat1_ALG_YMW16 -o plots/cat1_posteriors.eps -title "CHIME Cat1" -up 1 -bo 1
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import argparse
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pymultinest

plt.style.use("classic")
mpl.rcParams['font.size'] = 24
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['axes.labelsize'] = 22


def plot2dposterier_withconf(dat, likli, indx=[], labels=[], rate=0.5,
                             par=np.array([]), parm=np.array([]),
                             rangedat=np.array([]), levels=[0.68]):
    """Plot 2D posteriors with confidence contours"""
    row, col = dat.shape
    if len(indx) == 0:
        indx = np.arange(0, col)
    if len(labels) == 0:
        labels = ['a'] * col

    dat = dat.copy()
    dat = dat[np.ix_(np.arange(int(row - row * rate), row), indx)]
    likli = likli[np.ix_(np.arange(int(row - row * rate), row))]
    row, col = dat.shape
    max_yticks = 4

    if len(rangedat) == 0:
        rangedat = np.zeros((col, 2))
        for i in range(0, col):
            rangedat[i, 0] = np.min(dat[:, i])
            rangedat[i, 1] = np.max(dat[:, i])

    npar = col

    for vari in range(0, npar):
        # Histogram
        ax = plt.subplot(npar, npar, vari * npar + vari + 1)
        ind = ((dat[:, vari] < rangedat[vari, 1]) & (dat[:, vari] > rangedat[vari, 0]))
        n, bins, patches = plt.hist(dat[ind, vari], 100, density=True,
                                     histtype='stepfilled',
                                     range=(rangedat[vari, 0], rangedat[vari, 1]))
        plt.setp(patches, 'facecolor', 'lightblue', 'alpha', 0.6)
        yloc = plt.MaxNLocator(max_yticks)
        ax.xaxis.set_major_locator(yloc)
        xloc = plt.MaxNLocator(max_yticks)
        ax.yaxis.set_major_locator(xloc)
        plt.xlabel(labels[vari])

        # Confidence intervals
        for ls in [0.68, 0.95]:
            histc, bin_edges = np.histogram(dat[:, vari], density=True, bins=100)
            vom = (bin_edges[:-1] + bin_edges[1:]) * 0.5
            ind_sort = np.argsort(histc)[::-1]
            v = histc.copy()
            v[0] = histc[ind_sort[0]]
            for i in range(1, len(ind_sort)):
                v[i] = v[i - 1] + histc[ind_sort[i]]
            v = v / float(np.sum(histc))
            thre_idx = np.where((v[:-1] < ls) & (v[1:] >= ls))[0]
            if len(thre_idx) > 0:
                thre = histc[ind_sort[thre_idx[0]]]
                lv = np.min(vom[histc >= thre])
                rv = np.max(vom[histc >= thre])
                print(f'{vari}-th parameter, sigma={ls}, lv={lv:.4f}, rv={rv:.4f}')
                if ls == 0.68:
                    plt.plot([lv, lv], [0, max(n)], ls='dashed', color='k', linewidth=1)
                    plt.plot([rv, rv], [0, max(n)], ls='dashed', color='k', linewidth=1)
                else:
                    plt.plot([lv, lv], [0, max(n)], ls='dotted', color='k', linewidth=1)
                    plt.plot([rv, rv], [0, max(n)], ls='dotted', color='k', linewidth=1)

        if par.size > 0:
            plt.plot([par[vari], par[vari]], [0, max(n)], ls='solid',
                     color='k', linewidth=2)

        if len(parm) > 0:
            for i in range(len(parm)):
                parm0 = parm[i]
                plt.plot([parm0[vari], parm0[vari]], [0, max(n)], ls='dashed',
                         color='k', linewidth=2)

        plt.xlim(rangedat[vari, :])

        # 2D contours
        for varj in range(vari + 1, npar):
            ax = plt.subplot(npar, npar, (vari) * npar + varj + 1)
            x = dat[:, varj]
            y = dat[:, vari]
            ind = ((x < rangedat[varj, 1]) & (x > rangedat[varj, 0]) &
                   (y < rangedat[vari, 1]) & (y > rangedat[vari, 0]))
            x = dat[ind, varj]
            y = dat[ind, vari]
            ngridx = 20
            ngridy = 30
            H, xedges, yedges = np.histogram2d(x, y, bins=(ngridx, ngridy),
                                               range=(rangedat[varj, :], rangedat[vari, :]))
            H = H.transpose()
            xedges = (xedges[:-1] + xedges[1:]) / 2
            yedges = (yedges[:-1] + yedges[1:]) / 2
            mxx, mxy = np.meshgrid(xedges, yedges)
            plt.contourf(mxx, mxy, H, 100, cmap='Blues')
            plt.title(f"{labels[vari]}-{labels[varj]}")

            # Confidence levels
            indx_g, indy_g = np.meshgrid(np.arange(0, ngridx), np.arange(0, ngridy))
            vmx2 = np.squeeze(np.reshape(mxx, (-1, 1)))
            vmy2 = np.squeeze(np.reshape(mxy, (-1, 1)))
            vm = np.squeeze(np.reshape(H, (-1, 1)))
            vx = np.squeeze(np.reshape(indx_g, (-1, 1)))
            vy = np.squeeze(np.reshape(indy_g, (-1, 1)))
            vm2 = np.sort(vm)[::-1]
            ix = np.argsort(vm, axis=0)[::-1]
            ix = np.ix_(ix)
            vx2 = vx[ix]
            vy2 = vy[ix]
            vmx2 = vmx2[ix]
            vmy2 = vmy2[ix]
            vm2 = np.cumsum(vm2 / np.sum(vm2))
            cmxx2 = H.copy()
            mxx2 = mxx.copy()
            mxy2 = mxy.copy()
            for ki in range(0, len(vm2)):
                mxx2[vy2[ki], vx2[ki]] = vmx2[ki]
                mxy2[vy2[ki], vx2[ki]] = vmy2[ki]
                cmxx2[vy2[ki], vx2[ki]] = vm2[ki]
            conls = plt.contour(mxx2, mxy2, cmxx2, levels, colors='k')
            plt.clabel(conls, inline=1, fontsize=10)
            ax.get_yaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)
            yloc = plt.MaxNLocator(max_yticks)
            ax.xaxis.set_major_locator(yloc)
            xloc = plt.MaxNLocator(max_yticks)
            ax.yaxis.set_major_locator(xloc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot posteriors from Bayesian E_iso function inference')
    parser.add_argument('-f', action='store', dest='fname', type=str,
                        help='Input MultiNest output prefix')
    parser.add_argument('-o', action='store', dest='fout', type=str,
                        help='Output plot file')
    parser.add_argument('-title', action='store', dest='title', type=str,
                        help='Plot title')
    parser.add_argument('-up', dest='bolupp', action='store_true', default=False,
                        help='Use flat prior parameter names')
    parser.add_argument('-bo', dest='bolout', action='store_true', default=False,
                        help='Save plot to file')
    args = parser.parse_args()

    fname = args.fname
    fout = args.fout
    tit = args.title
    bolupp = args.bolupp
    bolout = args.bolout

    # Read MultiNest output
    a = pymultinest.Analyzer(n_params=6, outputfiles_basename=fname)
    b = a.get_equal_weighted_posterior()

    vlik = b[:, -1]
    allres = b[:, :-1]
    mxchain = allres

    # Transform phis to log space
    mxchain[:, 0] = np.log10(mxchain[:, 0])

    # Parameter ranges for E_iso function (CHIME calibrated)
    if bolupp:
        rdat = np.array([[0, 4], [-3.0, 1.0], [39, 43], [37, 42], [-0.4, 0.4], [0, 0.8]])
    else:
        rdat = np.array([[-3, 4], [-3.0, 1.0], [39, 43], [37, 42], [-0.4, 0.4], [0, 0.8]])

    vpar = mxchain[np.argmax(vlik), :]
    print(f'Best fit: {vpar}')

    labels = [r'$\log\,\phi^*$', r'$\alpha$', r'$\log E^*_{\rm iso}$',
              r'$\log E_{\rm iso,0}$', r'$\mu_w$', r'$\sigma_w$']

    if bolout:
        import os
        os.makedirs(os.path.dirname(fout) if os.path.dirname(fout) else '.', exist_ok=True)
        plt.figure(figsize=(20, 16))
        plot2dposterier_withconf(mxchain, vlik, par=vpar, indx=[], rangedat=rdat,
                                 rate=1.0, levels=[0.68, 0.95], labels=labels)
        plt.suptitle(tit)
        plt.savefig(fout, dpi=150, bbox_inches='tight')
        print(f'Saved plot to {fout}')
    else:
        plt.figure(figsize=(20, 16))
        plot2dposterier_withconf(mxchain, vlik, par=vpar, indx=[], rangedat=rdat,
                                 rate=1.0, levels=[0.68, 0.95], labels=labels)
        plt.suptitle(tit)
        plt.show()
