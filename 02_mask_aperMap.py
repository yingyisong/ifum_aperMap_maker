#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy import constants
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord, SkyOffsetFrame, EarthLocation
from astropy.stats import median_absolute_deviation
from astropy.stats import biweight

from scipy.optimize import curve_fit
#from scipy import interpolate
#from scipy import optimize
#from scipy.optimize import least_squares
#from scipy.optimize import minimize
#from scipy.interpolate import UnivariateSpline
#import scipy.integrate as integrate

def readFloat_space(fileName, column1): #column number starts from 0
    x0 = []
    for line in open(fileName, 'r'):
        #skip over empty lines or lines starting with spaces or spaces+#
        tempString=line.strip()
        if (tempString[0] == '#' or tempString[0]==' '):
            continue

        line = line.split()
        if line[column1]=='-':
            line[column1]=float('nan')

        x = np.float64(line[column1])
        x0.append(x)
    return np.array(x0, dtype='float64')

def func_parabola_var(x, c, b, p, q):
    a = p*c+q
    return a*(x-b)**2+c

################################################################
################################################################
work_dir = "/Users/yysong/git2/ifum_aperMap_example"


map_sub = "_LSB_1x2_Halpha_Li"
#map_dir = os.path.join(work_dir,"aperMap")

para_sub = "0183"
para_dir = os.path.join(work_dir,"tab")

band_sub = 'full'

fits_dir = os.path.join(work_dir,"data_packed")
fig_dir = os.path.join(work_dir,"fig")
mask_dir = os.path.join(work_dir,"img_mask")

for side in ['b', 'r']:

    #### load bestfit parameters: b, p, q
    path_para = para_dir+"/cf_%s.txt"%(side+para_sub)
    for line in open(path_para, 'r'):
        tempString = line.strip()
        if (tempString[0] == '#' or tempString[0]==''):
            continue

        line = line.split(',')
        b_cf, p_cf, q_cf = np.float32(line[0]), np.float64(line[1]), np.float64(line[2])

    #### load band range: c1, c2
    path_band = para_dir+"/%s_%s.txt"%(band_sub, side+para_sub)
    for line in open(path_band, 'r'):
        tempString = line.strip()
        if (tempString[0] == '#' or tempString[0]==''):
            continue

        line = line.split(',')
        c1, c2 = np.float32(line[0]), np.float64(line[1])

    # #### load the AperMap fits file
    # path_map = map_dir+"/ap%s.fits"%(side+map_sub)
    # hdul_map = fits.open(path_map)
    # hdr_map  = hdul_map[0].header
    # data_map = np.int32(hdul_map[0].data)

    # N1 = hdr_map["NAXIS1"]
    # N2 = hdr_map["NAXIS2"]

    # ####
    # for yy in range(N2):
    #     x1 = func_parabola_var(np.float32(yy)+1., c1, b_cf, p_cf, q_cf)
    #     x2 = func_parabola_var(np.float32(yy)+1., c2, b_cf, p_cf, q_cf)
    #     data_map[yy,:np.int32(x1)] = 0
    #     data_map[yy,np.int32(x2)-1:] = 0

    # hdu_band = fits.PrimaryHDU(data_map, header=hdr_map)
    # hdu_band.writeto(map_dir+"/ap%s_%s.fits"%(side+map_sub,band_sub), overwrite=True)

    #### load a packed fits file for testing
    path_fits = fits_dir+"/%s.fits"%(side+para_sub)
    hdul_fits = fits.open(path_fits)
    hdr_fits  = hdul_fits[0].header
    data_fits = np.int32(hdul_fits[0].data)
    yy = np.arange(hdr_fits['NAXIS2'])

    x1 = func_parabola_var(yy, c1, b_cf, p_cf, q_cf)
    x2 = func_parabola_var(yy, c2, b_cf, p_cf, q_cf)

    #### write to img_mask
    foutput = open("%s/img_mask_%s%s"%(mask_dir,side,map_sub),'w')
    foutput.write("# The region to include in the final fits file\n")
    foutput.write("#X1   X2   Y1   Y2\n")
    for i in range(len(yy)):
        foutput.write("%.0f %.0f %.0f %.0f\n"%(x1[i],x2[i],yy[i],yy[i]))
    foutput.close()

    #### plot 1
    fig = plt.figure(1, figsize=(8,8))
    fig.show()
    fig.clf()
    ax = fig.add_subplot(111)
    ax.imshow(data_fits, origin='lower')

    ax.plot(x1, yy, 'r-')
    ax.plot(x2, yy, 'r-')

    fig.set_tight_layout(True)
    fig.savefig('%s/img_%s_%s.pdf'%(fig_dir,side+para_sub,band_sub), format='pdf', transparent=True)



'''
    #### plot 2
    fig = plt.figure(2, figsize=(8,6))
    fig.show()
    fig.clf()
    ax = fig.add_subplot(111)
    ax.set_xlabel("$x_\mathrm{cf}-x_\mathrm{cf,mean}$")
    ax.set_ylabel("$y$")
    ax.set_xlim([-6,6])
    ax.set_ylim([0,np.int32(hdr_fits['NAXIS2'])])

    ax.axvline(x=0,ls=':',c='k',label='$a=%.2e$\n$b=%.1f$'%(popt_mean[0], popt_mean[1]))
    ax.legend(loc='center right')

    for i in range(N_cf):
        n1 = i*8
        n2 = (i+1)*8
        xx = func_parabola(y_dp[n1:n2], popt_cf[i][0], popt_cf[i][1], 0)
        xx0 = func_parabola(y_dp[n1:n2], popt_mean[0], popt_mean[1], 0)
        #xx = func_4th(y_dp[n1:n2], popt_cf[i][0], popt_cf[i][1], popt_cf[i][2], popt_cf[i][3])
        #xx0 = func_4th(y_dp[n1:n2], popt_mean[0], popt_mean[1], popt_mean[2], popt_cf[i][3])
        ax.plot(xx-xx0, y_dp[n1:n2], 'x--')
    fig.set_tight_layout(True)
    fig.savefig('./fig/cf_mean_diff_%s.pdf'%(fname), format='pdf', transparent=True)

    #### plot 3
    fig = plt.figure(3, figsize=(8,6))
    fig.show()
    fig.clf()
    ax = fig.add_subplot(111)
    ax.set_xlabel("$x_\mathrm{cf}-x_\mathrm{cf,lf}$")
    ax.set_ylabel("$y$")
    ax.set_xlim([-1,1])
    ax.set_ylim([0,np.int32(hdr_fits['NAXIS2'])])

    ax.axvline(x=0,ls=':',c='k',label='$a=(%.2e) \cdot x+(%.2e)$\n$b=%.1f$'%(popt_lf[0],popt_lf[1], popt_mean[1]))
    ax.legend(loc='center right')

    for i in range(N_cf):
        n1 = i*8
        n2 = (i+1)*8
        xx = func_parabola(y_dp[n1:n2], popt_cf[i][0], popt_cf[i][1], 0)
        xx0 = func_parabola(y_dp[n1:n2], func_linear(popt_cf[i][2],popt_lf[0],popt_lf[1]), popt_mean[1], 0)
        ax.plot(xx-xx0, y_dp[n1:n2], 'x--')
    fig.set_tight_layout(True)
    fig.savefig('./fig/cf_lf_diff_%s.pdf'%(fname), format='pdf', transparent=True)
'''
