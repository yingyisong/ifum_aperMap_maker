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

def func_parabola(x, a, b, c):
    return a*(x-b)**2+c

def func_4th(x, a, b, c, d):
    return a*x**3+b*x**2+c*x+d
    #return a*((x+b)**2+c)*((x+b)**2+d)+e

def func_linear(x, p, q):
    return p*x+q

################################################################
################################################################
def main(argv, arc):
    fnum = argv[1] #'0084'
    work_dir = "/Users/yysong/git2/ifum_aperMap_example"
    path_tab = os.path.join(work_dir,"tab")
    path_fig = os.path.join(work_dir,"fig")
    path_data_packed = os.path.join(work_dir,"data_packed")

    for side in ['b', 'r']:
        fname = side + fnum

        #### curve fit: x=a(y+b)^2+c
        path_dp = "%s/dp_%s.txt"%(path_tab,fname)
        x_dp = (2048-readFloat_space(path_dp,0))*2
        y_dp = 2056-readFloat_space(path_dp,1)  ### caution!!!

        N_dp = 8
        N_cf = np.int32(len(x_dp)/N_dp)

        print("\n#### %s side ####"%side)
        print("Fitting x=a*(y+b)^2+c")
        print("[a, b ,c]")
        popt_cf = []
        popt_cf1 = []
        for i in range(N_cf):
            n1 = i*8
            n2 = (i+1)*8

            popt, pcov = curve_fit(func_parabola, y_dp[n1:n2], x_dp[n1:n2])
            print(popt)
            popt_cf.append(popt)

            #popt, pcov = curve_fit(func_4th, y_dp[n1:n2], x_dp[n1:n2])
            #print(popt)
            #popt_cf.append(popt)
        popt_cf = np.array(popt_cf)
        #popt_cf1 = np.array(popt_cf1)
        popt_mean = np.mean(popt_cf, axis=0)
        print("Mean:")
        print(popt_mean)

        #### linear fit: a=p*c+q
        popt_lf, pcov_lf = curve_fit(func_linear, popt_cf[:,2], popt_cf[:,0])
        print("\nFitting a=p*x+q")
        print("[p, q]")
        print(popt_lf)

        #### record fitting results to file
        foutput = open("%s/cf_%s.txt"%(path_tab,fname),'w')
        foutput.write("# x=a*(y-b)^2+c\n")
        foutput.write("# a=p*c+q\n")
        foutput.write("# b, p, q\n")
        foutput.write("%.1f, %.6e, %.6e\n"%(popt_mean[1], popt_lf[0], popt_lf[1]))
        foutput.close()

        #### load a packed fits file
        path_fits = "%s/%s.fits"%(path_data_packed,fname)
        hdul_fits = fits.open(path_fits)
        hdr_fits  = hdul_fits[0].header
        data_fits = np.float64(hdul_fits[0].data)

        ####
        yy = np.arange(hdr_fits['NAXIS2'])

        #### plot 1
        fig = plt.figure(1, figsize=(8,8))
        fig.show()
        fig.clf()
        ax = fig.add_subplot(111)
        ax.imshow(data_fits, origin='lower')
        ax.plot(x_dp,y_dp,'rx')
        for i in range(N_cf):
            xx = func_parabola(yy, popt_cf[i][0], popt_cf[i][1], popt_cf[i][2])
            #xx = func_4th(yy, popt_cf[i][0], popt_cf[i][1], popt_cf[i][2], popt_cf[i][3])
            ax.plot(xx, yy, 'r-')
            xx = func_parabola(yy, popt_mean[0], popt_mean[1], popt_cf[i][2])
            #xx = func_4th(yy, popt_mean[0], popt_mean[1], popt_mean[2], popt_cf[i][3])
            ax.plot(xx, yy, 'g-')
        fig.set_tight_layout(True)
        fig.savefig('%s/cf_img_%s.pdf'%(path_fig,fname), format='pdf', transparent=True)

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
        fig.savefig('%s/cf_mean_diff_%s.pdf'%(path_fig,fname), format='pdf', transparent=True)

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
        fig.savefig('%s/cf_lf_diff_%s.pdf'%(path_fig,fname), format='pdf', transparent=True)

if __name__ == '__main__':
    main(sys.argv, len(sys.argv))
