#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 1 14:06:19 2021

Loop over the cross section to track for all longitudinal pos:

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize
import scipy.integrate as integrate

#path = '/home/chris/Desktop/NERSC_Sep_Grad/'
#ind = 9
#start = 97 #um
#trunc = 0#15
#offset_arr = np.arange(-120,100.5,2)
path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g8e17/'
ind=5
offset_arr = np.arange(-140,90.5,2)
tranExtent = 200
start = 90
trunc = 0#20
npcase = 2e16
"""
rmajor_arr = np.zeros(len(offset_arr))
rminor_arr = np.zeros(len(offset_arr))

for i in range(len(offset_arr)):
    vmax = 1e18*(npcase/3.7e17)
    vmin = 1e16*(npcase/3.7e17)
    central_off = offset_arr[i]
    
    
    params = {'vmin' : vmin,
              'vmax' : vmax,
              'plasma' : 'electrons',
              'dumpInd' : ind,
              'path' : path,
              'simName' : 'MatchedBeams',
              'zoom' : 4.0,
              'alphaCutoff' : 0.05,
              'centralOff' : central_off,
              'tranExtent' : tranExtent,
              'plot' : False,
              'threshold' : 150
              }
    theta,r = plot.wake_cross_section(params)
    #r, theta, rhoPol = plot.wake_cross_section(params)
    
    # Define model function to be used to fit to the data above:
    def coslike(x, *p):
        rp, b, re = p
        return rp-re*np.cos(x+b)
            
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 0., 1.]
            
    xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
    rmajor_arr[i] = xcoeff[0]
    rminor_arr[i] = xcoeff[2]
"""
dx=1.2#0.25
dx=1.1972789115646258


x = offset_arr[trunc:]*dx*-1+47
y = rminor_arr[trunc:]
#h = np.polyfit(x,y,1)

z = offset_arr[trunc:]*dx*-1+start
"""
rp = 84.81
rm = 76.36
halfloc = 100
"""
rp = 80.41
rm = 72.58
halfloc = 114
emp = 1/2*(rp-rm)/(halfloc)*z

fig, (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
fig.set_size_inches(5,5)

pl1 = ax0.plot(offset_arr*dx*-1+start,rmajor_arr+rminor_arr,c='k',label="Greater")
ax0.plot(offset_arr*dx*-1+start,rmajor_arr-rminor_arr,c='k',ls='--',label="Lesser")
#ax0.set_title("Wake radius in first bucket")
#sax0.set_xlabel("Distance behind drive beam "+r'$(\mu m)$')
ax0.set_ylabel("Sheath Radius "+r'$(\mu m)$')
ax0.text(-27,77,"(a)",color='black')
ax0.legend(loc = 8)
#ax0.grid()

pl2 = ax1.plot(z,rminor_arr[trunc:], label = 'Simulation',c='b')
#plt.plot(z,x*h[0]+h[1],ls='--', label = 'Fit')
ax1.plot(z,emp, ls='--',label ='Model',c='red')
ax1.set_xlabel("Distance Behind Drive Beam "+r'$(\mu m)$')
ax1.set_ylabel("Wake Vertical Offset "+r'$(\mu m)$')
ax1.text(-27,14.7,"(b)",color='black')
ax1.legend(loc=9)
#ax1.grid()

#ax0.get_shared_x_axes().join(ax0, ax1)
fig.subplots_adjust(hspace=0)

plt.show();











