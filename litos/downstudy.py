#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 13:52:24 2017

@author: litos
"""

## common libraries
#import sys
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

## custom functions
from fitdownramp import fitdownramp
from propbeam import propbeam
from plasma_ramp import plasma_ramp
from draw_ellipse import draw_ellipse
from calc_Bmag import calc_Bmag
from calc_M import calc_M
from calcwaist import calcwaist

## define constants
c = 3e8 # m/s, speed of light

## set plasma parameters
# use values for main accelerator plasma
np0 = 5e16 # cm^-3, plasma number density
wp0 = (5.64e4)*np.sqrt(np0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
shape = 'gauss'
#shape = 'gen_gauss'
#shape = 'trap'
#shape = 'sigmoid'
#shape = 'gomp'
#shape = 'xu3'
#shape = 'xu4'
#shape = 'xu5'

## set beam parameters
# use values for nominal vacuum waist
gb0    = 24219 # Lorentz factor
eps0   = 1.1*(5e-6) # m-rad, normalized emittance
kb0    = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number
beta0  = (1/kb0) # m
alpha0 = 0.0
gamma0 = (1+(alpha0**2))/beta0 # 1/m
T0     = [beta0,alpha0,gamma0] # Twiss vector

beam0 = [gb0,eps0,beta0,alpha0,gamma0]

z = np.linspace(0,2.0,2001)

## loop over ramp lengths
lp = np.linspace(0.01,0.25,25)
Tf = np.zeros((len(lp),3))
betaw = np.zeros(len(lp))
zw = np.zeros(len(lp))
for i in range(0,len(lp)):
    pargs = [0,lp[i],2,1e-3]
    npl   = plasma_ramp(np0,shape,z[0:len(z)-1],pargs,'down')
    
    beam = propbeam(beam0,z,npl[0:len(z)-1])
    [gb,eps,beta,alpha,gamma] = beam[len(beam)-1][:]
    Tf[i][:] = [beta,alpha,gamma]
    [betaw[i],zw[i]] = calcwaist(Tf[i])
    zw[i] = z[len(z)-1] - zw[i]

## plot
#ymax_plt = 1*beam[len(beam)-1,2]
#norm = ymax_plt/max(npl)
plt.figure()
plt.plot(lp,betaw,'.')
plt.xlabel(r'$l_p$ (m)')
plt.ylabel(r'$\beta^*$ (m)')
plt.title(r"Down ramp with %s profile" % shape)
#plt.text(0.5*max(z),0.5*ymax_plt,\
#plt.text(r"$\sigma$ = %.2f m, waist = %.2f m" % (lpfit,waistfit),\
#         verticalalignment='center',\
#         horizontalalignment='center')
#plt.ylim((0,ymax_plt))
#plt.axhline(targ_beta,color='r')
plt.show()


plt.figure()
plt.plot(lp,zw,'.')
plt.xlabel(r'$l_p$ (m)')
plt.ylabel(r'$z_w$ (m)')
plt.title(r"Down ramp with %s profile" % shape)
#plt.text(0.5*max(z),0.5*ymax_plt,\
#plt.text(r"$\sigma$ = %.2f m, waist = %.2f m" % (lpfit,waistfit),\
#         verticalalignment='center',\
#         horizontalalignment='center')
#plt.ylim((0,ymax_plt))
#plt.axhline(targ_beta,color='r')
plt.show()
