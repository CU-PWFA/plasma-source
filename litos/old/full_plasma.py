#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:12:40 2017

@author: litos
"""


## common libraries
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time as time

## custom functions
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
np0  = 5e16 # cm^-3, plasma number density
wp0  = (5.64e4)*np.sqrt(np0) # rad/s, plasma ang. freq.
kp0  = wp0/c # m^-1, plasma wave number
lp0  = 0.1189 # m, ramp length
dgds = 7827 # energy gain / m
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
gb0    = 20000 # Lorentz factor
eps0   = (5e-6) # m-rad, normalized emittance
kb0    = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number
beta0  = 0.10 # m, vacuum waist beta
alpha0 = 0.0
gamma0 = (1+(alpha0**2))/beta0 # 1/m
T0     = [beta0,alpha0,gamma0] # Twiss vector

beam0 = [gb0,eps0,beta0,alpha0,gamma0]

Lp0   = 2.0 # m, dist from start to flat-top
waist = -0.3850 # m, waist relative to flat-top
beam0 = [gb0,eps0,beta0,alpha0,gamma0]
beam0 = propbeam(beam0,[0,-(Lp0+waist)])
beam0 = beam0[len(beam0)-1][:]

z = np.linspace(0,2.0,201)

pargs = [Lp0,lp0,2.0,1e-3]
npl   = plasma_ramp(np0,shape,z[0:len(z)-1],pargs)

beam1  = propbeam(beam0,z,npl,0)

i = 0
i   = int(np.round(i))
gb  = beam1[i][0]
eps = beam1[i][1]
T   = beam1[i][2:]
[x,xp,xpn] = draw_ellipse(eps/gb,T,0,True)

fig  = plt.figure()
scat = plt.scatter(x,xpn,marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(x),1.1*max(x)))
plt.ylim((1.1*min(xpn),1.1*max(xpn)))
plt.show()

def update_ellipse(i,beam,scat):
    i   = int(np.round(i))
    gb  = beam[i][0]
    eps = beam[i][1]
    T   = beam[i][2:]
    [x,xp,xpn] = draw_ellipse(eps/gb,T,0,False,51)
    xxp = np.zeros((len(x),2))
    xxpn= np.zeros((len(x),2))
    for j in range(0,len(x)):
        xxp[j][:]  = [x[j],xp[j]]
        xxpn[j][:] = [x[j],xpn[j]]
        
    scat.set_offsets(xxpn)
#    scat.set_array(gb)
#    time.sleep(0.05)
    return scat,

nstep = beam1.shape[0]

fig  = plt.figure()
scat = plt.scatter([],[],marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(x),1.1*max(x)))
plt.ylim((1.1*min(xpn),1.1*max(xpn)))

ani = animation.FuncAnimation(fig, update_ellipse,\
                              frames=np.linspace(0,nstep-1,nstep),\
                              fargs=(beam1, scat),blit=True)
plt.show()