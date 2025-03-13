#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 11:04:17 2021

Wanted to plot the sheath density vs theta

This for fig 2 in paper 2

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#path = '/home/chris/Desktop/NERSC_Sep_Grad/'
#ind = 9
#central_off = 0

path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g8e17/'
#path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g2.5e16/'

ind=5
central_off = -21#-100
npcase = 2e16#
tranExtent = 200

#82 for 2e16, 122 for 1e16, 37 for 3e17
radmax = 78#40
yoff = 10

"""
path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n1e17_g2e18/'

ind = 5
central_off = -33
tranExtent = 95
npcase = 1e17
radmax = 50
yoff = 5
"""
params = {'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : False,
          'radmax' : radmax,
          'yoff' : yoff,
          'drive' : 'rhoBeam',
          'npcase' : npcase/1e17
          }

x,y,n = plot.wake_cross_section_sheath(params)
y2, nvy = plot.wake_cross_section_densityslice(params)
#r, theta, rhoPol = plot.wake_cross_section(params)
sfit = np.polyfit(y,n,1)
yfull = np.linspace(min(y),max(y),20)
print("n_s0 = ",sfit[1] , " (cm^-3)")
print("dn/dy_s = ",sfit[0]*1e4 , " (cm^-4)")
"""
fig, (ax0,ax1) = plt.subplots(ncols=2,sharey=True)
fig.set_size_inches(5,5)

pl1 = ax0.scatter(y,n,s=10,marker='o', facecolors='none', edgecolors='black',label='Sheath Density')
ax0.plot(yfull,yfull*sfit[0]+sfit[1],c='r',label='Linear Fit')
ax0.set_ylabel("Density "+r'$(cm^{-3})$')
ax0.set_xlabel("y Axis "+r'$(\mu m)$')
ax0.legend()

ion = [8e17,2e16]
pl1 = ax1.plot(y2,nvy,c='k',label='Electrons')
ax1.plot(y2,y2*ion[0]/1e4+ion[1],c='b',label='Ions')
ax1.plot(y2,y2*sfit[0]+sfit[1],c='r',label='Sheath Fit')
ax1.set_xlabel("y Axis "+r'$(\mu m)$')
ax1.legend()

fig.subplots_adjust(wspace=0)
plt.show()
"""
ion = [8e17,2e16]

fig = plt.figure(figsize=(5,3.3))



ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.scatter(y,n,s=10,marker='o', facecolors='none', edgecolors='grey',label='Sheath Density')
ax.plot(y2,nvy,c='k',label='Electrons')
ax.plot(y2,y2*ion[0]/1e4+ion[1],c='b',label='Ions')
ax.plot(y2,y2*sfit[0]+sfit[1],c='r',label='Sheath Fit')
ax.set_ylabel("Density "+r'$(cm^{-3})$')
ax.set_xlabel("y Axis "+r'$(\mu m)$')
ax.set_ylim([-1e15,1.2*max(n)])
ax.legend()

plt.show()

y=np.flip(y,0)*-1
n=np.flip(n,0)
nvy = np.flip(nvy,0)
ion[0]=ion[0]*-1
sfit[0]=sfit[0]*-1

scale = 1e17

#fig = plt.figure(figsize=(5,3.3))

fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1,sharex=True)
fig.set_size_inches(5,6)

#ax0 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax0.scatter(y,n/scale,s=10,marker='o', facecolors='none', edgecolors='grey',label='Sheath Density')
ax0.plot(y2,nvy/scale,c='k',label='Electrons')
ax0.plot(y2,(y2*ion[0]/1e4+ion[1])/scale,c='b',label='Ions')
ax0.plot(y2,(y2*sfit[0]+sfit[1])/scale,c='r',label='Sheath Fit')
ax0.set_ylabel(r'$(\rho-J_z/c)/e$'+" "+r'$\mathrm{(10^{17} cm^{-3})}$')
#ax0.set_xlabel("y "+r'$\mathrm{(\mu m)}$')
ax0.set_ylim([-1e15/scale,1.4*max(n)/scale])
ax0.text(-205,1.07,'(a)   '+r'$\xi=112.2\mathrm{\ \mu m}$',color='black')
ax0.legend()

###Below here I do the back of the wake

radmax = 40;        params['radmax'] = radmax
yoff = 9;          params['yoff'] = yoff
central_off = -95; params['centralOff'] = central_off

x,y,n = plot.wake_cross_section_sheath(params)
y2, nvy = plot.wake_cross_section_densityslice(params)
sfit = np.polyfit(y,n,1)
qfit = np.polyfit(y,n,2)
yfull = np.linspace(min(y),max(y),20)

y=np.flip(y,0)*-1
n=np.flip(n,0)
nvy = np.flip(nvy,0)
#ion[0]=ion[0]*-1
sfit[0]=sfit[0]*-1
qfit[1]=qfit[1]*-1

#ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax1.scatter(y,n/scale,s=10,marker='o', facecolors='none', edgecolors='grey',label='Sheath Density')
ax1.plot(y2,nvy/scale,c='k',label='Electrons')
ax1.plot(y2,(y2*ion[0]/1e4+ion[1])/scale,c='b',label='Ions')
ax1.plot(y2,(y2*sfit[0]+sfit[1])/scale,c='r',label='Sheath Fit')
ax1.plot(y2,(np.square(y2)*qfit[0]+y2*qfit[1]+qfit[2])/scale,c='g',ls='--',label='Quadratic Fit')
ax1.set_ylabel(r'$(\rho-J_z/c)/e$'+" "+r'$\mathrm{(10^{17} cm^{-3})}$')
ax1.set_xlabel("y "+r'$\mathrm{(\mu m)}$')
ax1.set_ylim([-1e15/scale,1.2*max(n)/scale])
ax1.text(-205,1.0,'(b)',color='black')
ax1.text(-205,0.9,'      '+r'$\xi=201.0\mathrm{\ \mu m}$',color='black')
ax1.legend()

fig.subplots_adjust(hspace=0)
plt.savefig("/home/chris/Desktop/figs/fig2.eps",bbox_inches='tight')
plt.show()













