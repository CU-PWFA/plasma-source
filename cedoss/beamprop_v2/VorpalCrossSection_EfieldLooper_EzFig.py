#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 11:35:03 2022

Reads in Ey over full wake, and plots Ey(0) and Ey'(0) longituninally.

This copy is dedicated to the Panofsky-Wenzel Example for paper 2


@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as const

#First I like to run VorpalCrossSection_Looper2
# to get "rcontr_arr" for the radius evolution
# and "offset_arr*dx*-1+start" for the zsi

superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'
simName = 'MatchedBeams'

superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
superpath = '/media/chris/New Volume/VSimRuns/Oct2022_FinalSims/'
path = superpath + '/NERSC_n2e16_g0/'
ind=5
case = 1
if case == 1:
    path = superpath + 'NERSC_n2e16_g8e17/'
    dndy = 8e17 * 100**4 # m-4
    alpha = 0.0346
    g = 0.305
if case == 2:
    path = superpath + 'NERSC_n2e16_g2e17/'
    dndy = 2e17 * 100**4 # m-4
    alpha = 0.00866
    g = 0.076
if case == 3:
    path = superpath + 'NERSC_n2e16_g2.5e16/'
    dndy = 2.5e16 * 100**4 # m-4
    alpha = 0.00107
    g = 0.01
if case == 4: 
    path = superpath + 'NERSC_n2e16_g2.5e15/'
    dndy = 2.5e15 * 100**4 # m-4
    alpha = 0.000129
    g = 0.001
    
tranExtent = 200
n0 = 2e16 * 100**3 #m-3
ximax = 112.2e-6
rmax = 76.15e-6

e = 1.6022e-19
eps = 8.8542e-12

dx = 1.2 #um
start = 87 #um
offset_arr = np.arange(-141,91.5,2)

vector = 1
xory = 0

e0_arr = np.zeros(len(offset_arr))
e0p_arr= np.zeros(len(offset_arr))
k = 0
plt.figure(figsize=(5,5))
for i in range(len(offset_arr)):
    central_off = offset_arr[i]
    params = {'plasma' : 'electrons',
              'dumpInd' : ind,
              'path' : path,
              'simName' : simName,
              'zoom' : 4.0,
              'alphaCutoff' : 0.05,
              'centralOff' : central_off,
              'tranExtent' : tranExtent,
              'plot' : False,
              'vector' : vector,
              'field' : 'edgeE',#'ElecFieldPlasma'
              }
    x, y, evx, evy, eYZ = plot.wake_cross_section_field(params)
    
    params['field'] = 'faceB'
    if vector == 1:
        params['vector'] = 2
    else:
        params['vector'] = 1
    xb, yb, bvx, bvy, bYZ = plot.wake_cross_section_field(params)
    
    if vector == 1:
        eYZ = eYZ - bYZ*3e8
    elif vector ==2:
        eYZ = eYZ + bYZ*3e8
    
    Nz = len(x)
    
    #evy = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
    
    axeflip = False
    
    """
    evy1 = np.flip(eYZ[:,int((Nz+1)/2)-1],0)
    evy2 = np.flip(eYZ[:,int((Nz+1)/2)-2],0)
    for j in range(len(evy1)):
        evy[j] = (evy1[j]+evy2[j])/2
    """
    evy1 = eYZ[:,int((Nz+1)/2)-1]
    evy2 = eYZ[:,int((Nz+1)/2)-2]
    if axeflip:
        evy1 = eYZ[int((Nz+1)/2)-1,:]
        evy2 = eYZ[int((Nz+1)/2)-2,:]
    for j in range(len(evy1)):
        evy[j] = (evy1[j]+evy2[j])/2
    
    """
    #For Jan
    x = x[220:425]
    y = y[220:425]
    evx = evx[220:425]
    evy = evy[220:425]
    """
    center = int(len(evy)/2)
    e0_arr[i] = (evy[center+1]+evy[center-1])/2
    #e0_arr[i] = (evx[center+1]+evx[center-1])/2
    e0p_arr[i]=(evy[center+1]-evy[center-1])/(y[center+1]-y[center-1])*1e4
    
    #colors = plt.cm.jet(np.linspace(0, 1, 9))
    colors = plt.cm.brg(np.linspace(0, 1, 9))
    if central_off*dx*-1+start == 112.2:
        wz_0 = evy[int(len(evy)/2)]
    if (i%10==0) & (i>10) & (i<110):
        if xory == 0:
            
            plt.plot(y,evy/1e9,c=colors[k],label=str(central_off*dx*-1+start)+r'$\ \mu m$')
            """#Basic sin and cos approximation
            wz_0alt = evy[int(len(evy)/2)]
            xi = (central_off*dx*-1+start)*1e-6
            linterm_app = 1/4/eps*e*(-2*n0*alpha+0.6*dndy*rmax*np.sin(np.pi*xi/ximax))
            #linterm_app = 1/4/eps*e*(-2*n0*alpha+1.21*dndy*rmax*np.power(np.sin(np.pi/2*xi/ximax),2)*np.cos(np.pi/2*xi/ximax))
            wz_full_app = -1/8/eps*e*alpha*dndy*np.power(y*1e-6,2)+linterm_app*y*1e-6+wz_0alt
            plt.plot(y,wz_full_app,c=colors[k],ls='dotted')
            """
            
            wz_0alt = evy[int(len(evy)/2)]
            xi = (central_off*dx*-1+start)*1e-6
            rad_xi = rcontr_arr[i]*1e-6
            drad_xi = (rcontr_arr[i+1]-rcontr_arr[i-1])/((offset_arr[i+1]*dx*-1+start)-(offset_arr[i-1]*dx*-1+start))
            #print(rad_xi,drad_xi)
            #linterm_app = 1/4/eps*e*(-2*n0*alpha+1.21*rad_xi*drad_xi*dndy)
            #for the full approximation:
            linterm_app = 1/4/eps*e*(alpha**2*dndy*xi-2*n0*alpha-alpha*dndy**2*rad_xi**2/n0+2*rad_xi*dndy*drad_xi*(1-.397+alpha/n0*dndy*xi))
            wz_full_app = -1/8/eps*e*alpha*dndy*np.power(y*1e-6,2)+linterm_app*y*1e-6+wz_0alt
            plt.plot(y[90:246],wz_full_app[90:246]/1e9,c=colors[k],ls='dotted')
            
        if xory == 1:
            plt.plot(x,evx/1e9,c=colors[k],label=str(central_off*dx*-1+start)+r'$\ \mu m$')
            
            wz_0alt = evy[int(len(evy)/2)]
            xi = (central_off*dx*-1+start)*1e-6
            wz_full_app = 1/8/eps*e*alpha*dndy*np.power(x*1e-6,2)+wz_0alt
            plt.plot(x[90:246],wz_full_app[90:246]/1e9,c=colors[k],ls='dotted')
            
        k = k+1
        

if xory == 0:
    linterm = 1/4/eps*e*(alpha**2*dndy*ximax-2*n0*alpha-alpha*dndy*g*rmax)
    wz_full = -1/8/eps*e*alpha*dndy*np.power(y*1e-6,2)+linterm*y*1e-6+wz_0
    plt.plot(y,wz_full/1e9,c='black',ls='dotted',label="Analytic")
    #plt.ylim([-1.5e10,1.5e10])
    plt.xlabel("y "+r'$(\mu m)$')
if xory == 1:
    plt.xlabel("y "+r'$(\mu m)$')
#plt.title("Transverse Profiles of Longitudinal Electric Field")
#plt.ylabel(r'$\mathrm{W_z \ (10^{9} \times V/m )}$')
plt.ylabel(r'$\mathrm{W_z \ (GV/m )}$')
plt.legend(loc=3,fontsize='small')#,title="Distance Behind Driver"); 
plt.savefig("/home/chris/Desktop/figs/fig9.eps",bbox_inches='tight')

plt.show()

#dx = 1.2 #um
#start = 97 #um
"""
plt.plot(offset_arr*dx*-1+start,e0_arr,label="Positive Deflection")
plt.plot(offset_arr*dx*-1+start,-e0_arr,label="Negative Deflection")
plt.yscale('log')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Electric Field at r=0 "+r'$(V/m)$')
plt.legend(); plt.grid() ;plt.show()
"""
"""
plt.plot(offset_arr*dx*-1+start,e0_arr)
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Electric Field at r=0 "+r'$(V/m)$')
plt.legend(); plt.grid() ;plt.show()

plt.plot(offset_arr*dx*-1+start,e0p_arr,label="Positive Gradient")
plt.yscale('log')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("E-field Gradient at r=0 "+r'$(V/m^2)$')
plt.legend(); plt.grid() ;plt.show()
"""
#did some off-script analysis with these

def toUm(off, dx, start):
    return off*dx*-1+start

def slope(x,y,i):
    return (y[i+1]-y[i-1])/(x[i+1]-x[i-1])

#offset_arr[60]*dx*-1+start
#offset_arr[60]
#e0_arr[60]
#slope(toUm(offset_arr,dx,start),e0_arr,60)
