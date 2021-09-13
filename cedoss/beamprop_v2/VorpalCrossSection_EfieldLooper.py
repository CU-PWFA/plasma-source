#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:25:05 2021

Reads in Ey over full wake, and plots Ey(0) and Ey'(0) longituninally.

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

superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'
simName = 'MatchedBeams'

"""
#path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_Sep_Control2/'
#ind = 5
path = superpath + 'NERSC_Sep_Control2/'
ind = 5
#path = superpath + 'NERSC_Sep_Grad/'
#ind = 9
#path = superpath + 'NERSC_Mar_Grad0.001/'
#ind = 10
#path = superpath + 'NERSC_Dec_Grad/'
#ind = 9
#start = 107

tranExtent = 200
dx = 1.2 #um
start = 97 #um
offset_arr = np.arange(-120,100.5,2)
"""
"""
#Jan Grads
path = superpath + 'NERSC_Jan_Grad_Control/'
ind = 8
#path = superpath + 'NERSC_Jan_Grad/'
grad = 4e17
ind = 8

offset_arr = np.arange(-80,220.5,5) #JanGrad Full
tranExtent = 400        #tranextent for the sim
threshold = 150         #for use in cross section algorithm, make larger than rp and smaller than plasma extent
npcase = 1e16           #central density for the sim
dx=1.2#0.25             #dx for the sim
start = 250 #um         #for correcting where drive beam is
trunc = 10#15              #for trimming the back of the wake
left = 20; right = 250  #for plotting r+ and r-
wanchor = -5            #for finding where r=0 to anchor yoff fit
"""
"""
#JunGrads (3e17)
path = superpath + 'NERSC_Jun_Grad_Control/'
ind = 6
#path = superpath + 'NERSC_Jun_Grad/'
#grad = 6e18
#ind = 6
#path = superpath + 'NERSC_Jun_Grad_High/'
#grad = 1.2e19
#ind = 6

offset_arr = np.arange(-70,140.5,2)
tranExtent = 75        #tranextent for the sim
threshold = 40          #for use in cross section algorithm, make larger than rp and smaller than plasma extent
npcase = 3e17           #central density for the sim
dx=0.5                  #dx for the sim
start = 57 #um         #for correcting where drive beam is
trunc = 13#15           #for trimming the back of the wake
left = 10; right = 80  #for plotting r+ and r-
wanchor = -14            #for finding where r=0 to anchor yoff fit
"""

"""
path = '/home/chris/Desktop/NERSC_LongiFieldFix1/'
ind = 6
tranExtent = 200
dx = 0.857
start = 97
offset_arr = np.arange(-120,100.5,2)
"""

"""
path = '/home/chris/Desktop/NERSC_Deflection_July/'
ind = 5
tranExtent = 200
dx = 1.12
start = 97
offset_arr = np.arange(-120,100.5,2)
simName = 'PTPL_Gradient'
"""
"""
#August Linear Gradients, n=2e16 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + '/NERSC_n2e16_g0/'
ind=5
path = superpath + 'NERSC_n2e16_g8e17/'
path = superpath + 'NERSC_n2e16_g2e17/'
path = superpath + 'NERSC_n2e16_g2.5e16/'
path = superpath + 'NERSC_n2e16_g2.5e15/'
tranExtent = 200

dx = 1.2 #um
start = 90 #um
offset_arr = np.arange(-140,90.5,2)
"""

#August Linear Gradients, n=1e16 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + '/NERSC_n1e16_g0/'
ind=5
tranExtent = 200

dx = 1.233 #um
start = 90 #um
offset_arr = np.arange(-140,90.5,2)

"""
#August Linear Gradients, n=3e17 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + '/NERSC_n1e17_g0/'
ind=5
tranExtent = 95#75

dx = 0.5 #um
start = 90 #um
offset_arr = np.arange(-210,140.5,3)
"""

vector = 0 # 0 for Ez, 1 for Ey, 2 for Ex

e0_arr = np.zeros(len(offset_arr))
e0p_arr= np.zeros(len(offset_arr))

plt.figure(figsize=(8,7))
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
    evy1 = np.flip(eYZ[:,int((Nz+1)/2)-1],0)
    evy2 = np.flip(eYZ[:,int((Nz+1)/2)-2],0)
    if axeflip:
        evy1 = np.flip(eYZ[int((Nz+1)/2)-1,:],0)
        evy2 = np.flip(eYZ[int((Nz+1)/2)-2,:],0)
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
    e0p_arr[i]=(evy[center+1]-evy[center-1])/(y[center+1]-y[center-1])*1e4
    

    
    if (i%5==0) & (i>10) :
        plt.plot(y,evy,label=str(central_off*dx*-1+start)+r'$\ \mu m$')
plt.title("Transverse Profiles of Longitudinal Electric Field")
plt.xlabel("y "+r'$(\mu m)$')
plt.ylabel(r'$\mathrm{E_z \ (V/m)}$')
plt.grid(); plt.legend(title="Dist. Behind Driver"); plt.show()

#dx = 1.2 #um
#start = 97 #um

plt.plot(offset_arr*dx*-1+start,e0_arr,label="Positive Deflection")
plt.plot(offset_arr*dx*-1+start,-e0_arr,label="Negative Deflection")
plt.yscale('log')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Electric Field at r=0 "+r'$(V/m)$')
plt.legend(); plt.grid() ;plt.show()

plt.plot(offset_arr*dx*-1+start,e0_arr)
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Electric Field at r=0 "+r'$(V/m)$')
plt.legend(); plt.grid() ;plt.show()

plt.plot(offset_arr*dx*-1+start,e0p_arr,label="Positive Gradient")
plt.yscale('log')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("E-field Gradient at r=0 "+r'$(V/m^2)$')
plt.legend(); plt.grid() ;plt.show()

#did some off-script analysis with these

def toUm(off, dx, start):
    return off*dx*-1+start

def slope(x,y,i):
    return (y[i+1]-y[i-1])/(x[i+1]-x[i-1])

#offset_arr[60]*dx*-1+start
#offset_arr[60]
#e0_arr[60]
#slope(toUm(offset_arr,dx,start),e0_arr,60)
