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

superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'
simName = 'MatchedBeams'

superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + '/NERSC_n2e16_g0/'
ind=5
path = superpath + 'NERSC_n2e16_g8e17/'
tranExtent = 200

dx = 1.2 #um
start = 90 #um
offset_arr = np.arange(-140,90.5,2)

vector = 0

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
    
    colors = plt.cm.jet(np.linspace(0, 1, 10))
    
    if (i%10==0) & (i>10) :
        plt.plot(y,evy,c=colors[k],label=str(central_off*dx*-1+start)+r'$\ \mu m$')
        k = k+1
#plt.title("Transverse Profiles of Longitudinal Electric Field")
plt.xlabel("y "+r'$(\mu m)$')
plt.ylabel(r'$\mathrm{E_z \ (V/m)}$')
plt.legend(loc=3,fontsize='small')#,title="Distance Behind Driver"); plt.show()

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
