#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 11:33:06 2021

Script to call a function that analyzes the electric field in the transverse
plane of a wake

Copy of crosssection_efield, but here I want to just mess around with analysis of the longitudinal fields.

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

"""
path = '/home/chris/Desktop/NERSC_Jan_Grad/'
central_off = 100
tranExtent = 400
ind=8
"""

#path = '/home/chris/Desktop/NERSC_Sep_Control2/'
path = superpath+'NERSC_Sep_Control2/'
ind = 5
#path = superpath + 'NERSC_Sep_Grad/'
#ind = 9
#path = superpath+'NERSC_Mar_Grad0.001/'
#ind = 10
#path = superpath+'NERSC_Mar_Grad0.01/'
#ind = 10
tranExtent = 200
central_off = -30#-100

"""
path = '/home/chris/Desktop/NERSC_Dec_Grad/'
tranExtent = 200
ind = 9
central_off = -4
"""
"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/LowRes_HighGradient_August/LinearGradient/'
tranExtent = 75
central_off = -10
ind = 5
"""
"""
path = '/home/chris/Desktop/VELA_Oct_Grad/'
npcase = 2e16#
tranExtent = 200
central_off = 0
ind = 5
"""
vector = 0 # 0 for Ez, 1 for Ey, 2 for Ex

params = {'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : True,
          'vector' : vector,
          'field' : 'edgeE'#'ElecFieldPlasma'
          }
x, y, evx, evy, eYZ = plot.wake_cross_section_field(params)
"""
params['field'] = 'faceB'
if vector == 1:
    params['vector'] = 2
else:
    params['vector'] = 1
xb, yb, bvx, bvy, bYZ = plot.wake_cross_section_field(params)

params['vector']=0
xc, yc, bvxz, bvyz, bYZz = plot.wake_cross_section_field(params)

if vector == 1:
    eYZ = eYZ - bYZ*3e8
elif vector == 2:
    eYZ = eYZ + bYZ*3e8
"""
Nz = len(x)

evx = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
evy = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))


pi = np.pi
e = const.physical_constants['elementary charge'][0]
ring=0
error= 0

"""
#Aug LinGrad
#radius = 29.17866591e-6 #m
radius = 30.13476625e-6 #m
ybar = 0.5137*1e-6
#ybar = 0#2.25574658*1e-6
n_cent = 3e17*100**3 #m-3
slope = -6e18*100**4 #m-4
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
"""
"""
#SepControl2
radius = 80.82e-6 #m
ybar = 0
n_cent = 2e16*100**3 #m-3
slope = 0
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
yoff=0
error=8.633e7
"""

#Sep Grad
radius = 80.70044626e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -8e17*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 4.2188366149411554e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
ring=-1.21e9 #V/m

"""
#Mar Grad 0.001
radius = 80.5654136e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -2.5e15*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 0.0135e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
ring=-1.8e7 #V/m
error = -1.2e8 #V/m
"""
"""
#Mar Grad 0.01
radius = 80.54784833e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -2.5e16*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 0.10223593e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
ring=-1.8e7 #V/m
error = -1.35e8 #V/m
"""
"""
#Dec Grad
radius = 82.77904831e-6 #m
n_cent = 2e16*100**3 #m-3
slope = -8e17*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 4.4568012e-6 #m*
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
error=0
"""
"""
#Oct Grad   /4
radius = 81.85e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -2e17*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 1.11566955e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
ring = -2.3e8
"""
"""
#Jan Grad ind = +100
radius = 103.933674e-6#m
#n_cent = 1e16*100**3#*0.95 #m-3
n_cent = 9.712e15*100**3#*0.95 #m-3 #at y_off
slope = -4e17*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent#*0.7
yoff = 7.19227654e-6#m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
ring = 0#-1231003612.89
#error= 575263315.6091034
"""
"""
#LowResHighGrad
radius = 30.57525162e-6 #m
n_cent = 3e17*100**3 #m-3
slope = -3e19*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 1.87193507e-6 #m*
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12

error = 0
error = -slope*1.5e-17*.12
"""
if vector == 0:
    elab = r'$E_s$'
elif vector == 1:
    elab = r'$E_y$'
    #Original Terrible ones
    #e_v_y_theory = 1/(2*eps)*e*n_cent*(ysi-ybar) + 1/(8*eps)*e*slope*(3*(ysi-ybar)**2-2*radius**2)
    #e_v_x_theory = 1/(8*eps)*e*slope*(xsi**2-2*radius**2)
    
    #This one aight, just with yoff in the eqn
    e_v_y_theory = 1/(2*eps)*e*n_cent*(ysi+2*ybar-yoff) + 1/(8*eps)*e*slope*(3*(ysi-yoff)**2) + ring + error
    e_v_x_theory = 1/(8*eps)*e*slope*(3*yoff**2+xsi**2) + 1/eps*e*n_cent*(ybar-1/2*yoff) + ring + error
    
    #e_v_y_theory = 1/(2*eps)*e*n_cent*(ysi+2*ybar-yoff) + 1/(8*eps)*e*slope*(3*(ysi-yoff)**2)
    #e_v_x_theory = 1/(8*eps)*e*slope*(xsi**2) + 1/eps*e*n_cent*(ybar)
else:
    elab = r'$E_x$'
    e_v_y_theory = np.zeros(len(ysi))
    e_v_x_theory = 1/(2*eps)*e*n_cent*xsi - 1/(4*eps)*e*slope*xsi*yoff

plt.title(elab+" vs y")
plt.plot(y,evy,label="Simulation")
if vector != 0:
    plt.plot(y,e_v_y_theory,ls='dotted',label="Theory")
    #plt.plot(y,e_v_y_theory*0.7,ls='dotted',label="Theory*0.7")
plt.ylim([min(evy)*1.2,max(evy)*1.2])
plt.xlabel("y axis "+r'$(\mu m)$')
plt.ylabel(elab+" (SI)")
plt.grid(); plt.legend(); plt.show()

plt.title(elab+" vs y")
plt.plot(y,evy,label="Simulation")
if vector != 0:
    plt.plot(y,e_v_y_theory,ls='dotted',label="Theory")
    #plt.plot(y,e_v_y_theory*0.7,ls='dotted',label="Theory*0.7")
plt.ylim([-1e9,3e9])
plt.xlim([-10,10])
plt.xlabel("y axis "+r'$(\mu m)$')
plt.ylabel(elab+" (SI)")
plt.grid(); plt.legend(); plt.show()

plt.title(elab+" vs x")
plt.plot(x,evx,label="Simulation")
if vector != 0:
    plt.plot(x,e_v_x_theory,label="Theory")
    #plt.plot(y,e_v_x_theory*0.6,label="Theory*0.6")
plt.ylim([min(evx)*1.2,max(evx)*1.2])
plt.xlabel("x axis "+r'$(\mu m)$')
plt.ylabel(elab+" (SI)")
plt.grid(); plt.legend(); plt.show()