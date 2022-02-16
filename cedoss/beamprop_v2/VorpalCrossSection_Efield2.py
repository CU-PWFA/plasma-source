#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 23:49:10 2021

Script to call a function that analyzes the electric field in the transverse
plane of a wake

The sequel that uses the August 2021 simulations.  Also made the case selection
very nice.  Such high quality coding.  Please Clap

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

"""
#August 2021 n=2e16 runs
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + 'NERSC_n2e16_g0/'
grad = 0
yoff = 0 #m
radius = 76.149e-6 #m
ind = 5
tranExtent = 200        #tranextent for the sim
npcase = 2e16           #central density for the sim
dx=1.2                  #dx for the sim
central_off = -20
simname = 'MatchedBeams'
efield = 'edgeE'
bfield = 'faceB'
setno = 4
if setno == 1:
    path = superpath + 'NERSC_n2e16_g8e17/'
    grad = 8e17
    yoff = 3.962e-6 #m
    radius = 76.435e-6 #m
if setno == 10:
    path = superpath + 'NERSC_Deflection_Aug/'
    grad = 8e17
    yoff = 3.962e-6 #m
    radius = 76.435e-6 #m
    central_off = -50
    ind = 6
    simname = 'PTPL_Gradient'
    efield = 'ElecFieldPlasma'
    bfield = 'MagFieldPlasma'
elif setno == 2:
    path = superpath + 'NERSC_n2e16_g2e17/'
    grad = 2e17
    yoff = 0.946e-6 #m
    radius = 76.143e-6 #m
elif setno == 3:
    path = superpath + 'NERSC_n2e16_g2.5e16/'
    grad = 2.5e16
    yoff = 0.106e-6 #m
    radius = 76.122e-6 #m
elif setno == 4:
    path = superpath + 'NERSC_n2e16_g2.5e15/'
    grad = 2.5e15
    yoff = 0.0167e-6 #m
    radius = 76.154e-6 #m
"""

"""
#August Linear Gradients, n=1e16 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + 'NERSC_n1e16_g0/'
grad = 0
yoff = 0 #m
radius = 91.857e-6 #m
setno = 1
if setno == 1:
    path = superpath + 'NERSC_n1e16_g4e17/'
    grad = 4e17
    yoff = 5.590e-6 #m
    radius = 92.423e-6 #m
    
ind = 5
tranExtent = 200
dx = 1.233 #um
npcase = 1e16
central_off = -20
"""


#August Linear Gradients, n=1e17 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + 'NERSC_n1e17_g0/'
grad = 0
yoff = 0 #m
radius = 48.676e-6 #m
setno = 2
if setno == 1:
    path = superpath + 'NERSC_n1e17_g4e18/'
    grad = 4e18
    yoff = 1.748e-6 #m
    radius = 48.781e-6 #m
elif setno == 2:
    path = superpath + 'NERSC_n1e17_g2e18/'
    grad = 2e18
    yoff = 0.864e-6 #m
    radius = 48.702e-6 #m
    
ind=5
tranExtent = 95
dx = 0.5 #um
npcase = 1e17
central_off = -33


vector = 1 # 0 for Ez, 1 for Ey, 2 for Ex

params = {'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : simname,#'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : True,
          'vector' : vector,
          'field' : efield
          }
x, y, evx, evy, eYZ = plot.wake_cross_section_field(params)

params['field'] = 'faceB'
if vector == 1:
    params['vector'] = 2
else:
    params['vector'] = 1
xb, yb, bvx, bvy, bYZ = plot.wake_cross_section_field(params)

params['vector']=0
xc, yc, bvxz, bvyz, bYZz = plot.wake_cross_section_field(params)

"""
if vector == 1:
    evx = evx - bvx*3e8
    evy = evy - bvy*3e8
else:
    evx = evx + bvx*3e8
    evy = evy + bvy*3e8
"""
if vector == 1:
    eYZ = eYZ - bYZ*3e8
else:
    eYZ = eYZ + bYZ*3e8

Nz = len(x)

evx1 = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
evx2 = np.array(np.flip(eYZ[int((Nz+1)/2)-2,:],0))
evx = (evx1+evx2)/2
evy1 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
evy2 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-2],0))
evy = (evy1+evy2)/2

slope = -grad*100**4 #m-4
n_cent =  npcase*100**3#*0.95 #m-3
ybar = -1/4*radius**2*slope/n_cent
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
pi = np.pi
e = const.physical_constants['elementary charge'][0]
ring = -0 # -1.21e9 #V/m
error = 0

if vector == 0:
    elab = r'$E_s$'
elif vector == 1:
    elab = r'$E_y$'
    #This one aight, just with yoff in the eqn
    e_v_y_theory = 1/(2*eps)*e*n_cent*(ysi+2*ybar-yoff) + 1/(8*eps)*e*slope*(3*(ysi-yoff)**2) + ring + error
    e_v_x_theory = 1/(8*eps)*e*slope*(3*yoff**2+xsi**2) + 1/eps*e*n_cent*(ybar-1/2*yoff) + ring + error
    
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

cent = int(len(y)/2)
print("Ey_sim = ",evx[cent], "V/m")
print("Ey: theory - sim = ",e_v_x_theory[cent]-evx[cent]," V/m")