#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:13:32 2022

Script to call a function that analyzes the electric field in the transverse
plane of a wake

The sequel that uses the August 2021 simulations.  Also made the case selection
very nice.  Such high quality coding.

The next version is mostly looking at Jz dumps and seeing if there is a connection with rho

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


#August 2021 n=2e16 runs
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + 'JzDump/'
grad = 0
yoff = 0 #m
radius = 76.149e-6 #m
ind = 5
tranExtent = 200        #tranextent for the sim
npcase = 2e16           #central density for the sim
dx=1.2                  #dx for the sim
central_off = -20
simname = 'MatchedBeams'
Jz_plasma = 'JPlasma'
rho_plasma = 'rhoPlasma'

grad = 8e17
yoff = 3.962e-6 #m
radius = 76.435e-6 #m

c = const.physical_constants['speed of light in vacuum'][0]
e = const.physical_constants['elementary charge'][0]

vector = 0 # 0 for Ez, 1 for Ey, 2 for Ex
flip = True
c = const.speed_of_light

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
          'field' : Jz_plasma
          }
x, y, evx, evy, JzYZ = plot.wake_cross_section_field(params)
params['field'] = rho_plasma
params['vector'] = 0
xb, yb, bvx, bvy, rhoYZ = plot.wake_cross_section_field(params)

sumDen = -(rhoYZ-JzYZ/c) / (e*npcase*100**3) - 1

plt.imshow(sumDen)
plt.show()

plt.plot(y,np.flip(sumDen[:,int(len(evx)/2)],0))
plt.plot(y,np.zeros(len(y)),c='k',ls='--')
plt.title("Normalized Density Along Y Axis")
plt.ylabel(r'$-(\rho-J_z/c) \ \mathrm{(en_p)}$')
plt.xlabel("y (um)")
plt.show()

plt.plot(x,np.flip(sumDen[int(len(evy)/2),:],0))
plt.plot(y,np.zeros(len(y)),c='k',ls='--')
plt.title("Normalized Density Along X Axis")
plt.ylabel(r'$-(\rho-J_z/c) \ \mathrm{(en_p)}$')
plt.xlabel("x (um)")
plt.show()

sys.exit()

"""
if vector == 1:
    evx = evx - bvx*3e8
    evy = evy - bvy*3e8
else:
    evx = evx + bvx*3e8
    evy = evy + bvy*3e8
"""
Nz = len(x)
if flip:
    if vector == 1:
        eYZ = -1*(eYZ - bYZ*c)
        eyvx_1 = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
        eyvy_1 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
        eyvx_2 = np.array(np.flip(eYZ[int((Nz+1)/2)-2,:],0))
        eyvy_2 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-2],0))
        evx = (eyvx_1 + eyvx_2)/2
        evy = (eyvy_1 + eyvy_2)/2
    else:
        eYZ = -1*(eYZ + bYZ*c)
        exvx_1 = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
        exvy_1 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
        exvx_2 = np.array(np.flip(eYZ[int((Nz+1)/2)-2,:],0))
        exvy_2 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-2],0))
        evx = (exvx_1 + exvx_2)/2
        evy = (exvy_1 + exvy_2)/2
else:
    if vector == 1:
        eYZ = eYZ - bYZ*c
        eyvx_1 = np.array(eYZ[int((Nz+1)/2)-1,:])
        eyvy_1 = np.array(eYZ[:,int((Nz+1)/2)-1])
        eyvx_2 = np.array(eYZ[int((Nz+1)/2)-2,:])
        eyvy_2 = np.array(eYZ[:,int((Nz+1)/2)-2])
        evx = (eyvx_1 + eyvx_2)/2
        evy = (eyvy_1 + eyvy_2)/2
    else:
        eYZ = eYZ + bYZ*c
        exvx_1 = np.array(eYZ[int((Nz+1)/2)-1,:])
        exvy_1 = np.array(eYZ[:,int((Nz+1)/2)-1])
        exvx_2 = np.array(eYZ[int((Nz+1)/2)-2,:])
        exvy_2 = np.array(eYZ[:,int((Nz+1)/2)-2])
        evx = (exvx_1 + exvx_2)/2
        evy = (exvy_1 + exvy_2)/2

slope = -grad*100**4 #m-4
n_cent =  npcase*100**3#*0.95 #m-3
ybar = -1/4*radius**2*slope/n_cent
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
pi = np.pi
e = const.physical_constants['elementary charge'][0]
facB = 2.482
delta = 0.134
sheathslope = slope*facB
L_sh = delta*radius
ring = -(slope - sheathslope)*(L_sh)*(radius)*e/2/eps
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