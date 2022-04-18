#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:27:53 2022

Script to call a function that analyzes the electric field in the transverse
plane of a wake

The sequel that uses the August 2021 simulations.  Also made the case selection
very nice.  Such high quality coding.  Please Clap

Made for analyzing the 40 dump simulation of a deflection NERSC run

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
path = superpath + 'NERSC_n2e16_g0/'
grad = 0
yoff = 0 #m
radius = 76.149e-6 #m
#ind = 21
tranExtent = 200        #tranextent for the sim
npcase = 2e16           #central density for the sim
dx=1.2                  #dx for the sim
#central_off = -20
simname = 'MatchedBeams'
efield = 'edgeE'
bfield = 'faceB'
setno = 10
if setno == 10:
    path = superpath + 'Deflection_HighDumps/'
    grad = 8e17
    yoff = 3.962e-6 #m
    radius = 76.435e-6 #m
    central_off = -90#-20#-80  #-90 to -12
    #ind = 39
    simname = 'PTPL_Gradient'
    efield = 'ElecFieldPlasma'
    bfield = 'MagFieldPlasma'

vector = 1 # 0 for Ez, 1 for Ey, 2 for Ex
flip = True
c = const.speed_of_light

ind_arr = np.arange(40)+1
ey_arr = np.zeros(len(ind_arr))
ex_arr = np.zeros(len(ind_arr))
for i in range(len(ind_arr)):
    params = {'plasma' : 'electrons',
              'dumpInd' : ind_arr[i],
              'path' : path,
              'simName' : simname,#'MatchedBeams',
              'zoom' : 4.0,
              'alphaCutoff' : 0.05,
              'centralOff' : central_off,
              'tranExtent' : tranExtent,
              'plot' : False,
              'vector' : vector,
              'field' : efield
              }
    x, y, evx, evy, eYZ = plot.wake_cross_section_field(params)

    params['field'] = bfield
    if vector == 1:
        params['vector'] = 2
    else:
        params['vector'] = 1
    xb, yb, bvx, bvy, bYZ = plot.wake_cross_section_field(params)
    
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
    ey_arr[i]=evy[int((Nz)/2)]
    ex_arr[i]=evx[int((Nz)/2)]

plt.plot(ind_arr,ey_arr,label="Two-Bunch")

#For this stuff, compare with DeflectionPlots.py after running that

start = 90
dx = 1.2
#z_arr = offset_arr*dx*-1+start
zpos = (central_off*dx*-1+start)
z_ind = np.argmin(np.abs(z_arr-zpos))

ey_0 = ey_simey0[z_ind]
ey_0_emp = ey_fullemp[z_ind]

plt.title("Tranvserse Field at z = "+str(zpos)+" um")
plt.ylabel("Field at y=0 (V/m)")
plt.xlabel("Simulation Dump Ind")

plt.plot([ind_arr[0],ind_arr[-1]],[ey_0,ey_0],ls='--',label="Single Bunch")
plt.plot([ind_arr[0],ind_arr[-1]],[ey_0_emp,ey_0_emp],ls='--',label="Empirical Model")
ey_0_int = np.trapz(ey_arr,ind_arr)*0.08
plt.plot([ind_arr[0],ind_arr[-1]],[ey_0_int,ey_0_int],ls='dotted',label="Calibrated Integral")
plt.plot([ind_arr[0],ind_arr[-1]],[0,0],c='k',ls='dotted')
plt.legend();plt.show()

"""
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
"""