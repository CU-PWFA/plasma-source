#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:04:49 2021

Script to call a function that analyzes the electric field in the transverse
plane of a wake

Version for fig 3 in paper 2, using sep grad

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
import scipy.constants as const

#path = '/home/chris/Desktop/NERSC_Sep_Grad/'
#ind = 9
#central_off = 0
path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g8e17/'
ind=5
central_off = -20
tranExtent = 200

vecarr = np.array([1,2])
for i in range(len(vecarr)):

    vector = vecarr[i] # 0 for Ez, 1 for Ey, 2 for Ex
    
    params = {'plasma' : 'electrons',
              'dumpInd' : ind,
              'path' : path,
              'simName' : 'MatchedBeams',
              'zoom' : 4.0,
              'alphaCutoff' : 0.05,
              'centralOff' : central_off,
              'tranExtent' : tranExtent,
              'plot' : False,
              'vector' : vector,
              'field' : 'edgeE'#'ElecFieldPlasma'
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
    
    Nz = len(x)
    if vector == 1:
        eYZ = eYZ - bYZ*3e8
        eyvx = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
        eyvy = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
    else:
        eYZ = eYZ + bYZ*3e8
        exvx = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
        exvy = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
 
    pi = np.pi
    e = const.physical_constants['elementary charge'][0]
    error= 0
    
    #Sep Grad
    radius = 80.70044626e-6 #m
    n_cent = 2e16*100**3#*0.95 #m-3
    slope = -8e17*100**4 #m-4
    ybar = -1/4*radius**2*slope/n_cent
    yoff = 4.2188366149411554e-6 #m
    ysi = y*1e-6
    xsi = x*1e-6
    eps = 8.854e-12
    ring=-1.10e9#-1.21e9 #V/m
    
    if vector == 0:
        elab = r'$E_s$'
    elif vector == 1:
        elab = r'$E_y$'
        ey_v_y_theory = 1/(2*eps)*e*n_cent*(ysi+2*ybar-yoff) + 1/(8*eps)*e*slope*(3*(ysi-yoff)**2) + ring + error
        ey_v_x_theory = 1/(8*eps)*e*slope*(3*yoff**2+xsi**2) + 1/eps*e*n_cent*(ybar-1/2*yoff) + ring + error
    else:
        elab = r'$E_x$'
        ex_v_y_theory = np.zeros(len(ysi))
        ex_v_x_theory = 1/(2*eps)*e*n_cent*xsi - 1/(4*eps)*e*slope*xsi*yoff

Escale = 1e9
simcol = "b"
thecol = "r"

fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(nrows=2,ncols=2,sharex=True)
fig.set_size_inches(10,7.5)

plxx = ax0.plot(y,exvx/Escale,c=simcol,label="Simulation")
#ax0.set_title(r'$E_x$'+" vs x")
if vector != 0:
    ax0.plot(x,ex_v_x_theory/Escale,c=thecol,ls='dotted',label="Theory")
ax0.set_ylim([min(exvx)*1.2/Escale,max(exvx)*1.2/Escale])
#ax0.set_xlabel("x axis "+r'$(\mu m)$')
ax0.set_ylabel(r'$E_x$'+" (GV/m)")
#ax0.grid()
ax0.legend(title="(a)   "+r'$E_x$'+" vs x")

plxy = ax1.plot(y,exvy/Escale,c=simcol,label="Simulation")
#ax1.set_title(r'$E_x$'+" vs y")
if vector != 0:
    ax1.plot(y,ex_v_y_theory/Escale,c=thecol,ls='dotted',label="Theory")
ax1.set_ylim([min(exvy)*1.2/Escale,max(exvy)*1.2/Escale])
#ax1.set_xlabel("y axis "+r'$(\mu m)$')
#ax1.set_ylabel(r'$E_x$'+" (SI)")
#ax1.grid()
ax1.legend(title="(b)   "+r'$E_x$'+" vs y")

plyx = ax2.plot(y,eyvx/Escale,c=simcol,label="Simulation")
#ax2.set_title(r'$E_y$'+" vs x")
if vector != 0:
    ax2.plot(x,ey_v_x_theory/Escale,ls='dotted',c=thecol,label="Theory")
ax2.set_ylim([min(eyvx)*1.2/Escale,max(eyvx)*1.2/Escale])
ax2.set_xlabel("x Axis "+r'$(\mu m)$')
ax2.set_ylabel(r'$E_y$'+" (GV/m)")
#ax2.grid()
ax2.legend(title="(c)   "+r'$E_y$'+" vs x")

plyy = ax3.plot(y,eyvy/Escale,c=simcol,label="Simulation")
#ax3.set_title(r'$E_y$'+" vs y")
if vector != 0:
    ax3.plot(y,ey_v_y_theory/Escale,ls='dotted',c=thecol,label="Theory")
ax3.set_ylim([min(eyvy)*1.2/Escale,max(eyvy)*1.2/Escale])
ax3.set_xlabel("y Axis "+r'$(\mu m)$')
#ax3.set_ylabel(r'$E_y$'+" (SI)")
#ax3.grid()
ax3.legend(title="(d)   "+r'$E_y$'+" vs y")
fig.subplots_adjust(hspace=0)
plt.show()