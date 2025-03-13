#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:04:49 2021

Script to call a function that analyzes the electric field in the transverse
plane of a wake

Version for fig 3 in paper 2, using sep grad (FIG 7)

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
#path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g8e17/'
#path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g0/'
path = '/media/chris/New Volume/VSimRuns/Oct2022_FinalSims/NERSC_n2e16_g8e17/'
#path = '/media/chris/New Volume/VSimRuns/Oct2022_FinalSims/NERSC_n2e16_g0/'
ind=5
central_off = -21#-100
tranExtent = 200
c = const.speed_of_light

flip = True

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
    
    #params['vector']=0
    #xc, yc, bvxz, bvyz, bYZz = plot.wake_cross_section_field(params)
    
    Nz = len(x)
    if flip:
        if vector == 1:
            eYZ = -1*(eYZ - bYZ*c)
            eyvx_1 = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
            eyvy_1 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
            eyvx_2 = np.array(np.flip(eYZ[int((Nz+1)/2)-2,:],0))
            eyvy_2 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-2],0))
            eyvx = (eyvx_1 + eyvx_2)/2
            eyvy = (eyvy_1 + eyvy_2)/2
        else:
            eYZ = -1*(eYZ + bYZ*c)
            exvx_1 = np.array(np.flip(eYZ[int((Nz+1)/2)-1,:],0))
            exvy_1 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-1],0))
            exvx_2 = np.array(np.flip(eYZ[int((Nz+1)/2)-2,:],0))
            exvy_2 = np.array(np.flip(eYZ[:,int((Nz+1)/2)-2],0))
            exvx = (exvx_1 + exvx_2)/2
            exvy = (exvy_1 + exvy_2)/2
    else:
        if vector == 1:
            eYZ = eYZ - bYZ*c
            eyvx_1 = np.array(eYZ[int((Nz+1)/2)-1,:])
            eyvy_1 = np.array(eYZ[:,int((Nz+1)/2)-1])
            eyvx_2 = np.array(eYZ[int((Nz+1)/2)-2,:])
            eyvy_2 = np.array(eYZ[:,int((Nz+1)/2)-2])
            eyvx = (eyvx_1 + eyvx_2)/2
            eyvy = (eyvy_1 + eyvy_2)/2
        else:
            eYZ = eYZ + bYZ*c
            exvx_1 = np.array(eYZ[int((Nz+1)/2)-1,:])
            exvy_1 = np.array(eYZ[:,int((Nz+1)/2)-1])
            exvx_2 = np.array(eYZ[int((Nz+1)/2)-2,:])
            exvy_2 = np.array(eYZ[:,int((Nz+1)/2)-2])
            exvx = (exvx_1 + exvx_2)/2
            exvy = (exvy_1 + exvy_2)/2
    
    pi = np.pi
    e = const.physical_constants['elementary charge'][0]
    error= 0
    
    #Sep Grad
    radius = 76.15e-6 #m #41.97e-6 for -100
    n_cent = 2e16*100**3#*0.95 #m-3
    slope = 8e17*100**4 #m-4
    if flip:
        slope = slope*-1
    zsi_max = 112.2e-6#114e-6
    facA = .3112#0.3239
    facB = 3.19089#4.059#2.482
    delta = 0.090421#0.0651#0.134
    rpl = radius*np.sqrt((n_cent-facA*radius*slope)/n_cent) #sign flip b/c slope is wack
    rmn = radius*np.sqrt((n_cent+facA*radius*slope)/n_cent)
    yoff = zsi_max*(rpl-rmn)/(2*zsi_max)
    #yoff = 210e-6*(rpl-rmn)/(2*zsi_max)
    #yoff = 8.43e-6# for -100
    n_cent = n_cent+slope*yoff
    #radius = 41.97e-6
    ybar = 1/4*radius**2*slope/n_cent
    #yoff = 4.2188366149411554e-6 #m
    ysi = y*1e-6
    xsi = x*1e-6
    eps = 8.854e-12
    #ring=-7.17e8#-8.6e8#-3.358e8#-1.10e9#-1.21e9 #V/m
    

    

    
    sheathslope = slope*facB
    #sheathslope = 6.739e18*-1e8
    L_sh = delta*radius
    ring = -(slope - sheathslope)*(L_sh)*(radius)*e/2/eps
    
    
    if vector == 0:
        elab = r'$E_s$'
    elif vector == 1:
        elab = r'$E_y$'
        ey_v_y_theory = 1/(2*eps)*e*n_cent*(ysi-2*ybar-yoff) + 1/(8*eps)*e*slope*(3*(ysi-yoff)**2) + ring + error
        ey_v_x_theory = 1/(8*eps)*e*slope*(3*yoff**2+xsi**2) + 1/eps*e*n_cent*(-ybar-1/2*yoff) + ring + error
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
    ax0.plot(x,ex_v_x_theory/Escale,c=thecol,ls='dotted',label="Model")
ax0.set_ylim([min(exvx)*1.2/Escale,max(exvx)*1.2/Escale])
#ax0.set_xlabel("x axis "+r'$(\mu m)$')
ax0.set_ylabel(r'$E_x+cB_y$'+" (GV/m)")
#ax0.grid()
ax0.legend(title="(a)   "+r'$W_x$'+" vs x")
"""
plxy = ax1.plot(y,exvy/Escale*1e13,c=simcol,label="Simulation"+r'$\ \times \ 10^{13}$')
#ax1.set_title(r'$E_x$'+" vs y")
if vector != 0:
    ax1.plot(y,ex_v_y_theory/Escale,c=thecol,ls='dotted',label="Model")
ax1.set_ylim([min(exvy)*1.2*1e13/Escale,max(exvy)*1.2/Escale*1e13])
#ax1.set_xlabel("y axis "+r'$(\mu m)$')
#ax1.set_ylabel(r'$E_x$'+" (SI)")
#ax1.grid()
ax1.legend(title="(b)  Magnified "+r'$W_x$'+" vs y")
"""
plxy = ax1.plot(y,exvy/Escale,c=simcol,label="Simulation")
#ax1.set_title(r'$E_x$'+" vs y")
if vector != 0:
    ax1.plot(y,ex_v_y_theory/Escale,c=thecol,ls='dotted',label="Model")
ax1.set_ylim([-.8,.8])
#ax1.set_xlabel("y axis "+r'$(\mu m)$')
#ax1.set_ylabel(r'$E_x$'+" (SI)")
#ax1.grid()
ax1.legend(title="(b)   "+r'$W_x$'+" vs y")

plyx = ax2.plot(y,eyvx/Escale,c=simcol,label="Simulation")
#ax2.set_title(r'$E_y$'+" vs x")
if vector != 0:
    ax2.plot(x,ey_v_x_theory/Escale,ls='dotted',c=thecol,label="Model")
ax2.set_ylim([min(eyvx)*1.2/Escale,max(eyvx)*1.2/Escale])
ax2.set_xlabel("x Axis "+r'$(\mu m)$')
ax2.set_ylabel(r'$E_y-cB_x$'+" (GV/m)")
#ax2.grid()
ax2.legend(title="(c)   "+r'$W_y$'+" vs x")

plyy = ax3.plot(y,eyvy/Escale,c=simcol,label="Simulation")
#ax3.set_title(r'$E_y$'+" vs y")
if vector != 0:
    ax3.plot(y,ey_v_y_theory/Escale,ls='dotted',c=thecol,label="Model")
ax3.set_ylim([min(eyvy)*1.2/Escale,max(eyvy)*1.2/Escale])
ax3.set_xlabel("y Axis "+r'$(\mu m)$')
#ax3.set_ylabel(r'$E_y$'+" (SI)")
#ax3.grid()
ax3.legend(title="(d)   "+r'$W_y$'+" vs y")
fig.subplots_adjust(hspace=0)
plt.savefig("/home/chris/Desktop/figs/fig7.eps",bbox_inches='tight')
plt.show()