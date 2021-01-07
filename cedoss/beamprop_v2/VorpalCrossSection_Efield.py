#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 11:23:30 2020

Script to call a function that analyzes the electric field in the transverse
plane of a wake

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


#path = '/home/chris/Desktop/WakeShape_LinGrad/'
#path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/LinearGradient/'
#path = '/home/chris/Desktop/NERSC_LIN_Aug/'
#path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/LinearGradient_ControlCase/'
#path = '/home/chris/Desktop/NERSC_Sep_Control/'
#path = '/home/chris/Desktop/NERSC_Sep_Control2/'

path = '/home/chris/Desktop/NERSC_Sep_Grad/'
path = '/home/chris/Desktop/VELA_Oct_Grad/'

central_off = -70#-20#-125#-20#-100
vector = 1 # 0 for Ez, 1 for Ey, 2 for Ex

#tranExtent = 75
tranExtent = 200
params = {'plasma' : 'electrons',
          'dumpInd' : 5,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : True,
          'vector' : vector,
          'field' : 'ElecFieldPlasma'
          }
x, y, evx, evy, eYZ = plot.wake_cross_section_field(params)

params['field'] = 'faceB'
params['vector'] = 2
xb, yb, evxb, evyb, eYZb = plot.wake_cross_section_field(params)

evx = evx + evxb*3e8
evy = evy - evyb*3e8

pi = np.pi
e = const.physical_constants['elementary charge'][0]
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
radius = 39.32077702e-6 #m
ybar = 0
n_cent = 2e16*100**3 #m-3
slope = 0
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
"""
"""
#Sep Grad
radius = 0#37.5e-6#79.6081268e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -8e17*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent#+5.07650248e-6
yoff = 5.07650248e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12
"""
#Sep Grad
radius = 37.5e-6#58.05997485e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -8e17*100**4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 7.36044052e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12

#Oct Grad   /4
radius = 56.10007199e-6 #m
n_cent = 2e16*100**3#*0.95 #m-3
slope = -8e17*100**4/4 #m-4
ybar = -1/4*radius**2*slope/n_cent
yoff = 1.5815867e-6 #m
ysi = y*1e-6
xsi = x*1e-6
eps = 8.854e-12

if vector == 0:
    elab = r'$E_s$'
elif vector == 1:
    elab = r'$E_y$'
    e_v_y_theory = 1/(2*eps)*e*n_cent*(ysi-ybar) + 1/(8*eps)*e*slope*(3*(ysi-ybar)**2-2*radius**2)
    e_v_x_theory = 1/(8*eps)*e*slope*(xsi**2-2*radius**2)
else:
    elab = r'$E_x$'
    e_v_y_theory = np.zeros(len(ysi))
    e_v_x_theory = 1/(2*eps)*e*n_cent*xsi

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
plt.ylim([-2e9,2e9])
plt.xlim([-10,10])
plt.xlabel("y axis "+r'$(\mu m)$')
plt.ylabel(elab+" (SI)")
plt.grid(); plt.legend(); plt.show()

plt.title(elab+" vs x")
plt.plot(x,evx,label="Simulation")
if vector != 0:
    plt.plot(x,e_v_x_theory,label="Theory")
    plt.plot(y,e_v_x_theory*0.6,label="Theory*0.6")
plt.ylim([min(evx)*1.2,max(evx)*1.2])
plt.xlabel("x axis "+r'$(\mu m)$')
plt.ylabel(elab+" (SI)")
plt.grid(); plt.legend(); plt.show()