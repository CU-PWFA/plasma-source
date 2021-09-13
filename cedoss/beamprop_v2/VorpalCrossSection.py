#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 11:23:30 2020

Script to call a function that analyzes the shape of a wake using a transverse
slice in a rhoPlasma file.

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#path = '/home/chris/Desktop/WakeShape_LinGrad/'
#path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/LinearGradient/'
#path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/LinearGradient_ControlCase/'
#path = '/home/chris/Desktop/NERSC_Sep_Control/'
#path = '/home/chris/Desktop/NERSC_LIN_Aug/'
#path = '/media/chris/New Volume/VSimRuns/NERSCResults/Paper1/FACETII_NERSC_Run4/'
#path = '/home/chris/Desktop/NERSC_Sep_Control2/'

#path = '/home/chris/Desktop/NERSC_Sep_Grad/'
#path = '/home/chris/Desktop/NERSC_Dec_Grad/'
#path = '/home/chris/Desktop/VELA_Oct_Grad/'

#npcase = 2e16#4e17

#central_off = 0
"""
vmax = 3.5e17
vmin = 2.5e17
central_off = 200
"""
#tranExtent = 200

"""
#path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_Sep_Control2/'
#ind=5
path = '/home/chris/Desktop/NERSC_Sep_Grad/'
ind=9
path = '/home/chris/Desktop/NERSC_Mar_Grad0.001/'
ind = 10
path = '/home/chris/Desktop/NERSC_Mar_Grad0.01/'
ind = 10
central_off = 0
#path = '/home/chris/Desktop/NERSC_Dec_Grad/'
#ind=9

npcase = 2e16#
tranExtent = 200

vmax = 1e18*(npcase/3.7e17)
vmin = 1e16*(npcase/3.7e17)
"""

"""
path = '/home/chris/Desktop/NERSC_Jan_Grad/'
npcase = 1e16#
tranExtent = 400
central_off = 100
ind = 8
"""
"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/LowRes_HighGradient_August/LinearGradient/'
npcase = 3e17#
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
"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_LIN_Aug/'
npcase = 3e17#
tranExtent = 75
ind = 7
central_off=-10
"""

path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_Jun_Grad_High/'
npcase = 3e17#
tranExtent = 75
central_off = 50
ind = 6
threshold = 40

path = '/home/chris/Desktop/NERSC_Deflection_Test3/'
central_off = 0
ind = 5
npcase = 2e16#
tranExtent = 200
threshold = 100

path = '/home/chris/Desktop/NERSC_LongiFieldFix1/'
ind = 3
tranExtent = 200
centrol_off = 0
threshold = 100

vmax = 1e18*(npcase/3.7e17)
vmin = 1e16*(npcase/3.7e17)


vmax = 1e18*(npcase/3.7e17)#*1.8
vmin = 1e16*(npcase/3.7e17)#*90

params = {'vmin' : vmin,
          'vmax' : vmax,
          'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : True,
          'drive' : 'rhoBeam',
          'threshold' : threshold #how far out in r to consider wake sheath
          }
plot.drive_witness_density_new(params)
#plot.wake_cross_section_efield(params)
#plot.wake_cross_section(params)

theta,r = plot.wake_cross_section(params)
#r, theta, rhoPol = plot.wake_cross_section(params)

# Define model function to be used to fit to the data above:
def coslike(x, *p):
    rp, b, re = p
    return rp-re*np.cos(x+b)
        
# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [1., 0., 1.]
        
xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
#xcoeff = np.array([  3.97932549e+01,  -4.01443535e-01,   2.79529067e-02])
# Get the fitted curve
hist_fit = coslike(theta, *xcoeff)

plt.plot(theta,r)
plt.plot(theta,hist_fit)
plt.title("Wake Shape in Polar Coordinates")
plt.xlabel("Angle (rad)")
plt.ylabel("Radius of Wake Edge (um)")
plt.grid(); plt.show()

print(xcoeff)
