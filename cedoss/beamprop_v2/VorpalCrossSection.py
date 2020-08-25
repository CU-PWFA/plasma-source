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
path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/LinearGradient_ControlCase/'
path = '/home/chris/Desktop/NERSC_LIN_Aug/'

npcase = 4e17
vmax = 1e18*(npcase/3.7e17)
vmin = 1e16*(npcase/3.7e17)
central_off = -65#-100
"""
vmax = 3.5e17
vmin = 2.5e17
central_off = 200
"""
tranExtent = 75
params = {'vmin' : vmin,
          'vmax' : vmax,
          'plasma' : 'electrons',
          'dumpInd' : 2,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : True
          }
#plot.wake_cross_section_efield(params)
plot.wake_cross_section(params)

theta,r = plot.wake_cross_section(params)
#r, theta, rhoPol = plot.wake_cross_section(params)

# Define model function to be used to fit to the data above:
def coslike(x, *p):
    rp, b, re = p
    return rp-re*np.cos(x+b)
        
# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [1., 0., 1.]
        
xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
        
# Get the fitted curve
hist_fit = coslike(theta, *xcoeff)

plt.plot(theta,r)
plt.plot(theta,hist_fit)
plt.title("Wake Shape in Polar Coordinates")
plt.xlabel("Angle (rad)")
plt.ylabel("Radius of Wake Edge (um)")
plt.grid(); plt.show()
