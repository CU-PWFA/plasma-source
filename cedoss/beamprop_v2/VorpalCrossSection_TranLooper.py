#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 12:59:21 2021

Wanted to plot the sheath density vs theta

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


#path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_Sep_Control2/'
#ind=5
path = '/home/chris/Desktop/NERSC_Sep_Grad/'
ind=9
#path = '/home/chris/Desktop/NERSC_Mar_Grad0.001/'
#ind = 10
#path = '/home/chris/Desktop/NERSC_Mar_Grad0.01/'
#ind = 10
central_off = 0
#path = '/home/chris/Desktop/NERSC_Dec_Grad/'
#ind=9
npcase = 2e16#
vmax = 1e18*(npcase/3.7e17)
vmin = 1e16*(npcase/3.7e17)
tranExtent = 200


"""
path = '/home/chris/Desktop/NERSC_Jan_Grad/'
npcase = 1e16#
tranExtent = 400
central_off = 100
ind = 8
vmax = 1e18*(npcase/3.7e17)
vmin = 1e16*(npcase/3.7e17)
"""
"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/LowRes_HighGradient_August/LinearGradient/'
npcase = 3e17#
tranExtent = 75
central_off = -10
ind = 5
vmax = 1e18*(npcase/3.7e17)*1.8
vmin = 1e16*(npcase/3.7e17)*90
"""

"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_LIN_Aug/'
npcase = 3e17#
tranExtent = 75
ind = 7
central_off=-10
vmax = 1e18*(npcase/3.7e17)
vmin = 1e18
"""

"""
path = '/home/chris/Desktop/VELA_Oct_Grad/'
npcase = 2e16#
tranExtent = 200
central_off = 0
ind = 5
vmax = 1e18*(npcase/3.7e17)
vmin = 1e16*(npcase/3.7e17)
"""
#82 for 2e16, 122 for 1e16, 37 for 3e17
radmax = 86

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
          'plot' : False,
          'radmax' : radmax,
          'drive' : 'rhoBeam'
          }

x,y,n = plot.wake_cross_section_sheath(params)
#r, theta, rhoPol = plot.wake_cross_section(params)
sfit = np.polyfit(y,n,1)
yfull = np.linspace(min(y),max(y),20)

plt.scatter(y,n)
plt.plot(yfull,yfull*sfit[0]+sfit[1],c='r')
plt.grid(); plt.show()
print("n_s0 = ",sfit[1] , " (cm^-3)")
print("dn/dy_s = ",sfit[0]*1e4 , " (cm^-4)")
