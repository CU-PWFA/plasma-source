#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 10:08:34 2021

Loop over the cross section to track for all longitudinal pos:

This looper looks at the sheath density profile
    
@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize
import scipy.integrate as integrate
import scipy.constants as const

path = '/home/chris/Desktop/WakeShape_LinGrad/'
path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/LinearGradient/'
#offset_arr = np.arange(-200,201,5)
offset_arr = np.arange(-120,100.5,2) #SepGrad Full

"""
path = '/home/chris/Desktop/NERSC_Jan_Grad/'
offset_arr = np.arange(-80,220.5,5) #JanGrad Full
tranExtent = 400 #Jan
npcase = 1e16 #jan
ind = 8
"""

#path = '/home/chris/Desktop/NERSC_Sep_Control2/'
#ind = 5
path = '/home/chris/Desktop/NERSC_Sep_Grad/'
ind = 9
#path = '/home/chris/Desktop/NERSC_Mar_Grad0.001/'
#path = '/home/chris/Desktop/NERSC_Mar_Grad0.01/'
#ind = 10
#path = '/home/chris/Desktop/NERSC_Dec_Grad/'
#ind = 9
offset_arr = np.arange(-120,100.5,2)
tranExtent = 200
npcase = 2e16

"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/NERSC_LIN_Aug/'
offset_arr = np.arange(-74,40.5,2)
npcase = 3e17#
tranExtent = 75
ind = 7
"""
"""
path = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/LowRes_HighGradient_August/LinearGradient/'
offset_arr = np.arange(-60,50.5,2)
npcase = 3e17#
tranExtent = 75
ind = 5
"""
"""
path = '/home/chris/Desktop/VELA_Oct_Grad/'
npcase = 2e16#
tranExtent = 200
offset_arr = np.arange(-120,100.5,2)
ind = 5
"""



#offset_arr = np.arange(-80,90.5,2)
rmajor_arr = np.zeros(len(offset_arr))
rminor_arr = np.zeros(len(offset_arr))
n0sh_arr = np.zeros(len(offset_arr))
n0sh2_arr = np.zeros(len(offset_arr))
dndysh_arr = np.zeros(len(offset_arr))

for i in range(len(offset_arr)):
    vmax = 1e18*(npcase/3.7e17)
    vmin = 1e16*(npcase/3.7e17)
    central_off = offset_arr[i]
    
    
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
              'plot' : False
              }
    theta,r = plot.wake_cross_section(params)
    #r, theta, rhoPol = plot.wake_cross_section(params)
    
    # Define model function to be used to fit to the data above:
    def coslike(x, *p):
        rp, b, re = p
        return rp-re*np.cos(x+b)
            
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 0., 1.]
            
    xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
    rmajor_arr[i] = xcoeff[0]
    rminor_arr[i] = xcoeff[2]
    
    if xcoeff[0] > 0: #Otherwise there is no wake
        radmax = xcoeff[0]*1.05
    
        params_sh = {'vmin' : vmin,
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
        x,y,n = plot.wake_cross_section_sheath(params_sh)
        sfit = np.polyfit(y,n,1)
        n0sh_arr[i] = sfit[1]
        dndysh_arr[i] = sfit[0]*1e4

dx=1.2#0.25
start = 97 #um

plt.plot(offset_arr*dx*-1+start,rmajor_arr)
plt.title("Wake radius in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake radius "+r'$(\mu m)$')
plt.grid(); plt.show();

trunc = 15    #15 for sep grad,
x = offset_arr[trunc:]*dx*-1+47
y = rminor_arr[trunc:]
h = np.polyfit(x,y,1)

#plt.plot(offset_arr*dx*-1+97,rminor_arr2, label = 'Sim2')
plt.plot(offset_arr[trunc:]*dx*-1+start,rminor_arr[trunc:], label = 'Sim')
plt.plot(offset_arr[trunc:]*dx*-1+start,x*h[0]+h[1],ls='--', label = 'Fit')
#plt.plot(offset_arr*dx*-1+97,x*0.039+h[1],ls='--', label = 'Guess')
plt.title("Wake vertical offset in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake vertical offset "+r'$(\mu m)$')
plt.legend(); plt.grid(); plt.show();

plt.plot(offset_arr*dx*-1+start,rmajor_arr-rminor_arr)
plt.plot(offset_arr*dx*-1+start,rmajor_arr+rminor_arr)
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2-rminor_arr2)
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2+rminor_arr2)
plt.title("Wake radius in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake radius "+r'$(\mu m)$')
plt.grid(); plt.show();

#This doesnt really work to fix it.  Idk
n0sh2_arr = n0sh_arr + dndysh_arr*rminor_arr/1e4

plt.plot(offset_arr[trunc:]*dx*-1+start,n0sh_arr[trunc:],label="Simulation")
plt.plot([-20,200],[7e16,7e16],ls='--',label="Value at "+r'$\xi_0$')
plt.plot(offset_arr[trunc:]*dx*-1+start,n0sh2_arr[trunc:],ls='dotted',label="Incl. Offset")
plt.title("Sheath's central density longitudinally")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Sheath's central density "+r'$(cm^{-3})$')
plt.grid(); plt.legend(); plt.show();

plt.plot(offset_arr[trunc:]*dx*-1+start,dndysh_arr[trunc:],label="Simulation")
plt.plot([-20,200],[1.89e18,1.89e18],ls='--',label="Value at "+r'$\xi_0$')
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2-rminor_arr2)
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2+rminor_arr2)
plt.title("Sheath density gradient longitudinally")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Sheath density gradient "+r'$(cm^{-4})$')
plt.grid(); plt.legend();  plt.show();

pi = np.pi
e = const.physical_constants['elementary charge'][0]
dndy_ion = 8e17

ering_arr = pi*e*rmajor_arr*rmajor_arr*0.1*2*dndysh_arr*1e-4**2
eion_arr =  pi*e*rmajor_arr*rmajor_arr*dndy_ion*1e-4**2
eoff_arr = 2*pi*e*npcase*rminor_arr*1e-4+0.5*pi*e*dndy_ion*3*rminor_arr*rminor_arr*1e-4**2

plt.plot(offset_arr[trunc:]*dx*-1+start,ering_arr[trunc:],label="Sheath Model")
plt.plot(offset_arr[trunc:]*dx*-1+start,eion_arr[trunc:],label="Ion Model * -1")
plt.plot(offset_arr[trunc:]*dx*-1+start,eoff_arr[trunc:],label="Offset Contribution")
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2-rminor_arr2)
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2+rminor_arr2)
plt.title("Longitudinal Variation of Sheath Deflection")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Electric field "+r'$(cgs)$')
plt.grid(); plt.legend();  plt.show();

plt.plot(offset_arr[trunc:]*dx*-1+start,ering_arr[trunc:]+eoff_arr[trunc:]-eion_arr[trunc:],c='k',label="Net Kick")
plt.plot(offset_arr[trunc:]*dx*-1+start,ering_arr[trunc:],label="Sheath Model")
plt.plot(offset_arr[trunc:]*dx*-1+start,-1*eion_arr[trunc:],label="Ion Model")
plt.plot(offset_arr[trunc:]*dx*-1+start,eoff_arr[trunc:],label="Offset Contribution")
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2-rminor_arr2)
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2+rminor_arr2)
plt.title("Longitudinal Variation of Sheath Deflection")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Electric field "+r'$(cgs)$')
plt.grid(); plt.legend();  plt.show();