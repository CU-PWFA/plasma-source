#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 14:06:19 2020

Loop over the cross section to track for all longitudinal pos:

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

dx=1.2#0.25
start = 97 #um

plt.plot(offset_arr*dx*-1+start,rmajor_arr)
plt.title("Wake radius in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake radius "+r'$(\mu m)$')
plt.grid(); plt.show();

trunc = 15
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

sys.exit()

def Stupakov(p, x):
    """
    if x*p[1] < np.pi / 2:
        return 2*np.sqrt(2*p[0]*np.sin(x*p[1]))
    else:
        return 2*np.sqrt(2*p[0])*np.cos(1/np.sqrt(2)*(x*p[1]-np.pi/2))
    """
    return np.piecewise(x, [x*p[1] < np.pi / 2, x*p[1] >= np.pi / 2], 
                        [lambda x: 2*np.sqrt(2*p[0]*np.sin(x*p[1])), 
                         lambda x: 2*np.sqrt(2*p[0])*np.cos(1/np.sqrt(2)*(x*p[1]-np.pi/2))])
    
kp = 266.123 #cm-1    
fix = 1.0#125

cut = 0#14
start = 4
rtop = np.flip((rmajor_arr[cut:-start])+(rminor_arr[cut:-start]),0)*1e-4*kp/fix
xaxs = np.flip(offset_arr[cut:-start]*dx*-1-min(offset_arr[cut:-start]),0)*1e-4*kp/fix
xaxs = xaxs - xaxs[0]

errfunc = lambda p, x, y: Stupakov(p, x) - y
p0 = [1e3,1e-2]
p1, success = optimize.leastsq(errfunc,p0[:], args=(xaxs, rtop))

plt.plot(xaxs,rtop,label="VSim Simulation")
#plt.plot(xaxs,Stupakov([1e3,1e-2], xaxs))
plt.plot(xaxs,Stupakov(p1, xaxs),label="Fitted Theory")
plt.title("Wake radius (top) compared to Stupakov theory")
plt.xlabel("Distance behind drive beam "+r'$(s*k_{p,0})$')
plt.ylabel("Wake radius "+r'$(r*k_{p,0})$')
plt.grid(); plt.legend(); plt.show()

cut = 5
start = 4
rtop = np.flip((rmajor_arr[cut:-start])-(rminor_arr[cut:-start]),0)*1e-4*kp*fix
rmin = np.flip((rminor_arr[cut:-start]),0)*1e-4*kp*fix
xaxs = np.flip(offset_arr[cut:-start]*dx*-1-min(offset_arr[cut:-start]),0)*1e-4*kp*fix
xaxs = xaxs - xaxs[0]

errfunc = lambda p, x, y: Stupakov(p, x) - y
p2, success = optimize.leastsq(errfunc,p0[:], args=(xaxs, rtop))

plt.plot(xaxs,rtop,label="VSim Simulation")
#plt.plot(xaxs,Stupakov([1e3,1e-2], xaxs))
plt.plot(xaxs,Stupakov(p2, xaxs),label="Fitted Theory")
plt.title("Wake radius (bot) compared to Stupakov theory")
plt.xlabel("Distance behind drive beam "+r'$(s*k_{p,0})$')
plt.ylabel("Wake radius "+r'$(r*k_{p,0})$')
plt.grid(); plt.legend(); plt.show()

print(p1)
print(p2)

plt.plot(xaxs,rtop,label="VSim Simulation")
plt.title("Integration on Bottom Trajectory (uniform)")
plt.xlabel("Distance behind drive beam "+r'$(s*k_{p,0})$')
plt.ylabel("Wake radius "+r'$(r*k_{p,0})$')


result = np.trapz(rtop, xaxs)
print(result)

trh = 0.807
result = np.trapz(rtop[np.where(rtop>trh)]-trh, xaxs[np.where(rtop>trh)])
print(result)

plt.plot([min(xaxs),max(xaxs)],[trh,trh],ls='--',label = "Half-Volume (r= 0.807)")
plt.grid(); plt.legend(); plt.show()

t1 = 4
t2 = 45
plt.plot(xaxs[t1:t2],rtop[t1:t2],label="VSim Simulation")
plt.grid(); plt.legend(); plt.show()
"""
#plt.plot(offset_arr*dx*-1+97,rminor_arr2, label = 'Sim2')
plt.plot(xaxs,rmin, label = 'Sim')
plt.plot(xaxs,(Stupakov(p1, xaxs)-Stupakov(p2, xaxs))/2, label = "Top Fit - Bot Fit")
#plt.plot(offset_arr*dx*-1+97,x*h[0]+h[1],ls='--', label = 'Fit')
#plt.plot(offset_arr*dx*-1+97,x*0.039+h[1],ls='--', label = 'Guess')
plt.title("Wake vertical offset in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(s*k_{p,0})$')
plt.ylabel("Wake vertical offset "+r'$(r*k_{p,0})$')
plt.legend(); plt.grid(); plt.show();

# Get the fitted curve
hist_fit = coslike(theta, *xcoeff)

plt.plot(theta,r)
plt.plot(theta,hist_fit)
plt.title("Wake Shape in Polar Coordinates")
plt.xlabel("Angle (rad)")
plt.ylabel("Radius of Wake Edge (um)")
plt.grid(); plt.show()
"""