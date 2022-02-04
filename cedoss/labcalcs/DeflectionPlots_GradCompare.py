#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 23:46:11 2021

Plotting some data from the 2bunch sims in Aug-Oct 2021

This version seeks to compare between different density gradients

@author: chris
"""

import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
from modules import ThreeDimensionAnalysis as ThrDim

separ = np.array([114,130,150,170,190])
zfoc = np.array([11.12,11.12,11.12,11.12,10.82])
zfoc_c = 11.12
sigx = np.array([3.95,3.97,3.99,3.98,4.07])
sigx_c = 3.93
x0 = np.array([0.02,0.05,-0.003,-0.006,0.02])
x0_c = 0.02
sigy = np.array([3.92,3.95,3.94,3.92,4.0])
sigy_c = np.array([3.95])
y0 = np.array([1.96,1.28,0.01,-1.99,-5.05])
y0_c = -0.02
ef = np.array([3.12,3.12,3.13,3.13,3.13])
ef_c = 3.11
"""
plt.scatter(separ,zfoc,label = "Two-Bunch Sims")
plt.plot([separ[0],separ[-1]],[zfoc_c,zfoc_c],c='k',ls='--',label="No Gradient")
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("Focal Distance " + r'$(cm)$')
plt.legend();plt.show()

plt.scatter(separ,sigx,label = "Two-Bunch Sims")
plt.plot([separ[0],separ[-1]],[sigx_c,sigx_c],c='k',ls='--',label="No Gradient")
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("x Spot size " + r'$(\mu m)$')
plt.legend();plt.show()

plt.scatter(separ,x0,label = "Two-Bunch Sims")
plt.plot([separ[0],separ[-1]],[x0_c,x0_c],c='k',ls='--',label="No Gradient")
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("x Centroid " + r'$(\mu m)$')
plt.legend();plt.show()

plt.scatter(separ,sigy,label = "Two-Bunch Sims")
plt.plot([separ[0],separ[-1]],[sigy_c,sigy_c],c='k',ls='--',label="No Gradient")
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("y Spot Size " + r'$(\mu m)$')
plt.legend();plt.show()

plt.scatter(separ,ef,label = "Two-Bunch Sims")
plt.plot([separ[0],separ[-1]],[ef_c,ef_c],c='k',ls='--',label="No Gradient")
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("Emmitance " + r'$(mm-mrad)$')
plt.legend();plt.show()
"""

#First I like to run VorpalCrossSection_Looper2
# to get "rcontr_arr" for the radius evolution
# and "offset_arr*dx*-1+start" for the zsi

left = 25
right = 65
z_arr = (np.flip(offset_arr[left:right],0)*dx*-1+start)
radarr = np.flip(rcontr_arr[left:right],0)
"""
print(z_arr) #to get between 114 and 200
sys.exit()
"""
zsi = z_arr*1e-6
radius_arr = radarr*1e-6  #EVOLVING radius

ysi = 0#y*1e-6
xsi = 0#x*1e-6
eps = 8.854e-12
pi = np.pi
e = const.physical_constants['elementary charge'][0]
me = const.physical_constants['electron mass'][0]
c = const.physical_constants['speed of light in vacuum'][0]
"""
#n2e16, g8e17
rp = 76.15e-6 #m  #MAXIMUM radius
n_cent = 2e16*100**3#*0.95 #m-3
slope = -8e17*100**4 #m-4
prop_dist = 0.1112 #m
ltpl = 300e-6
gamma = 2e4
zsi_max = 114e-6 #m
sheath_slope_ave = -8.64e17 * 100**4
#adj_fac = 2.4
ring_adj = -883344572.88*0.94
is2 = True
"""

#n2e16, g2.5e16
rp = 76.15e-6 #m  #MAXIMUM radius
n_cent = 2e16*100**3#*0.95 #m-3
slope = -2.5e16*100**4 #m-4
prop_dist = 0.1112 #m
ltpl = 300e-6
gamma = 2e4
zsi_max = 114e-6 #m
sheath_slope_ave = -6.01e16 * 100**4
ring_adj = -123211303.056
is2 = False


facA = 0.3239
rpl = rp*np.sqrt((n_cent-facA*rp*slope)/n_cent) #sign flip b/c slope is wack
rmn = rp*np.sqrt((n_cent+facA*rp*slope)/n_cent)
yoff = zsi*(rpl-rmn)/(2*zsi_max)

slope_cgs = -slope/100**4
slope_cgs_sh = -4.748e-37*slope_cgs**3  + 2.750e-18*slope_cgs**2 - 0.39*slope_cgs + 5.429e16
sheath_slope_emp = -slope_cgs_sh * 100**4

L_sh = 0.08*rp

ring_fullemp = (sheath_slope_emp)*(L_sh)*(radius_arr)*e/2/eps

ybar = -1/4*np.power(radius_arr,2)*slope/n_cent
ey_fullemp = 1/(2*eps)*e*n_cent*(ysi+2*ybar-yoff) + 1/(8*eps)*e*slope*(3*np.power(ysi-yoff,2)) + ring_fullemp

angle_fullemp = e*ey_fullemp*ltpl/(gamma*me*c**2)

if is2:
    angle_fullemp2 = angle_fullemp

plt.figure(figsize=(7,4))
plt.title("Deflection Angle for Two Different Cases")
plt.plot(zsi*1e6,angle_fullemp2*1e6,c = 'blue',label = 'High Gradient Case, g = 0.31')
plt.plot(zsi*1e6,angle_fullemp*1e6,c = 'red',label = 'Realistic Gradient Case, g = 0.01')
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("Vertical Deflection Angle " + r'$(\mu rad)$')
plt.legend(loc=3);plt.grid();plt.show()
