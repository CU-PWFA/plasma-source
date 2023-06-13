#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 10:59:32 2021

Plotting some data from the 2bunch sims in Aug-Oct 2021

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
y0 = -1*np.array([1.96,1.28,0.01,-1.99,-5.05]) #Fitted
y0_c = -1*-0.02
ef = np.array([3.12,3.12,3.13,3.13,3.13])
ef_c = 3.11

#Averaged y0 is basically the same
#y0 = -1*np.array([1.91,1.26,-0.03,-2.05,-5.10]) #Average

#yp0 = -1*np.array([0.017,0.011,-0.0003,-0.018,-0.047]) #Average y'

flip = True

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
#Then run VorpalCrossSection_EfieldLooper with vec=1
# to get "e0_arr" for the electric field on axis

left = 25
right = 65
z_arr = (np.flip(offset_arr[left:right],0)*dx*-1+start)
radarr = np.flip(rcontr_arr[left:right],0)
ey_simey0 = np.flip(e0_arr[left:right],0)
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

#n2e16, g8e17
rp = 76.15e-6 #m  #MAXIMUM radius
n_cent0 = 2e16*100**3#*0.95 #m-3
slope = 8e17*100**4 #m-4
if flip:
    slope = slope*-1
    ey_simey0 = ey_simey0*-1
prop_dist = 0.1112 #m
ltpl = 300e-6
gamma = 19569.5
zsi_max = 112.2e-6 #m
sheath_slope_ave = 2.40e18*100**4#-8.64e17 * 100**4
#adj_fac = 2.4
#ring_adj = 883344572.88*0.94
isProd = True

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
isProd = False
"""

facA = 0.3112#0.3239
rpl = rp*np.sqrt((n_cent0-facA*rp*slope)/n_cent0) #sign flip b/c slope is wack
rmn = rp*np.sqrt((n_cent0+facA*rp*slope)/n_cent0)
yoff = zsi*(rpl-rmn)/(2*zsi_max)
#yoff = 4.2188366149411554e-6 #m

facB = 3.19089#4.059#2.482
#facB = 4.059

slope_cgs = slope/100**4
#slope_cgs_sh = -4.748e-37*slope_cgs**3  + 2.750e-18*slope_cgs**2 - 0.39*slope_cgs + 5.429e16
slope_cgs_sh = slope_cgs*facB
sheath_slope_emp = slope_cgs_sh * 100**4

facD = 0.090421#0.0651#0.134
#facD = 0.0651
L_sh = facD*rp

#ring_adj = (sheath_slope_ave)*(L_sh)*(radius_arr)*e/2/eps * adj_fac #To match the error
ring_emp = (sheath_slope_ave-slope)*(L_sh)*(radius_arr)*e/2/eps
ring_fullemp = (sheath_slope_emp-slope)*(L_sh)*(radius_arr)*e/2/eps
#ring=-8.833e8#V/m

n_cent = n_cent0 + slope*yoff

ybar = 1/4*np.power(radius_arr,2)*slope/n_cent
ey_emp = 1/(2*eps)*e*n_cent*(ysi-2*ybar-yoff) + 1/(8*eps)*e*slope*(3*np.power(ysi-yoff,2)) + ring_emp
#ey_adjring = 1/(2*eps)*e*n_cent*(ysi-2*ybar-yoff) + 1/(8*eps)*e*slope*(3*np.power(ysi-yoff,2)) + ring_adj
ey_noring = 1/(2*eps)*e*n_cent*(ysi-2*ybar-yoff) + 1/(8*eps)*e*slope*(3*np.power(ysi-yoff,2))
ey_fullemp = 1/(2*eps)*e*n_cent*(ysi-2*ybar-yoff) + 1/(8*eps)*e*slope*(3*np.power(ysi-yoff,2)) + ring_fullemp
ey_noyoff = 1/(2*eps)*e*n_cent*(-2*ybar)


angle_emp = -e*ey_emp*ltpl/(gamma*me*c**2)
angle_fullemp = -e*ey_fullemp*ltpl/(gamma*me*c**2)
#angle_adjring = -e*ey_adjring*ltpl/(gamma*me*c**2)
angle_noring = -e*ey_noring*ltpl/(gamma*me*c**2)
angle_simey0 = -e*ey_simey0*ltpl/(gamma*me*c**2)
angle_noyoff = -e*ey_noyoff*ltpl/(gamma*me*c**2)

deflection_emp = angle_emp * prop_dist 
deflection_fullemp = angle_fullemp * prop_dist 
#deflection_adjring = angle_adjring * prop_dist 
deflection_noring = angle_noring * prop_dist 
deflection_simey0 = angle_simey0 * prop_dist 
deflection_noyoff = angle_noyoff * prop_dist 
"""
plt.figure(figsize=(7,4))
plt.title("Offset of Witness Beam at Focus 11 cm Downstream")
plt.plot([zsi[0]*1e6,zsi[-1]*1e6],[y0_c,y0_c],c='k',ls='--',label="No Gradient")
#plt.plot(zsi*1e6,deflection_noyoff*1e6,c = 'gold',label = 'No E_ring, No y_off')
#plt.plot(zsi*1e6,deflection_noring*1e6,c = 'g',label = 'No E_ring')
plt.plot(zsi*1e6,deflection_fullemp*1e6,c = 'red',label = 'Using Empirical Model')
#plt.plot(zsi*1e6,deflection_emp*1e6,c = 'purple',label = 'Sim. Sheath, Emp. E_ring')
#plt.plot(zsi*1e6,deflection_adjring*1e6,c = 'orange',label = 'Sim. "E_ring"')
plt.plot(zsi*1e6,deflection_simey0*1e6,c = 'b',label = 'Using Single-Bunch Simulations')
if isProd:
    plt.scatter(separ,y0,label = "Two-Bunch Simulations")
plt.xlabel("Drive-Witness Separation " + r'$(\mu m)$')
plt.ylabel("y Centroid " + r'$(\mu m)$')
plt.legend(loc=3);plt.grid();plt.show()
"""
###Paper Figure

fig, ax1 = plt.subplots(figsize=(5,3))
ax1.plot([zsi[0]*1e6,zsi[-1]*1e6],[y0_c/prop_dist,y0_c/prop_dist],c='k',ls='--',label="No Gradient")
ax1.plot(zsi*1e6,deflection_fullemp*1e6/prop_dist,c = 'red',label = 'Empirical Model')
ax1.plot(zsi*1e6,deflection_simey0*1e6/prop_dist,c = 'b',label = 'Single-Bunch Simulation')
ax1.scatter(separ,y0/prop_dist,c='black',label = "Two-Bunch Simulations")
#ax1.scatter(separ,yp0*1e3,c='black',label = "Two-Bunch Simulations")
ax1.set_xlabel("Drive-Witness Separation " + r'$(\mu m)$')
ax1.set_ylabel("Deflection Angle "+r'$\mathrm{(\mu rad)}$')
#ax2=ax1.twinx()
#ax2.set_ylabel("y Centroid " + r'$(\mu m)$')
#ax2.plot([zsi[0]*1e6,zsi[-1]*1e6],[y0_c,y0_c],c='k',ls='--')
#ax2.plot(zsi*1e6,deflection_fullemp*1e6,c = 'red')
#ax2.plot(zsi*1e6,deflection_simey0*1e6,c = 'b')
#ax2.scatter(separ,y0,c='black')
#ax2.scatter(separ,yp0*1e3*prop_dist,c='black')
ax1.legend(loc=0)
plt.savefig("/home/chris/Desktop/figs/fig8.eps",bbox_inches='tight')
plt.show()
