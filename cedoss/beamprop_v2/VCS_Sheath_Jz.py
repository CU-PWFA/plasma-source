#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:44:29 2022

@author: chris
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 11:04:17 2021

Wanted to plot the sheath density vs theta

This for fig 2 in paper 2

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

#path = '/home/chris/Desktop/NERSC_Sep_Grad/'
#ind = 9
#central_off = 0

path = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/JzDump/'

ind=5
central_off = 0#-20
npcase = 2e16#
tranExtent = 200

#82 for 2e16, 122 for 1e16, 37 for 3e17
radmax = 78#40
yoff = 10

e = const.physical_constants['elementary charge'][0]
c = const.physical_constants['speed of light in vacuum'][0]

params = {'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : False,
          'radmax' : radmax,
          'yoff' : yoff,
          'drive' : 'rhoBeam',
          'npcase' : npcase/1e17
          }

x,y,n = plot.wake_cross_section_sheath(params)
y2, nvy = plot.wake_cross_section_densityslice(params)
#r, theta, rhoPol = plot.wake_cross_section(params)
sfit = np.polyfit(y,n,1)
yfull = np.linspace(min(y),max(y),20)
print("n_s0 = ",sfit[1] , " (cm^-3)")
print("dn/dy_s = ",sfit[0]*1e4 , " (cm^-4)")
ion = [8e17,2e16]
y=np.flip(y,0)*-1
n=np.flip(n,0)
nvy = np.flip(nvy,0)
ion[0]=ion[0]*-1
sfit[0]=sfit[0]*-1
scale = npcase#1e17

#Now copy from my Jz cross section code:

#August 2021 n=2e16 runs
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
path = superpath + 'JzDump/'
grad = 0
yoff = 0 #m
radius = 76.149e-6 #m
ind = 5
tranExtent = 200        #tranextent for the sim
npcase = 2e16           #central density for the sim
dx=1.2                  #dx for the sim
#central_off = -20
simname = 'MatchedBeams'
Jz_plasma = 'JPlasma'
rho_plasma = 'rhoPlasma'

grad = 8e17
yoff = 3.962e-6 #m
radius = 76.435e-6 #m

c = const.physical_constants['speed of light in vacuum'][0]
e = const.physical_constants['elementary charge'][0]

vector = 0 # 0 for Ez, 1 for Ey, 2 for Ex
flip = True
c = const.speed_of_light

params = {'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : simname,#'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : False,
          'vector' : vector,
          'field' : Jz_plasma
          }
xa, ya, evx, evy, JzYZ = plot.wake_cross_section_field(params)
params['field'] = rho_plasma
params['vector'] = 0
xb, yb, bvx, bvy, rhoYZ = plot.wake_cross_section_field(params)

sumDen = -(rhoYZ-JzYZ/c) / (e*npcase*100**3)


den_v_x = np.flip(sumDen[int(len(evx)*1/4):int(len(evx)*3/4),int(len(evx)/2)],0)
axis1 = np.where(den_v_x==0)[0]+int(len(evx)*1/4)
nvals1 = np.zeros(len(axis1)*2)
yvals1 = np.zeros(len(axis1)*2)
for i in range(len(axis1)):
    den_v_y = np.flip(sumDen[:,axis1[i]],0)
    left_ind = np.argmax(den_v_y[:int(len(den_v_y)/2)])
    right_ind = np.argmax(den_v_y[int(len(den_v_y)/2):])+int(len(den_v_y)/2)
    nvals1[2*i] = den_v_y[left_ind]
    nvals1[2*i+1] = den_v_y[right_ind]
    yvals1[2*i] = ya[left_ind]
    yvals1[2*i+1] = ya[right_ind]
    
den_v_x = np.flip(sumDen[int(len(evx)/2),int(len(evx)*1/4):int(len(evx)*3/4)],0)
axis2 = np.where(den_v_x==0)[0]+int(len(evx)*1/4)
nvals2 = np.zeros(len(axis2)*2)
yvals2 = np.zeros(len(axis2)*2)
for i in range(len(axis2)):
    den_v_y = np.flip(sumDen[axis2[i],:],0)
    left_ind = np.argmax(den_v_y[:int(len(den_v_y)/2)])
    right_ind = np.argmax(den_v_y[int(len(den_v_y)/2):])+int(len(den_v_y)/2)
    nvals2[2*i] = den_v_y[left_ind]
    nvals2[2*i+1] = den_v_y[right_ind]
    yvals2[2*i] = -ya[axis2[i]]
    yvals2[2*i+1] = -ya[axis2[i]]

nvals = np.append(nvals1,nvals2)
yvals = np.append(yvals1,yvals2)

sfit_jz = np.polyfit(yvals,nvals,1)


plt.plot(ya,-np.flip(rhoYZ[:,int(len(evx)/2)],0),c='b',label='-rho')
plt.plot(ya,np.flip(JzYZ[:,int(len(evx)/2)],0)/c,c='r',label='Jz/c')
plt.legend();plt.show()

fig = plt.figure(figsize=(7,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.scatter(y,n/scale,s=10,marker='o', facecolors='none', edgecolors='grey',label='Sheath Density')
ax.scatter(yvals,nvals,s=10,marker='o', facecolors='none', edgecolors='orange',label='Jz spots')
ax.plot(y2,nvy/scale,c='k',label='Electrons')
ax.plot(ya,np.flip(sumDen[:,int(len(evx)/2)],0),c='g',ls='--',label='rho-Jz/c')
ax.plot(y2,(y2*ion[0]/1e4+ion[1])/scale,c='b',label='Ions')
ax.plot(y2,(y2*sfit[0]+sfit[1])/scale,c='r',label='Sheath Fit')
ax.plot(y2,(y2*sfit_jz[0]+sfit_jz[1]),c='r',ls='--',label='rho-Jz/c Fit')
#ax.plot(y2,(y2*sfit[0]+sfit[1])/scale+1.1,c='r',ls='--',label='Fit Offset')
#ax.set_ylabel("Density "+r'$\mathrm{(10^{17}\times cm^{-3})}$')
ax.set_ylabel("Density "+r'$\mathrm{(en_p)}$')
ax.set_xlabel("y "+r'$\mathrm{(\mu m)}$')
ax.set_ylim([-1e15/scale,1.5*max(n)/scale])
ax.legend()

plt.show()

print("Factor to multiply by:")
print(sfit_jz[0]/(sfit[0]/scale))

fig = plt.figure(figsize=(7,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.scatter(y,n/scale,s=10,marker='o', facecolors='none', edgecolors='grey',label='Sheath Density')
ax.scatter(yvals,nvals,s=10,marker='o', facecolors='none', edgecolors='orange',label='Jz spots')
ax.plot(y2,nvy/scale,c='k',label='Electrons')
ax.plot(ya,np.flip(sumDen[:,int(len(evx)/2)],0),c='g',ls='--',label='rho-Jz/c')
ax.plot(y2,(y2*ion[0]/1e4+ion[1])/scale,c='b',label='Ions')
ax.plot(y2,(y2*sfit[0]+sfit[1])/scale,c='r',label='Sheath Fit')
ax.plot(y2,(y2*sfit_jz[0]+sfit_jz[1]),c='r',ls='--',label='rho-Jz/c Fit')
#ax.plot(y2,(y2*sfit[0]+sfit[1])/scale+1.1,c='r',ls='--',label='Fit Offset')
#ax.set_ylabel("Density "+r'$\mathrm{(10^{17}\times cm^{-3})}$')
ax.set_ylabel("Density "+r'$\mathrm{(en_p)}$')
ax.set_xlabel("y "+r'$\mathrm{(\mu m)}$')
ax.set_ylim([-1e15/scale,1.5*max(n)/scale])
ax.legend()

plt.show()

"""
den_v_y = np.flip(sumDen[:,int(len(evx)/2)],0)
left_ind = np.argmax(den_v_y[:int(len(den_v_y)/2)])
right_ind = np.argmax(den_v_y[int(len(den_v_y)/2):])+int(len(den_v_y)/2)
slope = (den_v_y[left_ind]-den_v_y[right_ind])/(ya[left_ind]-ya[right_ind])
"""






