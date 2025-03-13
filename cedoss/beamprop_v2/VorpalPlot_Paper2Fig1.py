#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 12:33:54 2022

To get all the good stuff, get 'electron' 'rhoWitness' 'rhoDrive for all the dumps
then just the final 'WitnessBeam' and 'DriveBeam' dump.

This version of vorpalplot is only for fig 1 in papaer 2.  Also for crazy shit with wake circles.

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def alpha_colormap(cmap, cutoff, flip=True):
    N = cmap.N
    cmapt = cmap(np.arange(N))
    alpha = np.ones(N)
    if flip:
        temp = alpha[:int(cutoff*N)]
        M = len(temp)
        alpha[:int(cutoff*N)] = np.linspace(0, 1, M)
    else:
        alpha[int((1-cutoff)*N):] = 0.0
    cmapt[:, -1] = alpha
    cmapt = colors.ListedColormap(cmapt)
    return cmapt

superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'

path =  '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g8e17/'

#August 2021 n=2e16 runs
#superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
superpath = '/media/chris/New Volume/VSimRuns/Oct2022_FinalSims/'
path = superpath + 'NERSC_n2e16_g0/'
grad = 0
yoff = 0 #m
radius = 76.149e-6 #m
ind = 5
tranExtent = 200        #tranextent for the sim
threshold = 100
npcase = 2e16           #central density for the sim
start=87
dx=1.2                  #dx for the sim
#central_off = -116#-20 for 114 um, -50 for 150 um, -92 for 200 um, -110 for 222 um
simname = 'MatchedBeams'
efield = 'edgeE'
bfield = 'faceB'
setno = 1
if setno == 1:
    path = superpath + 'NERSC_n2e16_g8e17/'
    grad = 8e17
    yoff = 3.962e-6 #m
    radius = 76.435e-6 #m
if setno == 10:
    path = superpath + 'NERSC_Deflection_Aug/'
    grad = 8e17
    yoff = 3.962e-6 #m
    radius = 76.435e-6 #m
    central_off = -50
    ind = 6
    simname = 'PTPL_Gradient'
    efield = 'ElecFieldPlasma'
    bfield = 'MagFieldPlasma'
elif setno == 2:
    path = superpath + 'NERSC_n2e16_g2e17/'
    grad = 2e17
    yoff = 0.946e-6 #m
    radius = 76.143e-6 #m
elif setno == 3:
    path = superpath + 'NERSC_n2e16_g2.5e16/'
    grad = 2.5e16
    yoff = 0.106e-6 #m
    radius = 76.122e-6 #m
elif setno == 4:
    path = superpath + 'NERSC_n2e16_g2.5e15/'
    grad = 2.5e15
    yoff = 0.0167e-6 #m
    radius = 76.154e-6 #m

vmax = 1e18*(npcase/3.7e17)#*1.8
vmin = 1e16*(npcase/3.7e17)#*90

params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'path' : path,
          'dumpInd' : 5,
          #'simName' : 'ThinPlasmaLens3D',
          'simName' : 'MatchedBeams',
          #'simName' : 'PTPL_Gradient',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'vmax' : 1e18,
          'vmin' : 1e16
          }
rhoBXY, rhoXY, xleft, xright, ybot, ytop = plot.drive_witness_density_paperplot(params)


### Next, grab the circles from vorpal cross section 

def coslike(x, *p):
    rp, b, re = p
    return rp-re*np.cos(x+b)
    
def getCircle(central_off):
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
              'drive' : 'rhoBeam',
              'threshold' : threshold #how far out in r to consider wake sheath
              }
    
    theta,r = plot.wake_cross_section(params)
    #r, theta, rhoPol = plot.wake_cross_section(params)
    
    # Define model function to be used to fit to the data above:

            
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 0., 1.]
            
    xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
    #xcoeff = np.array([  3.97932549e+01,  -4.01443535e-01,   2.79529067e-02])
    # Get the fitted curve
    
    x = r*np.sin(theta+xcoeff[1])
    x=np.append(x,x[0])
    y = -r*np.cos(theta+xcoeff[1])
    y=np.append(y,y[0])
    
    th_rad = xcoeff[0]
    th_off = xcoeff[2]
    
    x_th = th_rad*np.sin(theta+xcoeff[1])
    x_th=np.append(x_th,x_th[0])
    y_th = -th_rad*np.cos(theta+xcoeff[1])
    y_th=np.append(y_th,y_th[0])+th_off
    
    return x, y, x_th, y_th, th_off

#offset_list = np.array([-20,-50,-92,-116])
#colors_list = np.array(['green','yellow','orange','red'])
offset_list = np.array([-11,-95,-118])#-21.-94,-118
colors_list2 = np.array(['green','blue','red'])
colors_list = np.array(['lightgreen','cornflowerblue','tomato'])
x1, y1, xt1, yt1, yoff_1 = getCircle(offset_list[0])
x2, y2, xt2, yt2, yoff_2 = getCircle(offset_list[1])
x3, y3, xt3, yt3, yoff_3 = getCircle(offset_list[2])
#x4, y4, xt4, yt4 = getCircle(offset_list[3])

### Now combine these two!!

fig, (ax0,ax2) = plt.subplots(nrows=2,ncols=1)#,sharex=True)
fig.set_size_inches(5,7)

#fig = plt.figure(figsize=(5, 3.75), dpi=100, frameon=False)
#ax = fig.add_subplot(111)
#ax = plt.Axes(fig, [0., 0., 1., 1.])
#ax.set_axis_off()
#fig.add_axes(ax)
    
    # Plot the drive beam
ax0.imshow(rhoBXY, interpolation='gaussian', extent=[xleft, xright, ybot, ytop], cmap='copper')

cmapP = alpha_colormap(plt.cm.get_cmap('inferno'), 0.2, True)
cf = ax0.imshow(rhoXY, interpolation='gaussian', aspect='auto', extent=[xleft, xright, ybot, ytop],
               norm=colors.LogNorm(vmin=1e15, vmax=2e18), cmap=cmapP)#1e16,2e18
fig.colorbar(cf, ax=ax0,label=r'$\mathrm{Plasma \ Electron \ Number \ Density \ (cm^{-3})}$')
for i in range(len(offset_list)):
    ax0.vlines(x = offset_list[i]*dx*-1+start, ymin = -200, ymax = 200, linestyle = '--', colors = colors_list2[i])
#plt.vlines(x = 114, ymin = -200, ymax = 200, linestyle = '--', colors = 'green',label = '1')
#plt.vlines(x = 150, ymin = -200, ymax = 200, linestyle = '--', colors = 'yellow',label = '2')
#plt.vlines(x = 200, ymin = -200, ymax = 200, linestyle = '--', colors = 'red',label = '3')
ax0.set_xlabel(r'$\xi \ \mathrm{(\mu m)}$')
ax0.set_ylabel(r'$y \ \mathrm{(\mu m)}$')
ax0.text(260,175,'(a)',color='white')
    
#plt.savefig(path+'Title_Wake.png')

#ax2 = fig.add_subplot(121)
ax2.plot(x1,y1,c=colors_list[0])
ax2.plot(xt1,yt1,c=colors_list2[0],linestyle='--',label=r'$\xi=$'+str(offset_list[0]*dx*-1+start)+r'$\mathrm{\ \mu m}$')
ax2.plot([-100,180],[yoff_1,yoff_1],c=colors_list2[0],linestyle='dotted')
ax2.plot(x2,y2,c=colors_list[1])
ax2.plot(xt2,yt2,c=colors_list2[1],linestyle='--',label=r'$\xi=$'+str(offset_list[1]*dx*-1+start)+r'$\mathrm{\ \mu m}$')
ax2.plot([-100,180],[yoff_2,yoff_2],c=colors_list2[1],linestyle='dotted')
ax2.plot(x3,y3,c=colors_list[2])
ax2.plot(xt3,yt3,c=colors_list2[2],linestyle='--',label=r'$\xi=$'+str(offset_list[2]*dx*-1+start)+r'$\mathrm{\ \mu m}$')
ax2.plot([-100,180],[yoff_3,yoff_3],c=colors_list2[2],linestyle='dotted')
#ax.plot(x4,y4,c=colors_list[3])
#ax.plot(xt4,yt4,c=colors_list[3],linestyle='--')
ax2.hlines(y = 0, xmin = -200, xmax = 200, linestyle = '--', colors = 'black')
ax2.set_aspect('equal')
ax2.set_xlim([-100,180])
ax2.set_xlabel(r'$x \ \mathrm{(\mu m)}$')
ax2.set_ylabel(r'$y \ \mathrm{(\mu m)}$')
ax2.legend(loc=1)
ax2.text(-97,77,'(b)',color='black')
plt.savefig("/home/chris/Desktop/figs/fig1.eps",bbox_inches='tight')
plt.show()
