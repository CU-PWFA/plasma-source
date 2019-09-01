#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:03:06 2019

Since Paraview cannot take clips of rotationally extruded data (meaning I
cant use Wedge simulations to directly export to Fourier propagation) this
script uses the radially Gaussian and axially exponential fits to approximate
a gas jet profile.

@author: chris
"""

import sys
import os
import numpy as np

sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

save_data = 1
reducer = 0
"""This was with the WedgeRAS R2 simulation
a = 7.55200689789e+15 #cm-3
b = 0.000487307278485 #um-1
c = 1.48516032641e+14 #cm-3
f0 = a*np.exp(-b*5000)+c
n0 = 5e16/1e17
## sig(y) = sa + sb y + sc y**2     From ApproxGasJet_SigVSAxial.py
sa = 3536.31577167; sb = 0.553815660633; sc = -9.23072971528e-06
#"""
a = 3.31743175609e+16 #cm-3
b = 0.00049976694808  #um-1
c = 2.39194210599e+14 #cm-3
f0 = a*np.exp(-b*5000)+c
n0 = 3e16/1e17
## sig(y) = sa + sb y + sc y**2     From ApproxGasJet_SigVSAxial.py
sa = 2986.91262632; sb = 0.550205150854; sc = -5.09030500347e-06

folder = '/home/chris/Desktop/DataLoads/DensityFilesNp/'
directory = 'gasjet_Ar3e16_60x24x400/'

xsize = 60e2; ysize = 24e2; zsize = 4e4#32e4 #um
nx = 2**(9-reducer)
ny = 2**(9-reducer)
nz = 2**(8-reducer)

x = np.linspace(-xsize/2, xsize/2, nx)
y = np.linspace(-ysize/2, ysize/2, ny)
z = np.linspace(-zsize/2, zsize/2, nz)

xaxis = np.reshape(x, (len(x), 1, 1))
yaxis = np.reshape(y, (1, len(y), 1))
zaxis = np.reshape(z, (1, 1, len(z)))

y_expo = a*np.exp(-b*(yaxis+5000))+c
xz_gauss = np.exp(-(np.square(xaxis)+np.square(zaxis))/(2*np.square(sa+sb*yaxis+sc*np.square(yaxis))))
approx = y_expo * xz_gauss * (n0/f0)

#x_Gauss  = ThrDim.Gaussian(px, xaxis)
#yz_plane = ThrDim.DoubleTanh(pyz, ThrDim.EllipDist(yaxis, zaxis, scl))
#approx = x_Gauss * yz_plane / n_0

#Plot our 3D approximate density
ThrDim.ImageCut(approx, x, y, z, 0, 0, 0, 1e-3,
                r'$(mm)$','Gas Density (wrong axis labels)','e17(cm^-3)',1)

if save_data == 1:
    path = folder + directory
    if not os.path.exists(path):
        os.makedirs(path)
    
    np.save(path+'gasdensity.npy', approx)
    params = {'Nx' : nx,
              'Ny' : ny,
              'Nz' : nz,
              'X' : xsize,
              'Y' : ysize,
              'Z' : zsize
              }
    np.save(path+'params', params)
    