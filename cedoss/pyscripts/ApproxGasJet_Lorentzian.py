#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:58:10 2019

Since Paraview cannot take clips of rotationally extruded data (meaning I
cant use Wedge simulations to directly export to Fourier propagation) this
script uses the radially Gaussian and axially exponential fits to approximate
a gas jet profile.

This version uses Lorentzian profiles radially, and takes advantage of the
approximate quadratic fit axially for this small domain.

@author: chris
"""

import sys
import os
import numpy as np

sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

save_data = 1
reducer = 0

n0 = 1e17/1e17

folder = '/home/chris/Desktop/DataLoads/DensityFilesNp/'
directory = 'gasjet_Ar3e16_60x24x400_Lor/'

xsize = 20e2; ysize = 20e2; zsize = 20e4#4e4#32e4 #um
nx = 2**(9-reducer)
ny = 2**(9-reducer)
nz = 2**(8-reducer)

x = np.linspace(-xsize/2, xsize/2, nx)
y = np.linspace(-ysize/2, ysize/2, ny)
z = np.linspace(-zsize/2, zsize/2, nz)

xaxis = np.reshape(x, (len(x), 1, 1))
yaxis = np.reshape(y, (1, len(y), 1))
zaxis = np.reshape(z, (1, 1, len(z)))

na = 2.56260482473e+19
nb = -3.63035778292e+15
nc = 650914468889.0

ga = 4020.58922632
gb = 0.811406408508
gc = 4.48681341389e-05

y_n0 = na + nb*yaxis + nc*np.square(yaxis)
y_gam= ga + gb*yaxis + gc*np.square(yaxis)
approx = y_n0/(2*np.pi)*y_gam/(np.square(xaxis)+np.square(zaxis)+np.square(1/2*y_gam))
approx = approx / (1 + np.power(0.00015*np.sqrt(np.square(xaxis)+np.square(zaxis)),6))
f0 = approx[int(len(x)/2),int(len(y)/2),int(len(z)/2)]
approx = approx * (n0/f0)

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
    