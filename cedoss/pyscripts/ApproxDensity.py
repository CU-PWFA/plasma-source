#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 09:43:08 2017

Taking the fit parameters from the outputs of RefractionAnalysis and 
ThreeDimensionAnalysis, we create a 3D approximate density profile to
serve as a simple model for the thin plasma lens.  If save_data is set
to 1, then this script also saves a mock density and params file to
mimic Robert's propagation output.  This allows for the similar analysis
of this approximation to simulation data.

@author: chris
"""

import sys
import os
import numpy as np

sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

save_data = 1

#Parameters for the yz plane's elliptical tanh
a = 14.3764312172
b = 2.26796567463
n_0 = 5.73756199406
scl = 0.0660840157499

#Parameters for the x axis Gaussian
#n_0 from above
sig = 494.534987065
x_0 = -3.44452948425

folder = '/home/chris/Desktop/FourierPlots/ApproximateSol/'
directory = 'ETanhGauss1/'

pyz = [a, b, n_0, scl]
px  = [n_0, sig, x_0]

#Build our 3D approximate density
xwindow = 5000; xstep = 25 #microns
ywindow = 50  ; ystep = .25
zwindow = 200 ; zstep = 1

x = np.arange(-xwindow, xwindow, xstep)
y = np.arange(-ywindow, ywindow, ystep)
z = np.arange(-zwindow, zwindow, zstep)

xaxis = np.reshape(x, (len(x), 1, 1))
yaxis = np.reshape(y, (1, len(y), 1))
zaxis = np.reshape(z, (1, 1, len(z)))

x_Gauss  = ThrDim.Gaussian(px, xaxis)
yz_plane = ThrDim.DoubleTanh(pyz, ThrDim.EllipDist(yaxis, zaxis, scl))
approx = x_Gauss * yz_plane / n_0

#Plot our 3D approximate density
ThrDim.ImageCut(approx, x, y, z, 0, 0, 0, 1e-3,
                '(mm)','Plasma Density','e17(cm^-3)',1)

if save_data == 1:
    path = folder + directory
    if not os.path.exists(path):
        os.makedirs(path)
    
    approx = ThrDim.RobertRoll(ThrDim.RobertRoll(approx))
    np.save(path+'finalDensity.npy', approx)
    params = {'Nx' : int(2 * ywindow / ystep),
              'Ny' : int(2 * zwindow / zstep),
              'Nz' : int(2 * xwindow / xstep),
              'X' : 2 * ywindow,
              'Y' : 2 * zwindow,
              'Z' : 2 * xwindow
              }
    np.save(path+'params', params)
    