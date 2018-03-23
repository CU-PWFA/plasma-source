#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:46:21 2018

python can't import modules in a cell mid-script.  What?

'...testarr.npy' is sig = 0.1325 / zb* = -0.36
'...testarr2.npy'   sig = 0.16   / zb* = -0.4513
'...testarr3.npy'   sig = 0.12   / zb* = -0.3226

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp


bmag_image = np.load('/home/chris/Desktop/DataLoads/ContourDensityGamma/testarr.npy')
npl0_arr = np.load('/home/chris/Desktop/DataLoads/ContourDensityGamma/testx.npy')
gamma_arr = np.load('/home/chris/Desktop/DataLoads/ContourDensityGamma/testy.npy')

for i in range(len(npl0_arr)):
    for j in range(len(gamma_arr)):
        if bmag_image[i][j] < 1.0:
            bmag_image[i][j] = np.inf
        if bmag_image[i][j] > 1e3:
            bmag_image[i][j] = np.inf
        if not np.isfinite(bmag_image[i][j]):
            bmag_image[i][j] = 1e3

minloc = PProp.PlotContour(bmag_image, npl0_arr, gamma_arr, r'$ npl0 $ [cm^-3]', r'$ \gamma $', log=1)