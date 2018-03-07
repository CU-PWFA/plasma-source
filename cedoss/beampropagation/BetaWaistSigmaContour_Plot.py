#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:46:21 2018

plots the contour over beta and sigma

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp

betalabel = 10

path = '/home/chris/Desktop/DataLoads/ContourBetaWaistSigma_10GeV_HighRes_PostFix/'+str(betalabel)+'cm/'
print(path)
bmag_image = np.load(path+'bmagarr.npy')
sigma_arr = np.load(path+'sig.npy')
zbeta_arr = np.load(path+'beta.npy')
"""
for i in range(len(sigma_arr)):
    for j in range(len(beta_arr)):
        if bmag_image[i][j] < 1.0:
            bmag_image[i][j] = np.inf
        if bmag_image[i][j] > 1e3:
            bmag_image[i][j] = np.inf
        if not np.isfinite(bmag_image[i][j]):
            bmag_image[i][j] = 1e3
"""
minloc = PProp.PlotContour(bmag_image, zbeta_arr*100, sigma_arr*100, r'$ z_{\beta^*} [cm] $', r'$ \sigma_{hw} $ [cm]')