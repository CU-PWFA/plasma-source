#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:46:21 2018

plots the contour over beta and sigma

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp


bmag_image = np.load('/home/chris/Desktop/DataLoads/ContourBetaSigma/bmagarr.npy')
sigma_arr = np.load('/home/chris/Desktop/DataLoads/ContourBetaSigma/sig.npy')
beta_arr = np.load('/home/chris/Desktop/DataLoads/ContourBetaSigma/beta.npy')
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
minloc = PProp.PlotContour(bmag_image, sigma_arr*100, beta_arr*100, r'$ \sigma_{hw} $ [cm]', r'$ \beta^* [cm] $')