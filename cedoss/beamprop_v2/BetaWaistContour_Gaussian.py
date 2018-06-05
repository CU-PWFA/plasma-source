#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 10:57:47 2018

Finds optimal beta and waist for Gaussian ramps

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import sys

debug = 1
path = '/home/chris/Desktop/BeamProp/testGaussian'

sighw = 0.01 * 1e6
z0 = sighw/1e6*10
argon_params = PProp.ReturnDefaultPlasmaParams(path, sigma_hw = sighw, plasma_start = z0, scaledown = 100)
argon = PProp.GaussianRampPlasma(argon_params, debug)

n = argon.nez
z = np.linspace(0,argon_params['Z'],len(n))/1e6

dump = 10
cores = 4

ramp_end = argon_params['z0']/1e6
endindex = np.nonzero(z>ramp_end)[0][0]

beta_arr = np.linspace(0.008, 0.015, 200)
waist_arr = np.linspace(-0.02,-0.0075, 200)

bmag_image = np.zeros((len(beta_arr),len(waist_arr)))
for i in range(len(beta_arr)):
    if i%10 == 0: print(i/len(beta_arr)*100,"%");
    for j in range(len(waist_arr)):
        betastar = beta_arr[i]
        waist_loc = waist_arr[j]
        
        #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z0)
        
        bmag = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #bmag = PProp.Calc_Bmag(beam_params,n, z)
        bmag_image[i][j] = bmag

minloc = PProp.PlotContour(bmag_image, beta_arr, waist_arr, r'$\beta^*$ [m]', r'$z_{v}$ [m]')