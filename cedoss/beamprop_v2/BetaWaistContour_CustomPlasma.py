#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:07:02 2019

Given a plasma of some function, make a contour of the matching beta waist and
waist location.  I use this as a crude way of finding the matching condition.

Specifically this is for the lithium oven, if I find other use for this code
I'll probably just copy it over.

@author: chris
"""
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import sys
sys.path.insert(0, "../../python")
from lens import profile

debug = 0
path = '/home/chris/Desktop/BeamProp/testGaussian'

argon_params = PProp.ReturnDefaultPlasmaParams(path, scaledown = 1)

leftext = 0#1 #m
ovenext = 0.5

zarr = np.linspace((-ovenext-leftext)*1e6, ovenext*1e6, argon_params['Nz'])
lioven = profile.lithium_oven_profile(z=zarr, center=0, ne0=0.34)
start = lioven[0]/1e6
nfunc = lioven[2]
narr=nfunc(zarr)

if debug == 1:
    plt.plot(zarr,narr)
    plt.grid(); plt.show()
    
argon = PProp.CustomPlasma(argon_params, narr, debug)
n = argon.nez
zarr = zarr - zarr[0]
z = zarr/1e6

relset = ovenext + leftext + start
beta_arr = np.linspace(0.041, 0.043, 101)
waist_arr = np.linspace(0.03, 0.045, 101) + relset
bmag_image = np.zeros((len(beta_arr),len(waist_arr)))
for i in range(len(beta_arr)):
    if i%10 == 0: print(i/len(beta_arr)*100,"%");
    for j in range(len(waist_arr)):
        betastar = beta_arr[i]
        waist_loc = waist_arr[j]
        #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=0)
        endindex = np.argmax(n)
        bmag = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #bmag = PProp.Calc_Bmag(beam_params,n, z)
        bmag_image[i][j] = bmag
        if debug == 1:
            betaf, alphaf, gammaf, gbf = PProp.Plot_CSEvo(beam_params, n, z, 0, legend_loc = 10)

#now I had to do a bunch of shit to get this to work, so the waist array needs to be interpreted as such
waist_arr = waist_arr - relset
minloc = PProp.PlotContour(bmag_image, beta_arr, waist_arr, r'$\beta^*$ [m]', r'$z_{v}$ [m]',simple=True)