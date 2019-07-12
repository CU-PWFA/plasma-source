#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:17:46 2019

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

case = 11
if case == 1:#Ideal Matching
    tpl_l = 0
    tpl_offset = 0
    tpl_n = 0 
    betastar = 0.04214 #m
    #waistloc = 0.0357 #m
if case == 10:#Best case with 5cm
    tpl_l = 0
    tpl_offset = 0
    tpl_n = 0 
    betastar = 0.05
if case == 11:#Best case with 10cm
    tpl_l = 0
    tpl_offset = 0
    tpl_n = 0 
    betastar = 0.10

tpl_foc_dist = 0.6684
if case == 2:#10cm with lens far upstream
    tpl_l = 55.1273022719
    tpl_n = 0.1
    d = -1.08197312858
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.1
    #waistloc = tpl_offset - d
if case == 3:#10cm with lens far downstream
    tpl_l = 257.700528841
    tpl_n = 0.1
    d = 1.08197312858
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.1
    #waistloc = tpl_offset - d

argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = 0., scaledown = 1)

leftext = 0.2#1 #m
ovenext = 0.5
dz = 1 #um
zarr = np.arange((-ovenext-leftext)*1e6, ovenext*1e6+dz, dz)
argon_params['Nz'] = len(zarr)
argon_params['Z'] = zarr[-1]-zarr[0]

lioven = profile.lithium_oven_profile(z=zarr, center=0, ne0=0.34)
start = lioven[0]/1e6
nfunc = lioven[2]
narr=nfunc(zarr)

if debug == 1:
    plt.plot(zarr,narr)
    plt.grid(); plt.show()
    
relset = ovenext + leftext + start
tpl_offset = tpl_offset + relset

argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, narr, tpl_offset*1e6, tpl_n, tpl_l, debug)
n = argon.nez
endindex = np.argmax(n)
zarr = zarr - zarr[0]
z = zarr/1e6

waist_arr = np.linspace(-0.01, 0.07, 101)
bmag_arr = np.zeros(len(waist_arr))
for i in range(len(waist_arr)):
    waistloc = waist_arr[i] + relset
    #Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                           beta_offset=waistloc, plasma_start=0)

    bmag_arr[i] = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
plt.plot(waist_arr,bmag_arr)
plt.show()
print(min(bmag_arr))