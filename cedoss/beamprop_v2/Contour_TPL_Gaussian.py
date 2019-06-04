#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 14:05:11 2018

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

debug = 0
path = '/home/chris/Desktop/BeamProp/testGaussian'
gamma = PProp.def_gamma

case = 20
if case == 1:
    sighw = 0.08 * 1e6
    tpl_n = 0.5
    tpl_l = 174.5
    zvac = -0.2006
    z0 = sighw/1e6*5
    betastar = .10
if case == 2:
    sighw = 0.01 * 1e6
    tpl_n = 0.5
    tpl_l = 607.5
    zvac = -0.01372
    z0 = sighw/1e6*6
    betastar = .10
    
if case == 20: #442.1um to match 5cm beta into a ramp that requires 2.5cm beta
    tpl_n = 0.5
    tpl_l = 442.1
    sighw = 0.0273 * 1e6
    zvac = -0.054
    betastar = 0.05
    z0 = sighw/1e6*6

argon_params = PProp.ReturnDefaultPlasmaParams(path, sigma_hw = sighw, plasma_start = z0, scaledown = 1)
argon = PProp.GaussianRampPlasma(argon_params, debug)

n_arr = argon.nez
z_arr = np.linspace(0,argon_params['Z'],len(n_arr))/1e6

ramp_end = argon_params['z0']/1e6
endindex = np.nonzero(z_arr>ramp_end)[0][0]

focal = Foc.Calc_Focus_Square_SI(tpl_n*1e17, tpl_l/1e6, gamma)

beta_f = Foc.Calc_BetaStar(betastar, focal)
tpl_f = focal*(1-beta_f/betastar)

waist_loc = zvac - tpl_f
"""For the full 1.01 contour line
offset_arr = np.linspace(-0.01, 0.01, 201)
length_arr = np.linspace(tpl_l - 70., tpl_l + 90, 201)
#"""
offset_arr = np.linspace(-0.01, 0.01, 201)
length_arr = np.linspace(tpl_l - 160., tpl_l + 180, 201)
bmag_image = np.zeros((len(offset_arr),len(length_arr)))

for i in range(len(offset_arr)):
    if i%10 == 0: print(i/len(offset_arr)*100,"%");
    for j in range(len(length_arr)):
        n = np.copy(n_arr); z = np.copy(z_arr)
        tpl_offset = waist_loc + offset_arr[i]
        tpl_l_set = length_arr[j]
        
        #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z0)
        argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n, tpl_offset*1e6, tpl_n, tpl_l_set, debug)
        
        bmag_image[i][j] = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #print("Bmag CS: ",bmagc)
minloc = PProp.PlotContour(bmag_image, offset_arr*100, length_arr, r'$\Delta z_{TPL}$ [cm]', r'$TPL_{\rm L}$ [um]', simple=True)