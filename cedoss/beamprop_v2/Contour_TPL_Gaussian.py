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

case = 30
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
if case == 30: #736.9um to match 5cm beta into a ramp that requires 2.5cm beta in 3e16
    tpl_n = 0.3
    tpl_l = 736.8569
    sighw = 0.025415 * 1e6
    zvac = -0.04546# -0.0246
    betastar = 0.05
    z0 = sighw/1e6*6
if case == 31: #736.9um to match 5cm beta into a ramp that requires 2.5cm beta in 3e16
    tpl_n = 0.3
    tpl_l = 744.187
    sighw = 0.025415 * 1e6
    zvac = -0.04546# -0.0246
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

print("thin",tpl_f)
#tpl_f = Foc.Calc_ThickWaistPos_DeltaOff_UnNormalized(Foc.Calc_K(tpl_n*1e17, gamma) ,tpl_l*1e-4, betastar*100, 0)/100
print("thick",tpl_f)

waist_loc = zvac - tpl_f - tpl_l*1e-6
#For the full 1.01 contour line
#offset_arr = np.linspace(-0.01, 0.01, 201)
#length_arr = np.linspace(tpl_l - 70., tpl_l + 90, 201)

offset_arr = np.linspace(-0.01+.0003, 0.01+.0003, 201)
length_arr = np.linspace(tpl_l - 245., tpl_l + 315, 201)
bmag_image = np.zeros((len(offset_arr),len(length_arr)))

for i in range(len(offset_arr)):
    if i%10 == 0: print(i/len(offset_arr)*100,"%");
    for j in range(len(length_arr)):
        n = np.copy(n_arr); z = np.copy(z_arr)
        tpl_offset = waist_loc + offset_arr[i] + 1/2*tpl_l*1e-6
        tpl_l_set = length_arr[j]
        
        #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z0)
        argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n, tpl_offset*1e6, tpl_n, tpl_l_set, debug)
        
        bmag_image[i][j] = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #print("Bmag CS: ",bmagc)
"""
bmag_image = np.load("/home/chris/Desktop/file.npy")
tpl_l = 736.8569
offset_arr = np.linspace(-0.01+.0003, 0.01+.0003, 201)
length_arr = np.linspace(tpl_l - 245., tpl_l + 315, 201)

#np.save("/home/chris/Desktop/file.npy", bmag_image)
"""
minloc = PProp.PlotContour(bmag_image, offset_arr*100, length_arr-tpl_l, r'$\Delta z_{TPL}$ (cm)', r'$\Delta L_{TPL}\mathrm{\ (\mu m)}$', simple=True, swapx = 0.03, swapy = 23.8)