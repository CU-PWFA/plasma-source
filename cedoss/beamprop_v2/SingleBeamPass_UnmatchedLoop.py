#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 09:09:09 2019

For the cases where a TPL could match, set thickness to zero and loop over
initial betafunction waists to find the initial best matching value.

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
    betastar = .10
if case == 2:
    sighw = 0.01 * 1e6
    tpl_n = 0.5
    tpl_l = 607.5
    zvac = -0.01372
    betastar = .10

if case == 11:
    tpl_n = 0.5
    tpl_l = 0
    sighw = 0.08 * 1e6
    zvac = -0.2006
    betastar = .061608
    
if case == 20: #442.1um to match 5cm beta into a ramp that requires 2.5cm beta
    tpl_n = 0.5
    tpl_l = 442.1
    sighw = 0.0273 * 1e6
    zvac = -0.054
    betastar = 0.05

tpl_l=0
z0 = sighw/1e6*5
zvac_arr = np.linspace(-0.1,0,101)
bmag_arr = np.zeros(len(zvac_arr))

for i in range(len(zvac_arr)):
    zvac = zvac_arr[i]
    argon_params = PProp.ReturnDefaultPlasmaParams(path, sigma_hw = sighw, plasma_start = z0, scaledown = 10)
    argon = PProp.GaussianRampPlasma(argon_params, debug)
    
    n = argon.nez
    z = np.linspace(0,argon_params['Z'],len(n))/1e6
    
    dump = 10
    cores = 4
    
    ramp_end = argon_params['z0']/1e6
    endindex = np.nonzero(z>ramp_end)[0][0]
    
    #focal = Foc.Calc_Focus_Square_SI(tpl_n*1e17, tpl_l/1e6, gamma)
    #focal = 1 #set this to bypass case 11 being weird with l=0
    #beta_f = Foc.Calc_BetaStar(betastar, focal)
    #tpl_f = focal*(1-beta_f/betastar)
    #waist_loc = zvac - tpl_f
    #tpl_offset = waist_loc
    
    #Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=zvac, plasma_start=z0)
    beam = PProp.GaussianBeam(beam_params, debug)
    #argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n, tpl_offset*1e6, tpl_n, tpl_l, debug)
    
    #z_fine = np.copy(z)*1e6
    #PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)
    
    #m = int(len(z_fine)/dump)
    
    #PProp.PlotEmittance(beam,z_fine,m)
    #PProp.PlotGamma(beam,z_fine,m)
    
    #print("Bmag BP: ",PProp.GetBmag(beam,m))
    bmag_arr[i] = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
plt.plot(zvac_arr,bmag_arr)
plt.show()
print(min(bmag_arr))