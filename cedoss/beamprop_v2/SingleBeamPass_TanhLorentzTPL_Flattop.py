#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 15:39:05 2019

TPL matching into hard cuttof flattop plasma source

This version uses a PTPL with DoubleTanh-Lorentzian density
distribution, and will cater towards the PWFA example in my
PTPL Paper 1.

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

debug = 1
path = '/home/chris/Desktop/BeamProp/testGaussian'
gamma = PProp.def_gamma

case = 1
if case == 1:
    sighw = 0.00001
    tpl_n = 0.5
    tpl_l = 995
    zvac = -0.02084 - .5*tpl_l*1e-6
    betastar = .10
if case == 2:
    sighw = 0.00001
    tpl_n = 1.0
    tpl_l = 596
    zvac = -0.01770015 - .5*tpl_l*1e-6
    betastar = .10
    
    a = 298.465129594
    b = 78.1121221798
    n_0 = 1.00251813855

    A = 7.56e3
    gamma = 4812.9
    x_0 = 0

    fit_tanh = [a,b,n_0]
    fit_lorentz = [A, gamma, x_0]

z0 = 0.1


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

waist_loc = zvac# - tpl_f
tpl_offset = waist_loc

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                               beta_offset=waist_loc, plasma_start=z0)
beam = PProp.GaussianBeam(beam_params, debug)
argon = PProp.CustomPlasma_ThinPlasmaLens_TanhLorentz(argon_params, n, tpl_offset*1e6, fit_tanh, fit_lorentz, debug)

z_fine = np.copy(z)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)

PProp.PlotEmittance(beam,z_fine,m)
PProp.PlotGamma(beam,z_fine,m)

print("Bmag BP: ",PProp.GetBmag(beam,m))
bmagc = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
print("Bmag CS: ",bmagc)

set_sta = 8000; set_fin = 12000
PProp.Plot_CSEvo(beam_params, n, z, z0, legend_loc = 0, subset = [set_sta, set_fin])