#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:42:36 2019

The panultimate final final code for the pwfa figure in Paper 1.  This complete beast includes
two lenses-one for the entrance and one for the exit-as well as TanhLorentz fits for the lenses
to make the scenario as ultra-realistic as possible.

...although I did get rid of the 

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

debug = 1
path = '/media/chris/New Volume/BeamProp/testGaussian'
gamma = PProp.def_gamma

tpl_n = 0.3
tpl_l = 736.9
sighw = 0.02542 * 1e6
zvac = -0.0455
betastar = 0.05
zvac2 = 0.10 + 0.0455 + 0.0246 + 0.00003 #+ 0.5*tpl_l2*1e-6
    
a = 376.409454963
b = 68.0148653818
n_0 = 0.294112354406

#These might not be accurate, but in any case they should be large enough
#to not matter in the first place.  From ref. ana., tanh works fine
#A = 3.780e3*1000000
#gammal = 4812.9*1000000
#x_0 = 0

z0 = sighw/1e6*5
tpl_l2 = tpl_l

fit_tanh = [a,b,n_0]
#fit_lorentz = [A, gammal, x_0]

argon_params = PProp.ReturnDefaultPlasmaParams(path, sigma_hw = sighw, plasma_start = z0, scaledown = 10)
argon = PProp.GaussianRampPlasma(argon_params, debug)

n = argon.nez
z = np.linspace(0,argon_params['Z'],len(n))/1e6

dump = 10
cores = 4

ramp_end = argon_params['z0']/1e6
endindex = np.nonzero(z>ramp_end)[0][0]

focal = Foc.Calc_Focus_Square_SI(tpl_n*1e17, tpl_l/1e6, gamma)
#focal = 1 #set this to bypass case 11 being weird with l=0
beta_f = Foc.Calc_BetaStar(betastar, focal)
tpl_f = focal*(1-beta_f/betastar)
waist_loc = zvac - tpl_f
tpl_offset = waist_loc

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                               beta_offset=waist_loc, plasma_start=z0)
beam = PProp.GaussianBeam(beam_params, debug)
argon = PProp.CustomPlasma_ThinPlasmaLens_Tanh(argon_params, n, tpl_offset*1e6, fit_tanh, debug)
argon = PProp.CustomPlasma_ThinPlasmaLens_Tanh(argon_params, n, zvac2*1e6, fit_tanh, debug)
""" If you want to turn off the actaul beam prop
z_fine = np.copy(z)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)

PProp.PlotEmittance(beam,z_fine,m)
PProp.PlotGamma(beam,z_fine,m)

print("Bmag BP: ",PProp.GetBmag(beam,m))
#"""
bmagc = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
print("Bmag CS: ",bmagc)

betaf, alphaf, gammaf, gbf = PProp.Plot_CSEvo(beam_params, n, z, z0, legend_loc = 10)

print(1/gammaf[-1])