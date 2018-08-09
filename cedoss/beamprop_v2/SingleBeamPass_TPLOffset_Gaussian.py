#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 13:40:22 2018

Used with restricted TPL Matching to look at general matching with separation

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

sighw = 0.08 * 1e6
z0 = sighw/1e6*5
"""#More traditional way to frame the setup
vacwaist_flattop = -0.1806 #Distance from vacuum waist to start of flattop
matwaist_flattop = -0.2006 #Distance from matched waist to flattop
z_mat = matwaist_flattop - vacwaist_flattop #Distance from vacuum waist to goal matched waist
"""
#Easier way to match variables in RestrictedTPLMatching.py
z_mat = 0.10 #Distance from vacuum waist to goal matched waist
matwaist_flattop = -0.2006 #Distance from matched waist to flattop
vacwaist_flattop = matwaist_flattop - z_mat #Distance from vacuum waist to start of flattop

betastar = .061608 #Vacuum waist value

tpl_n = 0.5
tpl_l = 351.13250510448523
separ = 0.05

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
waist_loc = vacwaist_flattop
tpl_offset = waist_loc + separ

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                               beta_offset=waist_loc, plasma_start=z0)
beam = PProp.GaussianBeam(beam_params, debug)
argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n, tpl_offset*1e6, tpl_n, tpl_l, debug)

z_fine = np.copy(z)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)

PProp.PlotEmittance(beam,z_fine,m)
PProp.PlotGamma(beam,z_fine,m)

print("Bmag BP: ",PProp.GetBmag(beam,m))
bmagc = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
print("Bmag CS: ",bmagc)

PProp.Plot_CSEvo(beam_params, n, z, z0, legend_loc = 10)

beta_matched = 0.061608
beam_params_matched = PProp.ReturnDefaultElectronParams(path, beta_star=beta_matched,
                                               beta_offset=matwaist_flattop, plasma_start=z0)
PProp.Plot_CSEvo_MatchedCompare(beam_params, beam_params_matched, n, z, z0, legend_loc = 10)