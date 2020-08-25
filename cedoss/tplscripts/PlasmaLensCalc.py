#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 11:56:28 2018

Given initial TPL and beam parameters, spit out all the important info

Useful also as a future template for calculating things

@author: chris
"""

import numpy as np
import sys

sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2
from modules import OideCalc as Oide

tpl_n = 1e18    # cm^-3
tpl_l = 50    # um
tpl_offset = 0.0  # m

gam = Foc.gam_def
emit = 3e-6     # m-rad
beta_i = 5.0   # m

#nbeam  = 6e9
#OR
charge = 1.5e-9 #C
nbeam = charge / 1.6022e-19
sigz = 5.2e-6*100 #cm

#delta_E = 0.0025
#sigma_E = np.sqrt(1/3) * delta_E
#OR
sigma_E = 0.001

rho_max = nbeam/(2*np.pi)**(3/2)/sigz/(beta_i*100*emit*100/gam)

tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n, tpl_l, gam)/100
tpl_k = Foc.Calc_K(tpl_n, gam)*100*100

w2_thick = W2.ThickW2_UnNormalized(tpl_k, tpl_l*1e-6, beta_i, tpl_offset)
em_growth = W2.CalcEmit(w2_thick, sigma_E)
projbeta_thick = W2.ProjBeta_Thick_Gauss(tpl_k, tpl_l*1e-6, beta_i, tpl_offset, sigma_E)

beta_star = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(tpl_k, tpl_l*1e-6, beta_i, tpl_offset)
z_waist = Foc.Calc_ThickWaistPos_DeltaOff_UnNormalized(tpl_k, tpl_l*1e-6, beta_i, tpl_offset)

#For some reason needs to be in cgs (probably constants in SigOide)
KLls_set = [tpl_k/100/100, tpl_l*1e-4, Oide.Get_ls_thick(tpl_k/100/100, tpl_l*1e-4, beta_i*100, tpl_offset*100)]
F_val = Oide.F_Oide(KLls_set)
sig_oide = Oide.Calc_SigOide(F_val, emit*100, gam, beta_star*100)/100

print("Focal length [m]: ", tpl_f)
print("sqrt(K)*L: ", np.sqrt(tpl_k)*tpl_l*1e-6)
print("Beam Density: ",rho_max)
print()
print("The following from Thick lens equations")
print("Focusing strength K [m^-2]: ", tpl_k)
print("Emittance growth: ", em_growth)
print("Final Emittance [m-rad]: ", em_growth * emit)
print("Centroid Focus location [m]: ", z_waist)
print("Centorid Focus beta [m]: ", beta_star)
print("Focus rms size [m]: ", np.sqrt(projbeta_thick*emit/gam))
print()
print("Oide Limit Calculations")
print("SR beam size [m]: ", sig_oide)