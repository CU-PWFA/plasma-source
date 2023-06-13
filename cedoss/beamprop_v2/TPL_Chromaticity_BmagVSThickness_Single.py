#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 12:02:15 2023

Copy of TPL_Chromaticity_LargeOffset built specifically to measure emittance
growth vs lens thickness

This version is a single propagation, looking to plot sigma vs z for numerical and theory.

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2

path = '/media/chris/New Volume/BeamProp/Placeholder'
debug = 1
zmult=1

num = 101
#len_arr = np.linspace(500, 7000, num) #original is 200 to 1200 for 1e18 cm-3
tpl_l = 7000 #um

gammab = PProp.def_gamma
tpl_n = .3
betastar = .05 #0.00213065326633
emit_n = 3e-6
gam = 19569.5

delta = 0.05#1
    
tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100

leftext = 0.005 #1
rightext = 3 #3

z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult)
n_arr = np.zeros(len(z_arr))+1e-9

waist_loc = 0.
tpl_offset = waist_loc
z_offset = -z_arr[0]
z_arr = z_arr + z_offset

e_spec = np.array([0, -delta, delta]) + 1.0

skipProp = True

if not skipProp:
    beam_params0 = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset,
                                                       gamma=gammab)
    argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = tpl_n)
    argon_params['Z'] = z_arr[-1]*1e6
    argon_params['Nz']= len(z_arr)
    argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + tpl_l/2, tpl_n, tpl_l, debug)
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                    plasma_start=z_offset, gamma=gammab)
    
    gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr, betapro_arr  = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)
    
    betapromin = min(betapro_arr)
    
    beam = PProp.GaussianBeam(beam_params, debug)
    dump = 10; cores = 4
    z_fine = np.copy(z_arr)*1e6
    PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)
    
    m = int(len(z_fine)/dump)
    emit_f = PProp.GetBmag(beam,m)
    sigmar_min = PProp.GetSigmaMin(beam, m)
    beta_min = PProp.GetBetaMin(beam, m)

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)
lwid = 2.0

kl = 1/tpl_f
bw_arr = betastar * kl
dw = 0
sigmaE = 0.577 * delta
bmag_w2_arr = W2.CalcEmit_OLD(W2.ThinW2_Norm(bw_arr, dw),sigmaE)

#PProp.PlotSigmar(beam,z_arr,m)

k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100
sqrtk = np.sqrt(k)
b0 = betastar*sqrtk
g0 = 1/b0
a0 = 0
l0 = tpl_l*1e-6*sqrtk
zw = ((b0-g0)*np.sin(l0)*np.cos(l0)+a0*np.cos(2*l0))/(b0*np.square(np.sin(l0))+g0*np.square(np.cos(l0))+a0*np.sin(2*l0))
zstar = zw/sqrtk

s = Foc.Calc_ThickWaistPos_DeltaOff_UnNormalized(k,tpl_l*1e-6,betastar,0)
beta_ideal = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(k, tpl_l*1e-6,betastar,0)
sigma_ideal = np.sqrt(beta_ideal*emit_n/gam)
sigma_w2 = np.sqrt(beta_ideal*emit_n*emit_f/gam)

print("Beam Size from Numerical",sigmar_min)
print("Beam Size from Thick Ideal",sigma_ideal)
print("Beam Size from Thick Theory",sigma_w2)
print("")
oldideal = np.sqrt(Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f, 0)*3e-6/gammab)
print("Beam Size from 'Old Ideal'",oldideal)
"""
#Thick

bmag_w2_arr_thick = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    w2 = W2.ThickW2_UnNormalized(k, len_arr[x]*1e-6, betastar, 0)
    bmag_w2_arr_thick[x] = W2.CalcEmit(w2, sigmaE)

projbeta_arr = np.zeros(len(len_arr))
projbeta_arr2 = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    projbeta_arr[x] = W2.ProjBeta_UnNormalized(kl[x], betastar, 0, delta)
    b0 = betastar; a0 = 0; g0 = 1/betastar
    projbeta_arr2[x] = W2.ProjBetaCS_UnNormalized(k, len_arr[x]*1e-6, b0, a0, g0, delta)
    
projbetathick_arr = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    projbetathick_arr[x] = W2.ProjBeta_Thick(k, len_arr[x]*1e-6, betastar, 0, delta)
"""
