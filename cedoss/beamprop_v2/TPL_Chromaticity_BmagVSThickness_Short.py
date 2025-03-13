#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:51:00 2019

This version only propagates through the lens to measure emittance growth.  None
of that propagating-to-the-focus-to-measure-beam-size since we already know it
is close to theory.

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
debug = 0#1
zmult=2
"""
num = 101
len_arr = np.linspace(500, 7000, num) #original is 200 to 1200 for 1e18 cm-3
#len_arr = np.linspace(80, 1200, num)
emit_arr = np.zeros(num)
betamin_arr = np.zeros(num)
tpl_f_arr = np.zeros(num)
centloc_arr = np.zeros(num)
betapromin_arr = np.zeros(num)
centbet_arr = np.zeros(num)

sigm_arr = np.zeros(num)
betam_arr = np.zeros(num)

position_error = 0

gammab = PProp.def_gamma
tpl_n = .3
betastar = .05 #0.00213065326633

delta = 0.05#0.01

for k in range(len(len_arr)):
    if k%2 == 0: print(k/len(len_arr)*100,"%");
    
    tpl_l = len_arr[k]
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    tpl_f_arr[k] = tpl_f
    
    leftext = 0.001 #1
    rightext = 0.04 #3
    
    z_arr = np.linspace(-1e-5, tpl_l*1e-6+1e-5, int((tpl_l+20)+1)*zmult) + (position_error / 1e6)
    n_arr = np.zeros(len(z_arr))+1e-9
    
    waist_loc = 0.
    tpl_offset = waist_loc
    z_offset = -z_arr[0]
    z_arr = z_arr + z_offset
    
    e_spec = np.array([0, -delta, delta]) + 1.0
    
    beam_params0 = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset,
                                                       gamma=gammab)
    argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = tpl_n)
    argon_params['Z'] = z_arr[-1]*1e6
    argon_params['Nz']= len(z_arr)
    argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error + tpl_l/2, tpl_n, tpl_l, debug)
    
    ###############################################################################
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                    plasma_start=z_offset, gamma=gammab)
    
    gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr, betapro_arr  = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)
    
    #emit_arr[k] = bmag_arr[-1]
    betapromin_arr[k] = min(betapro_arr)

    beam = PProp.GaussianBeam(beam_params, debug)
    dump = 10; cores = 4
    z_fine = np.copy(z_arr)*1e6
    PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

    m = int(len(z_fine)/dump)
    emit_arr[k] = PProp.GetBmag(beam,m)


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)
lwid = 2.0

kl = 1/tpl_f_arr
bw_arr = betastar * kl
dw = 0
sigmaE = 0.577 * delta
bmag_w2_arr = W2.CalcEmit_OLD(W2.ThinW2_Norm(bw_arr, dw),sigmaE)

#Thick
k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100
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
plt.figure(figsize=(3.4*1.5, 2*1.5))
#plt.title("Emittance Growth vs Lens Thickness")
plt.plot(len_arr*1e-6*np.sqrt(k), emit_arr, 'r-', label = "Simulated",linewidth = lwid)
plt.plot(len_arr*1e-6*np.sqrt(k), bmag_w2_arr, 'g-.', label = "Thin Calculated",linewidth = lwid)
plt.plot(len_arr*1e-6*np.sqrt(k), bmag_w2_arr_thick, 'b--', label = "Thick Calculated",linewidth = lwid)
plt.ylabel(r'$\epsilon_f/\epsilon_0$')
plt.xlabel(r'$\mathrm{Lens \ Thickness \ }\sqrt{K}L $')
#plt.ylim(1.0,1.0015)#1.0,1.01 for 1e18
plt.xlim(0.1, 1.1)#0.1, 1.5 for 1e18
plt.grid(); plt.legend(); plt.show()