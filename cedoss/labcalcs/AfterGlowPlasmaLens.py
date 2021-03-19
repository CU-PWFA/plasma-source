#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:34:07 2020

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2

gasden_arr = np.array([1e16,5e16,8e16,1e17,2e17,3.7e17,5e17,8e17,1e18,2e18])
#emitn_arr = np.array([3.147,3.147,3.147,3.148,3.256,3.2,3.224,4.324,4.130,7.101])
emitn_arr = np.array([3.0916e-6, 3.0916e-6, 3.0917e-6, 3.0923e-6, 
                      3.2682e-6, 3.1634e-6, 3.1994e-6, 4.2866e-6,
                      4.1653e-6, 6.9898e-6])
sig_arr = np.array([2.93e-6, 2.86e-6, 2.82e-6, 2.74e-6, 2.40e-6, 2.11e-6, 1.76e-6, 1.30e-6, 1.17e-6, 0.82e-6])

denarr_pred = np.linspace(gasden_arr[0],gasden_arr[-1],100)
sigarr_pred = np.zeros(len(denarr_pred))
emtarr_pred = np.zeros(len(denarr_pred))

beta_i = 0.05
tpl_l = 71.9426
tpl_offset = 0
emit = 3.0914e-6
sigma_E = 0.001
gam = Foc.gam_def

for i in range(len(denarr_pred)):
    tpl_n = denarr_pred[i]
    tpl_k = Foc.Calc_K(tpl_n, gam)*100*100
    beta_star = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(tpl_k, tpl_l*1e-6, beta_i, tpl_offset)
    w2_thick = W2.ThickW2_UnNormalized(tpl_k, tpl_l*1e-6, beta_i, tpl_offset)
    projbeta_thick = W2.ProjBeta_Thick_Gauss(tpl_k, tpl_l*1e-6, beta_i, tpl_offset, sigma_E)
    em_growth = W2.CalcEmit(w2_thick, sigma_E)
    
    emtarr_pred[i] = emit*em_growth
    sigarr_pred[i] = np.sqrt(projbeta_thick*emit*em_growth/gam)
    
plt.title("Witness Beam Final Emittance vs Density")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("Normalized Emittance " + r'$(mm-mrad)$')
plt.plot(gasden_arr[2:],emitn_arr[2:]*1e6, c='blue', ls = 'solid')
#plt.plot(denarr_pred, emtarr_pred)
plt.legend(); plt.grid(); plt.show()

plt.title("Witness Beam Focus Spot Size vs Density")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("Spot Size " + r'$(\mu m)$')
plt.plot(gasden_arr,sig_arr*1e6, c='blue', ls = 'solid', label='Simulation')
plt.plot(denarr_pred, sigarr_pred*1e6, c='red', ls = '--', label='Ideal Nonlinear Theory')

plt.legend(); plt.grid(); plt.show()