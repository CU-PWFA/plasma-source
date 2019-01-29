#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 12:35:49 2018

Calculates Oide values for a range of inputs

@author: chris
"""
import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import OideCalc as Oide
import numpy as np
import matplotlib.pyplot as plt
from modules import CalcEmitGrowth as W2

gam = Foc.gam_def
emit = 7e-6 *100 #cm-rad
beta_i = 10. #cm

L_set = 100e-6 * 100 # set to 0 to calculate from L*k_beta = 0.1
den_arr = np.logspace(18.4, 19.1, 50)
delta=0.01

#L_set = 0 # set to 0 to calculate from L*k_beta = 0.1
#den_arr = np.logspace(18, 20, 50)

sigo_arr = np.zeros(len(den_arr))
sigi_arr = np.zeros(len(den_arr))
Fval_arr = np.zeros(len(den_arr))
betam_arr = np.zeros(len(den_arr))
sigm_arr = np.zeros(len(den_arr))
focal_arr = np.zeros(len(den_arr))
len_arr = np.zeros(len(den_arr))
betaf_arr = np.zeros(len(den_arr))

for i in range(len(den_arr)):
    n0 = den_arr[i]
    if L_set == 0:
        L = 750. * np.sqrt(gam/n0) *100#cm, So that L*k_b = 0.1
    else:
        L = L_set
    
    K = Foc.Calc_K(n0, gam)
    focal = Foc.Calc_Focus_KLength(K, L)
    KLls_set = [K, L, Oide.Get_ls_corrected(L,focal,beta_i)]
    #print(KLls_set[0]*KLls_set[1]*KLls_set[2])
    F_val = Oide.F_Oide(KLls_set)
    betam_arr[i] = Oide.Calc_BetaMin(F_val, emit, gam)
    betaf_arr[i] = Foc.Calc_BetaStar(beta_i, focal)
    
    sig_min = Oide.Calc_SigMin(F_val, emit)
    sig = Oide.Calc_SigOide(F_val, emit, gam, betaf_arr[i])
    sig_ideal = Oide.Calc_SigIdeal(F_val, emit, gam, betaf_arr[i])
    
    sigo_arr[i] = sig
    sigi_arr[i] = sig_ideal
    sigm_arr[i] = sig_min
    
    Fval_arr[i] = F_val
    
    len_arr[i] = L
    focal_arr[i] = focal
    
projbeta_arr = np.zeros(len(den_arr))
bmag_w2_arr = np.zeros(len(den_arr))
sigmaE = np.sqrt(1/3) * delta

k_arr = 1/(len_arr*focal_arr)

for x in range(len(projbeta_arr)):
    projbeta_arr[x] = W2.ProjBeta_Thick(k_arr[x], len_arr[x], beta_i, 0, delta)
    w2 = W2.ThickW2_UnNormalized(k_arr[x], len_arr[x], beta_i, 0)
    bmag_w2_arr[x] = W2.CalcEmit(w2, sigmaE)


sigc_arr = np.sqrt(projbeta_arr*emit*bmag_w2_arr/gam)
    
plt.loglog(den_arr, Fval_arr)
plt.title("Dimensionless F function vs TPL density")
plt.xlabel(r'$n_0 \mathrm{\,[cm^{-3}]}$')
plt.ylabel(r'$F(\sqrt{K}L,\sqrt{K}l^*)$')
plt.grid(); plt.show()
    
#plt.plot(gam_arr, betaf_arr, label=r'$\beta_f^*$')
plt.semilogx(den_arr, betam_arr, label=r'$\beta^*_{opt}$')
plt.semilogx(den_arr, betaf_arr, label=r'$\beta^*_{f}$')
plt.title("Beta function at waist vs TPL density")
plt.xlabel(r'$n_0 \mathrm{\,[cm^{-3}]}$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()

plt.semilogx(den_arr, sigo_arr*1e4, label="Oide")
plt.semilogx(den_arr, sigi_arr*1e4, label="Ideal")
plt.semilogx(den_arr, sigm_arr*1e4, label="Minimum")
plt.semilogx(den_arr, sigc_arr*1e4, label="Chromaticity")
plt.title("Beam sizes vs TPL density")
plt.xlabel(r'$n_0 \mathrm{\,[cm^{-3}]}$')
plt.ylabel(r'$\sigma^* \mathrm{\,[\mu m]}$')
plt.grid(); plt.legend(); plt.show()