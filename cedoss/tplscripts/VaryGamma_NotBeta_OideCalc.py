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

n0 = 1e18 #cm^-3
emit = 7e-6 *100 #cm-rad
beta_i = 10. #cm

#len_arr = np.linspace(50e-6,300e-6,50) #m
#betaf_arr = np.zeros(len(len_arr))

gam_arr = np.linspace(1e4, 5e4, 50)

sigo_arr = np.zeros(len(gam_arr))
sigi_arr = np.zeros(len(gam_arr))
Fval_arr = np.zeros(len(gam_arr))
betam_arr = np.zeros(len(gam_arr))
sigm_arr = np.zeros(len(gam_arr))
focal_arr = np.zeros(len(gam_arr))
betaf_arr = np.zeros(len(gam_arr))

for i in range(len(gam_arr)):
    gam = gam_arr[i]
    L = 750. * np.sqrt(gam/n0) *100#cm, So that L*k_b = 0.1
    
    K = Foc.Calc_K(n0, gam)
    focal = Foc.Calc_Focus_KLength(K, L)
    KLls_set = [K, L, Oide.Get_ls(L,focal)]
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
    
    focal_arr[i] = focal
    
plt.semilogy(gam_arr, Fval_arr)
plt.title("Dimensionless F function vs TPL thickness")
plt.xlabel(r'$\gamma_b$')
plt.ylabel(r'$F(\sqrt{K}L,\sqrt{K}l^*)$')
plt.grid(); plt.show()
    
#plt.plot(gam_arr, betaf_arr, label=r'$\beta_f^*$')
plt.plot(gam_arr, betam_arr, label=r'$\beta^*_{opt}$')
plt.plot(gam_arr, betaf_arr, label=r'$\beta^*_{f}$')
plt.title("Beta function at waist vs TPL thickness")
plt.xlabel(r'$\gamma_b$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()

plt.plot(gam_arr, sigo_arr*1e4, label="Oide")
plt.plot(gam_arr, sigi_arr*1e4, label="Ideal")
plt.plot(gam_arr, sigm_arr*1e4, label="Minimum")
plt.title("Beam sizes vs TPL thickness")
plt.xlabel(r'$\gamma_b$')
plt.ylabel(r'$\sigma^* \mathrm{\,[\mu m]}$')
plt.grid(); plt.legend(); plt.show()