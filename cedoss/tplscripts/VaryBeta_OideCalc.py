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

gam = Foc.gam_def
emit = 500 * 7e-6 *100 #cm-rad
beta_arr = np.logspace(-0.5,1.5, 50) #cm
n0 = 1e18
L = 100e-6 * 100 #cm

sigo_arr = np.zeros(len(beta_arr))
sigi_arr = np.zeros(len(beta_arr))
Fval_arr = np.zeros(len(beta_arr))
betam_arr = np.zeros(len(beta_arr))
sigm_arr = np.zeros(len(beta_arr))
focal_arr = np.zeros(len(beta_arr))
len_arr = np.zeros(len(beta_arr))
betaf_arr = np.zeros(len(beta_arr))

for i in range(len(beta_arr)):
    beta_i = beta_arr[i]
    
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
    
    len_arr[i] = L
    focal_arr[i] = focal
    
plt.semilogx(beta_arr, betam_arr, label=r'$\beta^*_{opt}$')
plt.semilogx(beta_arr, betaf_arr, label=r'$\beta^*_{f}$')
plt.title("Beta function at waist vs inital beta function")
plt.xlabel(r'$\beta^*\mathrm{\,[cm]}$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()

plt.semilogx(beta_arr, sigo_arr*1e4, label="Oide")
plt.semilogx(beta_arr, sigi_arr*1e4, label="Ideal")
plt.semilogx(beta_arr, sigm_arr*1e4, label="Minimum")
plt.title("Beam sizes vs initial beta function")
plt.xlabel(r'$\beta^*\mathrm{\,[cm]}$')
plt.ylabel(r'$\sigma^* \mathrm{\,[\mu m]}$')
plt.grid(); plt.legend(); plt.show()