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

n0 = 5e16 #cm^-3
emit = 7e-6 #um-rad
beta_i = .10 #m
gam = Foc.gam_def

len_arr = np.linspace(50e-6,1000e-6,50) #m
betaf_arr = np.zeros(len(len_arr))

sigo_arr = np.zeros(len(len_arr))
sigi_arr = np.zeros(len(len_arr))
Fval_arr = np.zeros(len(len_arr))
betam_arr = np.zeros(len(len_arr))
sigm_arr = np.zeros(len(len_arr))

for i in range(len(len_arr)):
    L = len_arr[i]
    
    focal = Foc.Calc_Focus_Square(n0, L*1e6) #gives focal length in cm
    beta_f = Foc.Calc_BetaStar(beta_i*100, focal)/100
    KLls_set = Oide.Get_KLls(L,focal)
    
    F_val = Oide.F_Oide(KLls_set)
    sig_min = Oide.Calc_SigMin(F_val, emit)
    sig = Oide.Calc_SigOide(F_val, emit, gam, beta_f)
    sig_ideal = Oide.Calc_SigIdeal(F_val, emit, gam, beta_f)
    beta_min = Oide.Calc_BetaMin(F_val, emit, gam)
    
    betaf_arr[i] = beta_f
    betam_arr[i] = beta_min
    
    sigo_arr[i] = sig
    sigi_arr[i] = sig_ideal
    sigm_arr[i] = sig_min
    
    Fval_arr[i] = F_val
    
plt.semilogy(len_arr*1e6, Fval_arr)
plt.xlabel(r'$L_{pl} \mathrm{\,[\mu m]}$')
plt.ylabel(r'$F(\sqrt{K}L,\sqrt{K}l^*)$')
plt.grid(); plt.show()
    
plt.plot(len_arr*1e6, betaf_arr*100, label=r'$\beta_f^*$')
plt.plot(len_arr*1e6, betam_arr*100, label=r'$\beta^*_{min}$')
plt.title("Beta function at waist vs TPL thickness")
plt.xlabel(r'$L_{pl} \mathrm{\,[\mu m]}$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()

plt.plot(len_arr*1e6, betaf_arr*100, label=r'$\beta_f^*$')
plt.plot(len_arr*1e6, betam_arr*100, label=r'$\beta^*_{min}$')
plt.title("Beta function at waist vs TPL thickness")
plt.xlabel(r'$L_{pl} \mathrm{\,[\mu m]}$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()