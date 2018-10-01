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
from modules import CalcEmitGrowth as W2
import numpy as np
import matplotlib.pyplot as plt



#FACET II
n0 = 1e18 #cm^-3
emit = 3e-6 *100#cm-rad
beta_i = 50 #cm
gam = Foc.gam_def
delta = 0.0025
len_arr = np.linspace(300,1400,101)/1e4 #cm
ymode = 1 #1 for nm, 0 for um
sigmaE = np.sqrt(1/3) * delta
lumi=0

"""
#ILC
n0 = 1e20 #cm^-3
emity = 35e-9 *100 #cm-rad
emitx = 10e-6 *100 #cm-rad
beta_i = 8 #cm
gam = Foc.gam_def*50
sigmaE = 0.00124
delta = np.sqrt(3)*sigmaE
len_arr = np.linspace(50,500,101)/1e4 #cm
ymode = 1 #1 for nm, 0 for um
lumi=1

emit = emitx
"""
d_set = 0


nlen = len(len_arr)

betaf_arr = np.zeros(len(len_arr))
sigo_arr = np.zeros(len(len_arr))
sigi_arr = np.zeros(len(len_arr))
Fval_arr = np.zeros(len(len_arr))
K_arr = np.zeros(len(len_arr))
f_arr = np.zeros(len(len_arr))
betam_arr = np.zeros(len(len_arr))
sigm_arr = np.zeros(len(len_arr))
sig_product_arr = np.zeros(len(len_arr))

for i in range(len(len_arr)):
    L = len_arr[i]
    K = Foc.Calc_K(n0, gam)
    focal = Foc.Calc_Focus_KLength(K, L)
    
    KLls_set = [K, L, Oide.Get_ls_thick(K,L,beta_i,d_set)]
    beta_f = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(K,L,beta_i,d_set)
    
    F_val = Oide.F_Oide(KLls_set)
    sig_min = Oide.Calc_SigMin(F_val, emit)
    sig = Oide.Calc_SigOide(F_val, emit, gam, beta_f)
    sig_ideal = Oide.Calc_SigIdeal(F_val, emit, gam, beta_f)
    beta_min = Oide.Calc_BetaMin(F_val, emit, gam)
    
    if lumi==1:
        sigx = Oide.Calc_SigOide(F_val, emitx, gam, beta_f)
        sigy = Oide.Calc_SigOide(F_val, emity, gam, beta_f)
        sig_product_arr[i]=sigx*sigy
    
    betaf_arr[i] = beta_f
    betam_arr[i] = beta_min
    
    sigo_arr[i] = sig
    sigi_arr[i] = sig_ideal
    sigm_arr[i] = sig_min
    
    Fval_arr[i] = F_val
    f_arr[i] = focal

projbeta_arr = np.zeros(nlen)
bmag_w2_arr = np.zeros(nlen)

k_arr = 1/(len_arr*f_arr)

for x in range(len(projbeta_arr)):
    projbeta_arr[x] = W2.ProjBeta_Thick(k_arr[x], len_arr[x], beta_i, d_set, delta)
    w2 = W2.ThickW2_UnNormalized(k_arr[x], len_arr[x], beta_i, d_set)
    bmag_w2_arr[x] = W2.CalcEmit(w2, sigmaE)

sigc_arr = np.sqrt(projbeta_arr*emit/gam)   

plt.semilogy(len_arr*1e4, Fval_arr)
plt.title("Dimensionless F function vs TPL thickness")
plt.xlabel(r'$L_{pl} \mathrm{\,[\mu m]}$')
plt.ylabel(r'$F(\sqrt{K}L,\sqrt{K}l^*)$')
plt.grid(); plt.show()
    
plt.plot(len_arr*1e4, betaf_arr, label=r'$\beta_f^*$')
plt.plot(len_arr*1e4, betam_arr, label=r'$\beta^*_{opt}$')
plt.title("Beta function at waist vs TPL thickness")
plt.xlabel(r'$L_{pl} \mathrm{\,[\mu m]}$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()
"""
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

plt.rc('font', **font)
"""
lwid = 2.0

if ymode == 0:
    scale = 1e4
    label = r'$\sigma^* \mathrm{\,[\mu m]}$'

if ymode == 1:
    scale = 1e7
    label = r'$\sigma^* \mathrm{\,[n m]}$'

plt.figure(figsize=(3.4*1.5, 2*1.5))
plt.plot(len_arr*np.sqrt(k_arr), sigo_arr*scale,'b-', label="Sync. Rad.", linewidth = lwid)
plt.plot(len_arr*np.sqrt(k_arr), sigm_arr*scale, 'c--', label="Oide Limit", linewidth = lwid)
plt.plot(len_arr*np.sqrt(k_arr), sigc_arr*scale, 'r-', label="Chromaticity", linewidth = lwid)
plt.plot(len_arr*np.sqrt(k_arr), sigi_arr*scale, 'g-', label="Ideal "+r'$\sigma_E=0$', linewidth = lwid)
#plt.title("Beam sizes vs TPL thickness")
plt.xlabel(r'$\sqrt{K}l$')
plt.ylabel(label)
plt.grid(); plt.legend(); plt.show()

plt.plot(len_arr*np.sqrt(k_arr), sigo_arr*scale, label="Oide", linewidth = lwid)
plt.plot(len_arr*np.sqrt(k_arr), sigi_arr*scale, label="Ideal", linewidth = lwid)
plt.plot(len_arr*np.sqrt(k_arr), sigm_arr*scale, label="Minimum", linewidth = lwid)
plt.plot(len_arr*np.sqrt(k_arr), sigc_arr*scale, label="Chromaticity", linewidth = lwid)
plt.title("Beam sizes vs TPL thickness")
plt.xlabel(r'$\sqrt{K}l$')
plt.ylabel(label)
plt.grid(); plt.legend(); plt.show()

if lumi ==1:
    plt.plot(len_arr, sig_product_arr)
    plt.show()
    print(len_arr[np.argmin(sig_product_arr)]*1e4)
    print(len_arr[np.argmin(sig_product_arr)]*np.sqrt(k_arr[0]))