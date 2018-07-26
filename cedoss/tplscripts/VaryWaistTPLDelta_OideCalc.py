#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 10:42:44 2018

Oide effect vs distance delta between the lens and the vacuum waist position

@author: chris
"""

import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import OideCalc as Oide
from modules import CalcEmitGrowth as W2
import numpy as np
import matplotlib.pyplot as plt

gam = Foc.gam_def*20
emit = 10e-6 * 100#3e-6 *100 #cm-rad
beta_i = 5. #cm
d_arr = np.linspace(-10,5,201)
n0 = 1e20
L = 100e-6 * 100 #cm

delta = 0.01

simple = 1

nlen = len(d_arr)

sigo_arr = np.zeros(nlen)
sigi_arr = np.zeros(nlen)
Fval_arr = np.zeros(nlen)
betam_arr = np.zeros(nlen)
sigm_arr = np.zeros(nlen)
focal_arr = np.zeros(nlen)
len_arr = np.zeros(nlen)
betaf_arr = np.zeros(nlen)
ls_arr = np.zeros(nlen)

for i in range(len(d_arr)):
    d = d_arr[i]
    
    K = Foc.Calc_K(n0, gam)
    focal = Foc.Calc_Focus_KLength(K, L)
    KLls_set = [K, L, Oide.Get_ls_corrected(L, focal, beta_i, d)]
    ls_arr[i] = KLls_set[2]
    #print(KLls_set[0]*KLls_set[1]*KLls_set[2])
    F_val = Oide.F_Oide(KLls_set)
    betam_arr[i] = Oide.Calc_BetaMin(F_val, emit, gam)
    betaf_arr[i] = Foc.Calc_BetaStar_DeltaOff(beta_i, focal, d)
    
    sig_min = Oide.Calc_SigMin(F_val, emit)
    sig = Oide.Calc_SigOide(F_val, emit, gam, betaf_arr[i])
    sig_ideal = Oide.Calc_SigIdeal(F_val, emit, gam, betaf_arr[i])
    
    sigo_arr[i] = sig
    sigi_arr[i] = sig_ideal
    sigm_arr[i] = sig_min
    
    Fval_arr[i] = F_val
    
    len_arr[i] = L
    focal_arr[i] = focal
    
plt.plot(d_arr, betam_arr, label=r'$\beta^*_{opt}$')
plt.plot(d_arr, betaf_arr, label=r'$\beta^*_{f}$')
plt.title("Beta function at waist vs delta")
plt.xlabel(r'$\delta\mathrm{\,[cm]}$')
plt.ylabel(r'$\beta \mathrm{\,[cm]}$')
plt.grid(); plt.legend(); plt.show()

kl = 1/focal
projbeta_arr = np.zeros(nlen)
for x in range(len(projbeta_arr)):
    projbeta_arr[x] = W2.ProjBeta_UnNormalized(kl, beta_i, d_arr[x], delta)
sigmaE = 0.60 * delta
bw = beta_i * kl
dw_arr = d_arr * kl
w2_arr = W2.ThinW2_Norm(bw, dw_arr)
bmag_w2_arr = W2.CalcEmit(w2_arr,sigmaE)

sigc_arr = np.sqrt(projbeta_arr*emit*bmag_w2_arr/gam)

plt.plot(d_arr, sigo_arr*1e7, label="Oide")
plt.plot(d_arr, sigi_arr*1e7, label="Ideal")
plt.plot(d_arr, sigm_arr*1e7, label="Minimum")
plt.plot(d_arr, sigc_arr*1e7, label="Chromaticity")
plt.title("Beam sizes vs lens-waist separation")
if simple != 1:
    plt.xlabel(r'$\delta\mathrm{\,[cm]}$')
else:
    plt.xlabel(r'$\mathrm{Lens-waist \ separation \ [cm]}$')
plt.ylabel(r'$\sigma^* \mathrm{\,[n m]}$')
plt.grid(); plt.legend(); plt.show()

print("minimum sigma oide: ",min(sigo_arr*1e7)," nm")
print("minimum sigma chr.: ",min(sigc_arr*1e7)," nm")

"""
plt.plot(delta_arr, ls_arr)
plt.title("Distance from lens to waist vs delta")
plt.xlabel(r'$\delta\mathrm{\,[cm]}$')
plt.ylabel(r'$l^*\mathrm{\,[cm]}$')
plt.grid(); plt.show()

zmin_pred = focal - beta_i
zmin_calc = delta_arr[np.argmin(ls_arr)]
zmin_sigo = delta_arr[np.argmin(sigo_arr)]
print("Predicted min at: ",zmin_pred)
print("Calculated min  : ",zmin_calc)
print("Sigma's min at  : ",zmin_sigo)

betainit_arr = np.zeros(nlen)
betainit_arr = beta_i + np.square(delta_arr)/beta_i

fig, ax1 = plt.subplots()
plt.title("Beam beta and sigma at lens vs delta")
ax1.plot(delta_arr,betainit_arr,color='b')
ax1.set_xlabel(r'$\delta\mathrm{\,[cm]}$')
ax1.set_ylabel(r'$\beta \mathrm{\,[cm]}$', color = 'b')
ax1.tick_params('y', colors = 'b')
ax2 = ax1.twinx()
ax2.plot(delta_arr,Oide.Calc_SigIdeal(0, emit, gam, betainit_arr)*1e4,color='r')
ax2.set_ylabel(r'$\sigma^* \mathrm{\,[\mu m]}$', color = 'r')
ax2.tick_params('y', colors = 'r')
ax1.grid(); plt.show()
"""