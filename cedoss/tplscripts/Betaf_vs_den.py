#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:46:12 2020

Given initial TPL and beam parameters, spit out all the important info

This one is a template for looping and generating plots

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2
from modules import OideCalc as Oide

#tpl_n = 1e18    # cm^-3
x_arr = np.logspace(16,18.5,num=100)

den_set = np.array([1e16,5e16,8e16,1e17,2e17,3.7e17,5e17,8e17,1e18,2e18])
loc_set = np.array([0, 0.663, 1.083, 1.323, 1.323, 2.1, 2.403, 1.803, 1.683, 1.023])
fit_set = np.array([2.89, 2.82, 2.77, 2.73, 2.69, 2.12, 1.81, 1.48, 1.23, 0.83])
emn_set = np.array([3.147, 3.147, 3.147, 3.148, 3.256, 3.2, 3.224, 4.324, 4.130, 7.101])

calc_beta_set = np.square(fit_set*1e-6)*19569.5/(emn_set*1e-6)*100

tpl_l = 100    # um
tpl_offset = 0.0  # m

gam = Foc.gam_def
emit = 3e-6     # m-rad
beta_i = 0.050   # m

#nbeam  = 6e9
#OR
charge = 1.5e-9 #C
nbeam = charge / 1.6022e-19
sigz = 5.2e-6*100 #cm

#delta_E = 0.0025
#sigma_E = np.sqrt(1/3) * delta_E
#OR
sigma_E = 0.001

y1_arr = np.zeros(len(x_arr))
y2_arr = np.zeros(len(x_arr))

for i in range(len(x_arr)):
    tpl_n = x_arr[i]
    
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
    
    y1_arr[i] = beta_star
    y2_arr[i] = z_waist

plt.loglog(x_arr,y1_arr*100, color='b',label="Theory")
plt.scatter(den_set, calc_beta_set,marker="4",c='k',s=100,label="Fitted from Simulation")
plt.title("Beam betafunction at focus VS plasma lens density")
plt.ylabel(r'$\beta_f^*$'+ " (cm)")
plt.xlabel(r'$n_p\mathrm{\ (cm^{-3})}$')
plt.grid(b=True, which='major', color='0.65', linestyle='-')
plt.grid(b=True, which='minor', color='0.65', linestyle='dotted')
plt.legend(); plt.show()

plt.semilogx(x_arr,y2_arr*100, color='r',label="Theory")
plt.scatter(den_set, loc_set,marker="4",c='k',s=100,label="Simulation")
plt.title("Focus waist location VS plasma lens density")
plt.ylabel(r'$z_w$'+ " (cm)")
plt.xlabel(r'$n_p\mathrm{\ (cm^{-3})}$')
plt.grid(b=True, which='major', color='0.65', linestyle='-')
plt.grid(b=True, which='minor', color='0.65', linestyle='dotted')
plt.legend(); plt.show()