#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 12:35:49 2018

Calculates values from Oide paper

@author: chris
"""
import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import OideCalc as Oide
import numpy as np
import matplotlib.pyplot as plt
"""
n0 = 1e18 #cm^-3
emit = 7e-6 * 100# 7 cm-rad
beta_f = .1 #cm
gam = Foc.gam_def*5
L = 110e-6 *100 #cm
"""

n0 = 4.9e12 #cm^-3
emit = 100e-6 * 100# cm-rad
#beta_f = 43e-6 *100 #cm
beta_f = .0007 *100 #cm
gam = 30
L = 20000e-6 *100 #cm

K = Foc.Calc_K(n0, gam)
focal = Foc.Calc_Focus_KLength(K, L)
#focal = Foc.Calc_Focus_Square_SI(n0, L, gam)/100 #gives focal length in m
KLls_set = [K, L, Oide.Get_ls(L,focal)]
KLls_set = [K, L, Oide.Get_ls_corrected(L,focal,10)]
"""
f = 0.6 * 100
L = 2/3*f
K = 2.05/(L*f)
emit = 2.5e-8 * 100
gam = Foc.gam_def * 50
KLls_set = [K, L, f-1/2*L]
beta_f = 43e-6 * 100
"""
F_val = Oide.F_Oide(KLls_set)
print("F val : ",F_val)

sig_min = Oide.Calc_SigMin(F_val, emit)*1e4
print("sigmin: ",sig_min," [um]")

sig = Oide.Calc_SigOide(F_val, emit, gam, beta_f)*1e4
print("sig   : ",sig," [um]")

sig_ideal = Oide.Calc_SigIdeal(F_val, emit, gam, beta_f)*1e4
print("sig_ideal : ",sig_ideal," [um]")

beta_min = Oide.Calc_BetaMin(F_val, emit, gam)
print("beta min: ",beta_min," [cm]")

#Section to plot Oide's Fig 3, probably dont need ever again
"""
x_arr = np.linspace(.02,5,100) #L/l array
y = [1,1.5,2] #values of K*L*l
F1 = np.zeros(len(x_arr)); F2 = np.zeros(len(x_arr)); F3 = np.zeros(len(x_arr))
lf = 1
for i in range(len(x_arr)):
    L = lf * x_arr[i]
    K1 = y[0]/lf**2/x_arr[i]
    F1[i] = Oide.F_Oide([K1, L, lf - L/2])
    K2 = y[1]/lf**2/x_arr[i]
    F2[i] = Oide.F_Oide([K2, L, lf - L/2])
    K3 = y[2]/lf**2/x_arr[i]
    F3[i] = Oide.F_Oide([K3, L, lf - L/2])
    
plt.semilogy(x_arr, F3, label = "KLf = 2")
plt.semilogy(x_arr, F2, label = "KLf = 1.5")
plt.semilogy(x_arr, F1, label = "KLf = 1")
plt.xlim(0,5)
plt.ylim(.01,100)
plt.ylabel("F")
plt.xlabel("L/f")
plt.grid(); plt.legend(); plt.show()
"""