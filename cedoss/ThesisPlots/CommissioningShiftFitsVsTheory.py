#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 23:11:36 2023

Goal is to plot the fit parameters and theory vs density

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2


betafit_eslice_arr = np.array([[2.31,  1.74,  2.1,   1.57,  1.43,  2.58,  000],
                               [1.23,  1.11,  1.04,  1.46,  1.92,  3.01,  1.57],
                               [0.992, 0.335, 0.591, 0.494, 0.527, 0.976, 1.65],
                               [1.18,  0.626, 0.564, 0.357, 0.521, 1.31,  0.902],
                            [3.34E-01,1.43E+00,0.687,0.436, 0.372, 0.779, 0.894],
                             [0.409, 4.67E-01, 0.595,0.729, 0.589, 0.711, 1.03]])

emitn_eslice_arr = np.array([[1.90E-05, 2.43E-05, 1.88E-05, 2.71E-05, 3.95E-05, 3.56E-05, 0000],
                             [5.08E-05, 4.46E-05, 4.23E-05, 2.99E-05, 2.78E-05, 3.19E-05, 1.31E-04],
                             [8.03E-05, 3.22E-04, 1.11E-04, 8.30E-05, 9.47E-05, 1.02E-04, 7.93E-05],
                             [1.37E-04, 1.47E-04, 9.40E-05, 7.60E-05, 7.03E-05, 5.29E-05, 1.88E-05],
                             [4.72E-05, 1.19E-04, 1.50E-04, 1.06E-04, 9.51E-05, 1.16E-04, 1.67E-04],
                             [3.06E-04, 2.49E-04, 1.54E-04, 1.06E-04, 1.28E-04, 1.41E-04, 1.56E-04]])

betafit_fullproj_arr = np.array([151.5, 104.5, 84.19, 80.12, 55.56, 63.11])*1e-2

emitn_fullproj_arr = np.array([42.81, 71.02, 128.63, 106.79, 190.83, 190.46])

density_arr = np.array([0.001, 0.27, 1.62, 6.48, 15.3, 31.5])
nbeam_est = 1.7

e_slice_arr = np.array([9.85, 9.90, 9.95, 10.00, 10.05, 10.10, 10.15])
colors = plt.cm.autumn(np.linspace(1, 0, len(e_slice_arr)))

tpl_l_arr = np.array([0.00001, 10.07, 10.07, 13.01, 15.49, 17.55])*1e3
beta_i = 1.51
tpl_offset = 0
gam = 19569.5
beta_theory_nonlin = np.zeros(len(density_arr))
for j in range(len(density_arr)):
    tpl_k = Foc.Calc_K(density_arr[j]*1e16, gam)*100*100
    beta_theory_nonlin[j] = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(tpl_k, tpl_l_arr[j]*1e-6, beta_i, tpl_offset)

density_arr_smooth = np.logspace(-3.,1.4,10000)
beta_theory_smooth = np.zeros(len(density_arr_smooth))
tpl_l_set = 10.07*1e3
for j in range(len(density_arr_smooth)):
    tpl_k = Foc.Calc_K(density_arr_smooth[j]*1e16, gam)*100*100
    beta_theory_smooth[j] = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(tpl_k, tpl_l_set*1e-6, beta_i, tpl_offset)


plt.figure(figsize=(8,4.5))
for i in range(len(e_slice_arr)):
    if i == len(e_slice_arr)-1:
        plt.semilogx(density_arr[1:],betafit_eslice_arr[1:,i],label="%.2f"%e_slice_arr[i]+" GeV",lw=1,c=colors[i])
        plt.scatter(density_arr[1:],betafit_eslice_arr[1:,i],s=10,c=colors[i])
    else:
        plt.semilogx(density_arr,betafit_eslice_arr[:,i],label="%.2f"%e_slice_arr[i]+" GeV",lw=1,c=colors[i])
        plt.scatter(density_arr[1:],betafit_eslice_arr[1:,i],s=10,c=colors[i])

plt.semilogx(density_arr,betafit_fullproj_arr,lw=2,c="black",label="Full Projection")
#plt.semilogx(density_arr,beta_theory_nonlin,lw=2,ls='dashed',c="orange")
plt.semilogx(density_arr_smooth,beta_theory_smooth,lw=2,ls='dotted',c="blue",label="Theory")

plt.plot([nbeam_est,nbeam_est],[-1,5],ls='dotted',c="black")

plt.legend()
plt.xlim([0.001,32])
plt.ylim([0,4.1])
plt.xlabel("Plasma Density Estimate "+r'$(\mathrm{\times 10^{16} \ cm^{-3}})$')
plt.ylabel("Waist Betafunction "+r'$\beta^*$'+r'$(\mathrm{m})$')
plt.show()

###NOw for emittance

emit_i = 42.81 #um-rad
sigma_E = 0.02
emitn_theory_smooth = np.zeros(len(density_arr_smooth))
for j in range(len(density_arr_smooth)):
    tpl_k = Foc.Calc_K(density_arr_smooth[j]*1e16, gam)*100*100
    w2_thick = W2.ThickW2_UnNormalized(tpl_k, tpl_l_set*1e-6, beta_i, tpl_offset)
    em_growth = W2.CalcEmit(w2_thick, sigma_E)
    emitn_theory_smooth[j] = em_growth*emit_i
"""
scattering_arr = np.zeros(len(density_arr_smooth))
eleq = 1.602e-19
light = 3e8
ep0 = 8.85e-12
me = 9.11e-31
ra = 1e-10
zion=2
qs=1
re = 2.818e-15
for j in range(len(density_arr_smooth)):
    den = density_arr_smooth[j]*1e16*100**3
    kp = den*eleq**2/(me*ep0*light**2)
    rb = 5.31e5/np.sqrt(den/100**3)*1e-2
    #kp = 1/rb
    scatter = qs * (np.log(rb/ra)+(1.78*zion*(zion+1))/qs*np.log(297/np.sqrt(zion)))
    delta_emit = kp * re * scatter * tpl_l_set*1e-6 /(np.sqrt(2*gam))
    scattering_arr[j] = emit_i+delta_emit*1e6
"""
emitn_eslice_arr = emitn_eslice_arr*1e6
plt.figure(figsize=(8,4.5))
for i in range(len(e_slice_arr)):
    if i == len(e_slice_arr)-1:
        plt.semilogx(density_arr[1:],emitn_eslice_arr[1:,i],label="%.2f"%e_slice_arr[i]+" GeV",lw=1,c=colors[i])
        plt.scatter(density_arr[1:],emitn_eslice_arr[1:,i],s=10,c=colors[i])
    else:
        plt.semilogx(density_arr,emitn_eslice_arr[:,i],label="%.2f"%e_slice_arr[i]+" GeV",lw=1,c=colors[i])
        plt.scatter(density_arr[1:],emitn_eslice_arr[1:,i],s=10,c=colors[i])

plt.plot([nbeam_est,nbeam_est],[1e-6,500],ls='dotted',c="black")

plt.semilogx(density_arr,emitn_fullproj_arr,lw=2,c="black",label="Full Projection")
#plt.semilogx(density_arr,beta_theory_nonlin,lw=2,ls='dashed',c="orange")
plt.semilogx(density_arr_smooth,emitn_theory_smooth,lw=2,ls='dotted',c="blue",label="Theory")
#plt.semilogx(density_arr_smooth,scattering_arr,lw=2,ls='dotted',c="red",label="Theory")


plt.xlabel("Plasma Density Estimate "+r'$(\mathrm{\times 10^{16} \ cm^{-3}})$')
plt.ylabel("Normalized Emittance "+r'$\epsilon_N \ $'+r'$(\mathrm{\mu m-rad})$')

plt.ylim([0,330])
plt.xlim([-0.1,32])

plt.legend()
plt.show()






















