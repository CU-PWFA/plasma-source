#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  19 13:02:41 2017

Loops over flattop density for tolerance measurements

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp
import matplotlib.pyplot as plt

debug = 1
prl_case=1

gamma_set = 19569.5
dens = 5e16

n_center = 0.138

sig_arr = np.linspace(0.125, 0.155, num = 101)
bmag_image = np.zeros((len(sig_arr)))
nparts = 0

for i in range(len(sig_arr)):
    hwup_set = sig_arr[i]
    
    L_up_set = PProp.def_L_up
    beta_set = PProp.def_beta
    waist_set = PProp.def_waist
    if prl_case == 1:#Using my center calcs
        beta_set = 0.10
        waist_set = -0.3816
        L_up_set = 5*hwup_set
    
    #Make beam and bulk plasma just as in single_pass
    params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set,
                                       beta_change = beta_set, waist_change = waist_set,
                                       L_up_change = L_up_set, hwup_change = hwup_set)
    params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
    params['L_ft'] = 0; params['L_dn'] = 0
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    
    Bmag = PProp.CalcBmag(ebeam0, plasma0)
    if debug == 1: print('hwup: ',hwup_set,'Bmag: ',Bmag)
    bmag_image[i] = Bmag

for j in range(len(bmag_image)):
    if bmag_image[j] < 1.0:
        bmag_image[j] = np.inf
    if bmag_image[j] > 1.175:
        bmag_image[j] = np.inf
    if not np.isfinite(bmag_image[j]):
        bmag_image[j] = 1.0

minloc = np.argmin(bmag_image)
n0min = sig_arr[minloc]
tols = [1.10,1.05,1.01]
for tol in tols:
    print("Tolerance for "+str(tol)+":")
    result = PProp.Calc1DToleranceRange(bmag_image, sig_arr[1]-sig_arr[0], minloc, tol)
    print(" in +/- n:",result[0]/2)
    print(" rel tol :",result[0]/2/n0min)
    
dn = sig_arr[1]-sig_arr[0]
left = sig_arr[0]+dn*result[1]
right = sig_arr[0]+dn*result[2]

plt.plot(sig_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([sig_arr[0],sig_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title("B-mag vs Ramp Half-Length")
plt.xlabel(r'$ \sigma_{hw} $ [m]')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()