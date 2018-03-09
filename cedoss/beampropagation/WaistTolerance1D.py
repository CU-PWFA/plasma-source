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

waist_arr = np.linspace(-0.40, -0.36, num = 201)
bmag_image = np.zeros((len(waist_arr)))
nparts = 0

for i in range(len(waist_arr)):
    waist_set = waist_arr[i]
    
    L_up_set = PProp.def_L_up
    beta_set = PProp.def_beta
    hwup_set = PProp.def_hwup
    if prl_case == 1:#Using my center calcs
        beta_set = 0.10
        hwup_set = 0.1384
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
    if debug == 1: print('waist: ',waist_set,'Bmag: ',Bmag)
    bmag_image[i] = Bmag

for j in range(len(bmag_image)):
    if bmag_image[j] < 1.0:
        bmag_image[j] = np.inf
    if bmag_image[j] > 1.175:
        bmag_image[j] = np.inf
    if not np.isfinite(bmag_image[j]):
        bmag_image[j] = 1.0

minloc = np.argmin(bmag_image)
n0min = waist_arr[minloc]
tols = [1.10,1.05,1.01]
for tol in tols:
    print("Tolerance for "+str(tol)+":")
    result = PProp.Calc1DToleranceRange(bmag_image, waist_arr[1]-waist_arr[0], minloc, tol)
    print(" in +/- n:",result[0]/2)
    print(" rel tol :",result[0]/2/n0min)
    
dn = waist_arr[1]-waist_arr[0]
left = waist_arr[0]+dn*result[1]
right = waist_arr[0]+dn*result[2]

plt.plot(waist_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([waist_arr[0],waist_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title("B-mag vs Waist Location")
plt.xlabel(r'$ z_{\beta^*} $ [m]')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()