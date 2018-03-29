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

dens = 5e16

gamma_arr = np.linspace(1e4, 3.5e4, num=301)
bmag_image = np.zeros((len(gamma_arr)))
nparts = 0

for i in range(len(gamma_arr)):
    gamma_set = gamma_arr[i]
    L_up_set = PProp.def_L_up
    hwup_set = PProp.def_hwup
    beta_set = PProp.def_beta
    waist_set = PProp.def_waist
    
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
    if debug == 1: print('gamma: ',gamma_set,'Bmag: ',Bmag)
    bmag_image[i] = Bmag

for j in range(len(bmag_image)):
    if bmag_image[j] < 1.0:
        bmag_image[j] = np.inf
    if bmag_image[j] > 1.175:
        bmag_image[j] = np.inf
    if not np.isfinite(bmag_image[j]):
        bmag_image[j] = 1.0

minloc = np.argmin(bmag_image)
n0min = gamma_arr[minloc]
tols = [1.10,1.05,1.01]
for tol in tols:
    print("Tolerance for "+str(tol)+":")
    result = PProp.Calc1DToleranceRange(bmag_image, gamma_arr[1]-gamma_arr[0], minloc, tol)
    print(" in +/- n:",result[0]/2)
    print(" rel tol :",result[0]/2/n0min)
    
dn = gamma_arr[1]-gamma_arr[0]
left = gamma_arr[0]+dn*result[1]
right = gamma_arr[0]+dn*result[2]

plt.plot(gamma_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([gamma_arr[0],gamma_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.scatter(n0min,min(bmag_image),color='k',label='min n0 = '+str(n0min))
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title('B-mag vs Initial Beam '+r'$\gamma$')
plt.xlabel(r'$\gamma$')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()