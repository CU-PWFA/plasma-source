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
import EmpiricalApprox as Approx

debug = 1

gamma_set = 19569.5
kb_option = 2

beta_arr = np.linspace(0.025, 0.40, num = 301)
bmag_image = np.zeros((len(beta_arr)))

nparts = 0

for i in range(len(beta_arr)):
    beta = beta_arr[i]
    
    vals = Approx.CalcApprox(gamma_set, beta, kb_option)
    
    dens = vals[2]
    hwup_set = vals[3]
    waist_set = vals[4]
    L_up_set = vals[5]
    
    params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set,
                                           beta_change = beta, waist_change = waist_set,
                                           L_up_change = L_up_set, hwup_change = hwup_set)
    params['npart'] = nparts
    params['L_ft'] = 0; params['L_dn'] = 0
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    
    Bmag = PProp.CalcBmag(ebeam0, plasma0)
    if debug == 1: print('beta: ',beta,'Bmag: ',Bmag)
    bmag_image[i] = Bmag
    
for j in range(len(beta_arr)):
    if bmag_image[j] < 1.0:
        bmag_image[j] = np.inf
    if bmag_image[j] > 1e3:
        bmag_image[j] = np.inf
    if not np.isfinite(bmag_image[j]):
        bmag_image[j] = 1.0

plt.plot(100*beta_arr, bmag_image,label="Empirical Approximation")
plt.plot([100*beta_arr[0],100*beta_arr[-1]],[1.01,1.01],ls='--',label='B-mag = 1.01')
plt.title("B-mag vs Empirical Approximation; kb = "+str(Approx.Getkb(kb_option))+r'$m^{-1}$')
plt.xlabel(r'$ \beta^* $ [cm]')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()