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
beta_set = 0.10

kb_arr = np.linspace(50, 2050, num = 401)
bmag_image = np.zeros((len(kb_arr)))

nparts = 0

for i in range(len(kb_arr)):
    kb_set = kb_arr[i]
    
    vals = Approx.CalcApprox(gamma_set, beta_set, kb_set)
    
    dens = vals[2]
    hwup_set = vals[3]
    waist_set = vals[4]
    L_up_set = vals[5]
    
    params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set,
                                           beta_change = beta_set, waist_change = waist_set,
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
    if debug == 1: print('kb: ',kb_set,'Bmag: ',Bmag)
    bmag_image[i] = Bmag
    
for j in range(len(bmag_image)):
    if bmag_image[j] < 1.0:
        bmag_image[j] = np.inf
    if bmag_image[j] > 1e3:
        bmag_image[j] = np.inf
    if not np.isfinite(bmag_image[j]):
        bmag_image[j] = 1.0

plt.plot(kb_arr, bmag_image,label="Empirical Approximation")
plt.plot([kb_arr[0],kb_arr[-1]],[1.01,1.01],ls='--',label='B-mag = 1.01')
plt.title("B-mag vs Empirical Approximation; "+r'$\beta^*$' +" = "+str(100*beta_set)+' cm')
plt.xlabel(r'$ k_b\,\mathrm{[m^{-1}]}$')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()