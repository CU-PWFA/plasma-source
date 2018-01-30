#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  18 12:02:41 2017

Loops over vacuum betafunction waist position
for a ramp shape to match

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp
import matplotlib.pyplot as plt

debug = 1

waist_arr = np.linspace(-0.6, -0.2, num = 200)
bmag_image = np.zeros(len(waist_arr))

#ramp_shape = 'xu4'
ramp_shape = 'gauss'
nparts = 0
hwup_set = .12

for i in range(len(waist_arr)):
    waist = waist_arr[i]
    
    #Make beam and bulk plasma just as in single_pass
    params = PProp.ReturnDefaultParams(waist_change = waist, hwup_change = hwup_set, ramp_change=ramp_shape)
    params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
    params['L_ft'] = 0; params['L_dn'] = 0
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    
    Bmag = PProp.CalcBmag(ebeam0, plasma0)
    if debug == 1: print('Bmag: ',Bmag,'wb: ',waist)
    bmag_image[i] = Bmag

plt.plot(waist_arr, bmag_image, label=r'$\sigma_{hw}$ = '+ str(hwup_set) +' [m]')
plt.title("B-mag vs Waist Position")
plt.xlabel(r'$z_{\rm \beta *}$ [m]')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()

minloc = np.argmin(bmag_image)

levels = [1.01, 1.05, 1.10]
for k in levels:
    print("Tolerance for "+str(k)+":")
    s = PProp.Calc1DTolerance(bmag_image, waist_arr[1]-waist_arr[0], minloc, k)
    print(" in +/- x:",s)
    
print(waist_arr[minloc])