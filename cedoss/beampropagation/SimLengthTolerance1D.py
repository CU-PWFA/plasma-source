#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  18 12:02:41 2017

Loops over length of simulation to find convergence
for the required length of a VSim simulation

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp
import matplotlib.pyplot as plt

debug = 1

L_up_arr = np.linspace(0.7, 0.5, num = 20)
bmag_image = np.zeros(len(L_up_arr))

#ramp_shape = 'xu4'
nparts = 0

for i in range(len(L_up_arr)):
    L_up = L_up_arr[i]
    
    #Make beam and bulk plasma just as in single_pass
    params = PProp.ReturnDefaultParams(L_up_change = L_up)
    params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
    params['L_ft'] = 0; params['L_dn'] = 0
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    
    Bmag = PProp.CalcBmag(ebeam0, plasma0)
    if debug == 1: print('Bmag: ',Bmag,'L_up: ',L_up)
    bmag_image[i] = Bmag

plt.plot(L_up_arr, bmag_image)
plt.plot(L_up_arr, bmag_image-bmag_image+bmag_image[0])
plt.title("B-mag vs Total Upramp Length")
plt.xlabel(r'$L_{\rm up}$ [m]')
plt.ylabel("B-mag")
plt.grid(); plt.show()

minloc = np.argmin(bmag_image)