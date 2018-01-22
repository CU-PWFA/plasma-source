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
"""
debug = 1

npl0_arr = np.linspace(1e16, 2e17, num = 200)
bmag_image = np.zeros((len(npl0_arr)))

nparts = 2

for i in range(len(npl0_arr)):
    dens = npl0_arr[i]
    
    #Make beam and bulk plasma just as in single_pass
    params = PProp.ReturnDefaultParams(npl0_change = dens)
    params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    
    Bmag = PProp.CalcBmag(ebeam0, plasma0)
    if debug == 1: print('npl0: ',dens,'Bmag: ',Bmag)
    bmag_image[i] = Bmag
"""
plt.plot(npl0_arr, bmag_image)
plt.title("B-mag vs Flattop Density")
plt.xlabel("Flattop Density [cm^-3]")
plt.ylabel("B-mag")
plt.grid(); plt.show()

minloc = np.argmin(bmag_image)

levels = [1.01, 1.05, 1.10]
for k in levels:
    print("Tolerance for "+str(k)+":")
    s = PProp.Calc1DTolerance(bmag_image, npl0_arr[1]-npl0_arr[0], minloc, k)
    print(" in +/- x:",s)