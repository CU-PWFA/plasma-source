#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  19 12:02:41 2017

Loops over flattop density and centroid relativistic lorentz factor

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp

debug = 1

npl0_arr = np.linspace(1e17, 2e17, num = 20)
gamma_arr = np.linspace(2.5e4, 12.5e4, num = 20)
bmag_image = np.zeros((len(npl0_arr),len(gamma_arr)))

nparts = 10

for i in range(len(npl0_arr)):
    for j in range(len(gamma_arr)):
        dens = npl0_arr[i]
        gamma = gamma_arr[j]
        
        #Make beam and bulk plasma just as in single_pass
        params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma)
        params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
        
        twiss = PProp.CallMakeTwiss(params)
        parts = PProp.CallMakeParts(twiss, params)
        ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
        ebeam0 = PProp.PropagateBackwards(ebeam0, params)
        
        plasma0 = PProp.MakeBulkPlasma(params)
        ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
        
        Bmag = PProp.CalcBmag(ebeam0, plasma0)
        if debug == 1: print('Bmag: ',Bmag,'npl0: ',dens,'gbC: ',gamma)
        bmag_image[i][j] = Bmag

minloc = PProp.PlotContour(bmag_image, npl0_arr, gamma_arr, r'$ npl0 $ [cm^-3]', r'$ \gamma $')

levels = [1.01, 1.05, 1.10]
for k in levels:
    print("Tolerance for "+str(k)+":")
    PProp.Calc2DTolerance(bmag_image, npl0_arr, gamma_arr, k)

dens = minloc[0]; gamma = minloc[1]

#Very condensed single pass
params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma)
params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
twiss = PProp.CallMakeTwiss(params)
parts = PProp.CallMakeParts(twiss, params)
ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
ebeam0 = PProp.PropagateBackwards(ebeam0, params)
plasma0 = PProp.MakeBulkPlasma(params)
vbeam0 = ebeam0.copy()
ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0 = PProp.PropagateVirtual(vbeam0, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, plasma0)