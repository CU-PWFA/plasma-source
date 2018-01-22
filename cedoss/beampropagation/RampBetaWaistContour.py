#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  18 12:02:41 2017

Loops over ramp length and vacuum betafunction waist position
for a given ramp profile

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp

debug = 1

#waist_arr = np.linspace(-0.2, 0.2, num = 4)
#width_arr = np.linspace(0.00000, 0.0001, num = 4)
waist_arr = np.linspace(-0.5, -0.1, num = 100)
width_arr = np.linspace(0.05, 0.20, num = 100)
bmag_image = np.zeros((len(waist_arr),len(width_arr)))

#ramp_shape = 'xu4'
ramp_shape = 'gauss'
nparts = 100

for i in range(len(waist_arr)):
    for j in range(len(width_arr)):
        waist = waist_arr[i]
        hwup_set = width_arr[j]
        
        #Make beam and bulk plasma just as in single_pass
        params = PProp.ReturnDefaultParams(waist_change = waist, hwup_change = hwup_set, ramp_change=ramp_shape)
        params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
        
        twiss = PProp.CallMakeTwiss(params)
        parts = PProp.CallMakeParts(twiss, params)
        ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
        ebeam0 = PProp.PropagateBackwards(ebeam0, params)
        
        plasma0 = PProp.MakeBulkPlasma(params)
        ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
        
        Bmag = PProp.CalcBmag(ebeam0, plasma0)
        if debug == 1: print('Bmag: ',Bmag,'sig: ',hwup_set,'wb: ',waist)
        bmag_image[i][j] = Bmag

minloc = PProp.PlotContour(bmag_image, waist_arr, width_arr, r'$z_{\rm \beta *}$ [m]', r'$\sigma_{hw}$ [m]')

levels = [1.01, 1.05, 1.10]
for k in levels:
    print("Tolerance for "+str(k)+":")
    PProp.Calc2DTolerance(bmag_image, waist_arr, width_arr, k)

waist = minloc[0]; hwup_set = minloc[1]

#Very condensed single pass
params = PProp.ReturnDefaultParams(waist_change = waist, hwup_change = hwup_set, ramp_change=ramp_shape)
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