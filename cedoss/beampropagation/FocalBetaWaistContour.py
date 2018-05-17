#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  13 12:02:41 2017

Loops over TPL location and vacuum betafunction waist position

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp

debug = 1


#waist_arr = np.linspace(-0.14, 0.14, num = 20) - 0.36
waist_arr = np.linspace(-0.8, -0.4, num = 101)#100)
offset_arr = np.linspace(-0.65, -0.45, num = 101)#40)
bmag_image = np.zeros((len(offset_arr),len(waist_arr)))
"""
#Optimized for a TPL fixing Bmag 1.10 and position -0.48
tpl_n = 5e16
tpl_l = 3.9795918367346933e-05
"""
"""
#Optimized for a TPL fixing Bmag 1.68 and position -0.33
tpl_n = 5e16
tpl_l = 0.000175555555556
"""
"""
#Optimized for a TPL fixing Bmag 2.0 w/ hw_up 0.18 and position -0.447
tpl_n = 5e16
tpl_l = 0.0008526
"""
"""
#Optimized for a TPL fixing Bmag 2.0 w/ hw_up 0.094 and position -0.310
tpl_n = 5e16
tpl_l = 0.0002526
"""

tpl_n = 5e16
#tpl_offset = -0.602
tpl_l = 104.2e-6

for j in range(len(waist_arr)):
    #Make beam and bulk plasma just as in single_pass
    #See PlasmaLensContour for param setups
    gamma_set = 19569.5
    beta_set = 0.3
    dens = 5e+16
    hwup_set = 0.140314285925
    waist_set = waist_arr[j]# -0.602
    L_up_set = 2.0
    
    params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set,
                                           beta_change = beta_set, waist_change = waist_set,
                                           L_up_change = L_up_set, hwup_change = hwup_set)
    params['npart'] = 0
    params['L_ft'] = 0; params['L_dn'] = 0
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)

    for i in range(len(offset_arr)):
        tpl_offset = offset_arr[i]
        plasma = plasma0.copy(); ebeam = ebeam0.copy(); #vbeam = vbeam0.copy()
        plasma = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma)
        
        ebeam = PProp.PropagatePlasma(ebeam, plasma)
        
        Bmag = PProp.CalcBmag(ebeam, plasma)
        if debug == 1: print('Bmag: ',Bmag,'z: ',tpl_offset,'wb: ',waist_set)
        bmag_image[i][j] = Bmag

minloc = PProp.PlotContour(bmag_image, offset_arr, waist_arr, r'$z_{TPL}$ [m]', r'$z_{\rm \beta}$ [m]')
"""
levels = [1.01, 1.05, 1.10]
for k in levels:
    print("Tolerance for "+str(k)+":")
    PProp.Calc2DTolerance(bmag_image, offset_arr, waist_arr, k)
"""
tpl_offset = minloc[0]; waist = minloc[1]
print(tpl_offset, waist)
"""
#Very condensed single pass
params = PProp.ReturnDefaultParams(waist_change = waist, hwup_change = hwup_set)
params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
twiss = PProp.CallMakeTwiss(params)
parts = PProp.CallMakeParts(twiss, params)
ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
ebeam0 = PProp.PropagateBackwards(ebeam0, params)
plasma0 = PProp.MakeBulkPlasma(params)
vbeam0 = ebeam0.copy()
plasma0 = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma0)
ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0 = PProp.PropagateVirtual(vbeam0, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, plasma0)
"""