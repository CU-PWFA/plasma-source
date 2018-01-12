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
waist_arr = np.linspace(-0.5, 0.5, num = 100)
offset_arr = np.linspace(-0.60, -0.20, num = 40)
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

#Optimized for a TPL fixing Bmag 2.0 w/ hw_up 0.18 and position -0.447
tpl_n = 5e16
tpl_l = 0.0008526

"""
#Optimized for a TPL fixing Bmag 2.0 w/ hw_up 0.094 and position -0.310
tpl_n = 5e16
tpl_l = 0.0002526
"""

#Overrides!!
tpl_l = 0.000100

for j in range(len(waist_arr)):
    waist = waist_arr[j]
    #Make beam and bulk plasma just as in single_pass
    #See PlasmaLensContour for param setups
    params = PProp.ReturnDefaultParams(waist_change = waist, hwup_change=0.094)
    params['npart'] = 100 #Want 10,000 to 100,000 for a final figure
    
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
        if debug == 1: print('Bmag: ',Bmag,'z: ',tpl_offset,'wb: ',waist)
        bmag_image[i][j] = Bmag

minloc = PProp.PlotContour(bmag_image, offset_arr, waist_arr, r'$z_{TPL}$ [m]', r'$z_{\rm \beta}$ [m]')
"""
vbeam0 = ebeam0.copy()
tpl_offset = minloc[0]; tpl_l = minloc[1]

plasma0 = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma0)
ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0 = PProp.PropagateVirtual(vbeam0, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, plasma0)
"""