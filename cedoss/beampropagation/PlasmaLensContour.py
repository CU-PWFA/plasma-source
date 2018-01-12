#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 12:02:41 2017

Loops over TPL thickness and TPL location

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp

debug = 1

#Make beam and bulk plasma just as in single_pass
#For B_M=1.10
#params = PProp.ReturnDefaultParams(beta_change = 0.15)
#For B_M=1.68
#params = PProp.ReturnDefaultParams(beta_change = 0.15, hwup_change=0.10)

#For B_M=2.00 and a longer ramp
#params = PProp.ReturnDefaultParams(hwup_change=0.18)
params = PProp.ReturnDefaultParams(hwup_change=0.18, waist_change = 0.10)
#For B_M=2.00 and a shorter ramp
#params = PProp.ReturnDefaultParams(hwup_change=0.094)

params['npart'] = 100 #Want 10,000 to 100,000 for a final figure

twiss = PProp.CallMakeTwiss(params)
parts = PProp.CallMakeParts(twiss, params)
ebeam0 = PProp.CallMakeBeam(twiss, parts, params)

ebeam0 = PProp.PropagateBackwards(ebeam0, params)

plasma0 = PProp.MakeBulkPlasma(params)

#num should be much higher, but 10 can give a good rough estimation
offset_arr = np.linspace(-0.60, -0.20, num = 10)
length_arr = np.linspace(0.0, 300e-6, num = 10)

#offset_arr = np.linspace(-0.50, -0.20, num = 20)
#length_arr = np.linspace(0e-6, 300e-6, num = 20)

bmag_image = np.zeros((len(offset_arr),len(length_arr)))
tpl_n      = 5e16 #raise in discrete simulation sets

ebeam_init = ebeam0.copy(); plasma_init = plasma0.copy()
ebeam_init = PProp.PropagatePlasma(ebeam_init, plasma_init)
Bmag_init = PProp.CalcBmag(ebeam_init, plasma_init)
print('initial Bmag: ',Bmag_init)

for i in range(len(offset_arr)):
    for j in range(len(length_arr)):
        tpl_offset = offset_arr[i]
        tpl_l      = length_arr[j]
        
        if tpl_l == 0.0:
            bmag_image[i][j] = Bmag_init
        else:
            plasma = plasma0.copy(); ebeam = ebeam0.copy(); #vbeam = vbeam0.copy()
            plasma = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma)
            
            ebeam = PProp.PropagatePlasma(ebeam, plasma)
            
            Bmag = PProp.CalcBmag(ebeam, plasma)
            if debug == 1: print('Bmag: ',Bmag,'z: ',tpl_offset,'l: ',tpl_l)
            bmag_image[i][j] = Bmag

minloc = PProp.PlotContour(bmag_image, offset_arr, length_arr, r'$z_{TPL}$ [m]', r'$TPL_{\rm L}$ [m]')

vbeam0 = ebeam0.copy()
tpl_offset = minloc[0]; tpl_l = minloc[1]

plasma0 = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma0)
ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0 = PProp.PropagateVirtual(vbeam0, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, plasma0)
