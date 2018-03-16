#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  19 13:02:41 2017

Loops over flattop density for tolerance measurements

@author: chris
"""

import numpy as np
import os
import PlasmaPropagation as PProp
import EmpiricalApprox as Approx

debug = 1

gamma_set = 19569.5

kb_arr = np.linspace(50, 350, num = 201)
beta_arr = np.linspace(0.05, 0.65, num = 201)
bmag_image = np.zeros((len(kb_arr),len(beta_arr)))

nparts = 0

for i in range(len(kb_arr)):
    for j in range(len(beta_arr)):
        kb_set = kb_arr[i]
        beta_set = beta_arr[j]
        
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
        if debug == 1: print('kb: ',kb_set, 'beta', beta_set,'Bmag: ',Bmag)
        bmag_image[i][j] = Bmag

path = '/home/chris/Desktop/DataLoads/ContourEmpiricalApprox1/'
if not os.path.exists(path):
    print("Creating new directory")
    os.makedirs(path)
np.save(path + 'bmagarr.npy',bmag_image)
np.save(path + 'kb.npy',kb_arr)
np.save(path + 'beta.npy',beta_arr)