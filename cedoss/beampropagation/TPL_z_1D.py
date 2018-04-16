#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:02:41 2017

Single pass through plasma with optional TPL using Mike's code

@author: chris
"""

import PlasmaPropagation as PProp
import sys
sys.path.insert(0, "../../litos")
import particle_beam as pb
import numpy as np
import matplotlib.pyplot as plt

offs_arr=np.linspace(-.9,-.7,100)
bmag_arr=np.zeros(len(offs_arr))
for i in range(len(offs_arr)):
    params = PProp.ReturnDefaultParams()
    gamma_set = 19569.5
    beta_set = 0.3
    dens = 5e+16
    hwup_set = 0.140314285925
    waist_set = offs_arr[i]# -0.602
    L_up_set = 2.0
    
    params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set,
                                           beta_change = beta_set, waist_change = waist_set,
                                           L_up_change = L_up_set, hwup_change = hwup_set)
    params['L_ft'] = 0; params['L_dn'] = 0
    params['npart'] = 0

    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam_w = PProp.CallMakeBeam(twiss, parts, params)
    
    ebeam0 = PProp.PropagateBackwards(ebeam_w, params)
    
    twisspars = pb.get_twiss(ebeam0,0)
    i_twiss = twisspars[len(twisspars)-1].copy()
    beta = i_twiss["beta"]
    alpha = i_twiss["alpha"]
    
    plasma0 = PProp.MakeBulkPlasma(params)
    
    tpl_n      = 5e16
    tpl_offset = -0.602
    tpl_l = 104.2e-6
    
    plasma0 = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma0)
    
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    bmag_arr[i] = PProp.CalcBmag(ebeam0, plasma0)
    print('waist: ',waist_set,'Bmag: ',bmag_arr[i])

minloc = np.argmin(bmag_arr)
zmin = offs_arr[minloc]
print("min b-mag: ",str(min(bmag_arr))," tpl_l: ",zmin)

plt.plot(offs_arr, bmag_arr)