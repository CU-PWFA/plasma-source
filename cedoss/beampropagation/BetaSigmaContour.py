#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  19 12:02:41 2017

Loops over sigma_hw and beam beta

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp

debug = 1

#sigma_arr = np.logspace(15, 18.7, num = 30)
#beta_arr = np.logspace(0, 2, num = 30)
sigma_arr = np.linspace(.05, .25, num = 30)
beta_arr = np.linspace(.01, 0.2, num = 30)
bmag_image = np.zeros((len(sigma_arr),len(beta_arr)))

nparts = 0
dens = 4.6e16

for i in range(len(sigma_arr)):
    for j in range(len(beta_arr)):
        hwup_set = sigma_arr[i]
        beta_set = beta_arr[j]
        
        params = PProp.ReturnDefaultParams(beta_change=beta_set, hwup_change=hwup_set,
                                           npl0_change=dens,)
        params['npart'] = nparts
        params['L_ft'] = 0; params['L_dn'] = 0
        
        twiss = PProp.CallMakeTwiss(params)
        parts = PProp.CallMakeParts(twiss, params)
        ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
        ebeam0 = PProp.PropagateBackwards(ebeam0, params)
        
        plasma0 = PProp.MakeBulkPlasma(params)
        ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
        
        Bmag = PProp.CalcBmag(ebeam0, plasma0)
        if debug == 1:
            print('Bmag: ',Bmag,'sig: ',hwup_set,'beta: ',beta_set)
        if np.isnan(Bmag):
            Bmag = np.inf
        bmag_image[i][j] = Bmag

np.save('/home/chris/Desktop/DataLoads/ContourBetaSigma/bmagarr.npy',bmag_image)
np.save('/home/chris/Desktop/DataLoads/ContourBetaSigma/sig.npy',sigma_arr)
np.save('/home/chris/Desktop/DataLoads/ContourBetaSigma/beta.npy',beta_arr)
