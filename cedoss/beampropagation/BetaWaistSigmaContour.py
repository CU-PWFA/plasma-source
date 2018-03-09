#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  19 12:02:41 2017

Loops over sigma_hw and beam beta waist

@author: chris
"""

import numpy as np
import os
import PlasmaPropagation as PProp

debug = 1

sigma_arr = np.linspace(0.09, 0.18, num = 100)
zbeta_arr = np.linspace(-0.50, -0.25, num = 100)
bmag_image = np.zeros((len(zbeta_arr),len(sigma_arr)))

betalabel = 10

gamma_set = 19569.5

dens = 5.0e16 # First, normal run
#dens = 1.0e17# Second, high density run
#dens = 1.0e16# Third, low density run

beta_const = betalabel/100

for i in range(len(zbeta_arr)):
    for j in range(len(sigma_arr)):
        hwup_set = sigma_arr[j]
        zbeta_set = zbeta_arr[i]
        
        params = PProp.ReturnDefaultParams(beta_change=beta_const, hwup_change=hwup_set,
                                           npl0_change=dens,     waist_change=zbeta_set,
                                           L_up_change=5*hwup_set, gbC_change=gamma_set)
        params['npart'] = 0
        params['L_ft'] = 0; params['L_dn'] = 0
        
        twiss = PProp.CallMakeTwiss(params)
        parts = PProp.CallMakeParts(twiss, params)
        ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
        ebeam0 = PProp.PropagateBackwards(ebeam0, params)
        
        plasma0 = PProp.MakeBulkPlasma(params)
        ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
        
        Bmag = PProp.CalcBmag(ebeam0, plasma0)
        if debug == 1:
            print('Bmag: ',Bmag,'sig: ',hwup_set,'zbeta: ',zbeta_set)
        if np.isnan(Bmag):
            Bmag = np.inf
        bmag_image[i][j] = Bmag

#path = '/home/chris/Desktop/DataLoads/ContourBetaWaistSigma_First/'+str(betalabel)+'cm/'
path = '/home/chris/Desktop/DataLoads/ContourBetaWaistSigma_10GeV_HighRes_PostFix/'+str(betalabel)+'cm/'

if not os.path.exists(path):
    print("Creating new directory")
    os.makedirs(path)
np.save(path + 'bmagarr.npy',bmag_image)
np.save(path + 'sig.npy',sigma_arr)
np.save(path + 'beta.npy',zbeta_arr)
