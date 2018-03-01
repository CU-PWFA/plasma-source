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

npl0_arr = np.logspace(15, 18.7, num = 300)
gamma_arr = np.logspace(2.5, 6.4, num = 300)
bmag_image = np.zeros((len(npl0_arr),len(gamma_arr)))

nparts = 0
hwup_set= 0.12
waist_set = -0.3226

for i in range(len(npl0_arr)):
    for j in range(len(gamma_arr)):
        dens = npl0_arr[i]
        gamma = gamma_arr[j]
        if np.sqrt(dens)/gamma > 5e5:
            print("spagoot")
        params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma,
                                           hwup_change = hwup_set, waist_change = waist_set)
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
            print('Bmag: ',Bmag,'npl0: ',dens,'gbC: ',gamma)
        if np.isnan(Bmag):
            Bmag = np.inf
        bmag_image[i][j] = Bmag
"""
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testarr.npy',bmag_image)
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testx.npy',npl0_arr)
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testy.npy',gamma_arr)
"""
"""
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testarr2.npy',bmag_image)
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testx2.npy',npl0_arr)
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testy2.npy',gamma_arr)
"""

np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testarr3.npy',bmag_image)
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testx3.npy',npl0_arr)
np.save('/home/chris/Desktop/DataLoads/ContourDensityGamma/testy3.npy',gamma_arr)
