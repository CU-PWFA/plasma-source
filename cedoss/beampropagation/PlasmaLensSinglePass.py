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

#Make beam and bulk plasma just as in single_pass
"""params = PProp.ReturnDefaultParams(beta_change = 0.15)#, hwup_change=0.10)"""
params = PProp.ReturnDefaultParams()
#params = PProp.ReturnDefaultParams(npl0_change = 5.79e18, gbC_change = 2.34e6, 
#                                   waist_change = -0.4513, hwup_change = 0.16)

gamma_set = 19569.5
beta_set = 0.1
dens = 5e+16
hwup_set = 0.140314285925
waist_set = -0.390533448461
L_up_set = 0.701571429625

params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set,
                                       beta_change = beta_set, waist_change = waist_set,
                                       L_up_change = L_up_set, hwup_change = hwup_set)

params['npart'] = 10#1000
print(params['kb'])
twiss = PProp.CallMakeTwiss(params)
parts = PProp.CallMakeParts(twiss, params)
ebeam_w = PProp.CallMakeBeam(twiss, parts, params)

ebeam0 = PProp.PropagateBackwards(ebeam_w, params)

twisspars = pb.get_twiss(ebeam0,0)
i_twiss = twisspars[len(twisspars)-1].copy()
beta = i_twiss["beta"]
alpha = i_twiss["alpha"]

plasma0 = PProp.MakeBulkPlasma(params)
"""
ebeam_init = ebeam0.copy(); plasma_init = plasma0.copy()
ebeam_init = PProp.PropagatePlasma(ebeam_init, plasma_init)
Bmag_init = PProp.CalcBmag(ebeam_init, plasma_init)
print('initial Bmag: ',Bmag_init)
"""
vbeam0 = ebeam0.copy()
"""
tpl_n      = 5e16
tpl_offset = -0.447#-0.48775510204081629
tpl_l = 852.6e-6#3.9795918367346933e-05

plasma0 = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma0)
"""
ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0 = PProp.PropagateVirtual(vbeam0, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, plasma0)
