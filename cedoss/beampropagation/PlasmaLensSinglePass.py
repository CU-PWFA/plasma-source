#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:02:41 2017

@author: chris
"""

import PlasmaPropagation as PProp

#Make beam and bulk plasma just as in single_pass
params = PProp.ReturnDefaultParams(beta_change = 0.15)#, hwup_change=0.10)
params['npart'] = 1000

twiss = PProp.CallMakeTwiss(params)
parts = PProp.CallMakeParts(twiss, params)
ebeam0 = PProp.CallMakeBeam(twiss, parts, params)

ebeam0 = PProp.PropagateBackwards(ebeam0, params)

plasma0 = PProp.MakeBulkPlasma(params)

ebeam_init = ebeam0.copy(); plasma_init = plasma0.copy()
ebeam_init = PProp.PropagatePlasma(ebeam_init, plasma_init)
Bmag_init = PProp.CalcBmag(ebeam_init, plasma_init)
print('initial Bmag: ',Bmag_init)

vbeam0 = ebeam0.copy()

tpl_n      = 5e16
tpl_offset = -0.48775510204081629
tpl_l = 3.9795918367346933e-05

plasma0 = PProp.InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma0)
ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0 = PProp.PropagateVirtual(vbeam0, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, plasma0)
