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
#params = PProp.ReturnDefaultParams()
params = PProp.ReturnDefaultParams(npl0_change = 1e17, gbC_change = 45000)
#params = PProp.ReturnDefaultParams(hwup_change = 0.00010, ramp_change='xu4', waist_change = 0.3)#, hwup_change=0.10)

params['npart'] = 1000

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
