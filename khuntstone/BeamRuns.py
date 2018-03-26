#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 19:09:17 2018

@author: keenan
"""

import KH_PlasmaProp as PProp
import sys
sys.path.insert(0, "../litos")

# Function to get 
def getWaists(exit_beam, plasma):
    return 0
     
params = PProp.ReturnDefaultParams()


params['npart'] = 1000 #Want 10,000 to 100,000 for a final figure

twiss   = PProp.CallMakeTwiss(params)
parts   = PProp.CallMakeParts(twiss, params)
ebeam0  = PProp.CallMakeBeam(twiss, parts, params)
ebeam0  = PProp.PropagateBackwards(ebeam0, params)
plasma0 = PProp.MakeBulkPlasma(params)
ebeam0  = PProp.PropagatePlasma(ebeam0, plasma0)
# Get waist locations at entrance and exit
#s_f0 = plasma0['s']; s_f0 = s_f0[len(s_f0) - 1]
#a_i0 = ebeam0[0]['alpha']
#g_i0 = ebeam0[0]['gamma']
#s_entr_waist0 = a_i0/g_i0;
#a_f0 = ebeam0[len(ebeam0)-1]['alpha']
#g_f0 = ebeam0[len(ebeam0)-1]['gamma']
#s_exit_waist0 = s_f0 + a_f0/g_f0

# Unmatched case
params['npl0'] = 1e18; 
twiss   = PProp.CallMakeTwiss(params)
parts   = PProp.CallMakeParts(twiss, params)
ebeam1  = PProp.CallMakeBeam(twiss, parts, params)
ebeam1  = PProp.PropagateBackwards(ebeam1, params)
plasma1 = PProp.MakeBulkPlasma(params)
ebeam1  = PProp.PropagatePlasma(ebeam1, plasma1)
#s_f1 = plasma1['s']; s_f1 = s_f1[len(s_f1) - 1]
#a_i1 = ebeam1[0]['alpha']
#g_i1 = ebeam1[0]['gamma']
#s_entr_waist1 = a_i1/g_i1;
#a_f1 = ebeam1[len(ebeam1)-1]['alpha']
#g_f1 = ebeam1[len(ebeam1)-1]['gamma']
#s_exit_waist1 = s_f1 + a_f1/g_f1

