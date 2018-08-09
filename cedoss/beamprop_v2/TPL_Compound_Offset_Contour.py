#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 13:11:36 2018

Compound Plasma Lenses, and with offset between first lens and vacuum waist

This one is for a 2d contour of x and d
@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

path = '/home/chris/Desktop/BeamProp/GasCellTest'
debug = 0

zmult=1

gammab = PProp.def_gamma

tpl_n = 10.

tpl_l1 = 200
tpl_l2 = 400
tpl_f1 = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l1, gammab)/100
tpl_f2 = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l2, gammab)/100
leftext = 1 #1
rightext = 3 #3

d_arr = np.linspace(-0.25, -0.10, 101)
x_arr = np.linspace(0.0002, 0.006, 101)
betamin_image = np.zeros((len(d_arr),len(x_arr)))
for k in range(len(d_arr)):
    if k%10 == 0: print(k/len(d_arr)*100,"%"); 
    for m in  range(len(x_arr)):
        position_error = d_arr[k] * 1e6
        tpl_x = x_arr[m] #m
        
        z_arr = np.linspace(-leftext*tpl_f1, rightext*tpl_f1, int((leftext+rightext)*tpl_f1*1e6+1)*zmult) + (position_error / 1e6)
        n_arr = np.zeros(len(z_arr))
        
        betastar = .10 #0.00213065326633
        waist_loc = 0.
        tpl_offset = waist_loc
        z_offset = -z_arr[0]
        z_arr = z_arr + z_offset
        
        tpl_1loc = tpl_offset
        tpl_2loc = tpl_offset + tpl_x
        
        e_spec = np.array([0, -0.01, 0.01]) + 1.0
        
        beam_params0 = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                           beta_offset=waist_loc, plasma_start=z_offset,
                                                           gamma=gammab)
        argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                           nset = tpl_n)
        argon_params['Z'] = z_arr[-1]*1e6
        argon_params['Nz']= len(z_arr)
        argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
        argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_1loc*1e6 + position_error + tpl_l1/2, tpl_n, tpl_l1, debug)
        argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_2loc*1e6 + position_error + tpl_l2/2, tpl_n, tpl_l2, debug)
        
        maxbetacomp = np.zeros(len(z_arr))
        crange = 120*zmult
        for i in range(len(e_spec)):
        #Make beam and bulk plasma just as in single_pass
            beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                           beta_offset=waist_loc, plasma_start=z_offset,
                                                           gamma=gammab * e_spec[i])
            beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
            for j in range(len(maxbetacomp)):
                if beta[j] > maxbetacomp[j]:
                    maxbetacomp[j] = beta[j]
        betamin_image[k][m] = min(maxbetacomp)

minloc = PProp.PlotContour_General(betamin_image*1e6, d_arr*1e2, x_arr*1e3, r'$d$ [cm]', r'$x$ [mm]', r'$\beta_{min} \ \mathrm{[\mu m]}$')