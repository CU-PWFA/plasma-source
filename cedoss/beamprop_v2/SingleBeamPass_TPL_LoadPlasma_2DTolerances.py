#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 12:09:25 2018

Using gac cell data from Ken and optimal TPL placement for matching

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit
import pandas as pd
import scipy.interpolate as interp

path = '/home/chris/Desktop/BeamProp/GasCellTest'
excel_filename = '/home/chris/Desktop/Density for tall gas cell lower pinhole .8 density 5.25 us.xlsx'

debug = 0

xl = pd.ExcelFile(excel_filename)
sheet = xl.parse('Density Data')
z_arr = np.array(sheet['Length [m]'])
n_arr = np.array(sheet['Density [cm-3]'])/1e17

extend = 0.0005; num_iter = 40
n_arr_add = np.copy(n_arr); z_arr_add = np.copy(z_arr)
for x in range(num_iter):
    n_arr_add = np.concatenate(([.5*n_arr_add[0]],n_arr_add,[.5*n_arr_add[-1]]))
    z_arr_add = np.concatenate(([z_arr_add[0]-extend],z_arr_add,[z_arr_add[-1]+extend]))

spl = interp.InterpolatedUnivariateSpline(z_arr_add, n_arr_add)
z_arr_ext = np.linspace(z_arr_add[0], z_arr_add[-1], 400000)#8000)
n_arr_ext = spl(z_arr_ext)
z_offset = -z_arr_ext[0]
z_arr_ext = z_arr_ext + z_offset

dump = 10
cores = 4

ramp_end = z_arr[118] + z_offset
endindex = np.nonzero(z_arr_ext>ramp_end)[0][0]

tpl_n = 10.
#tpl_l = 74.91
tpl_f = 0.01475

betastar = .10 #0.00213065326633
waist_loc = -0.000206030150754 - tpl_f

#tpl_offset = waist_loc

off_arr = np.linspace(-0.016,-0.013, 101)
len_arr = np.linspace(65, 85, 51)
bmag_image = np.zeros((len(off_arr),len(len_arr)))
for i in range(len(off_arr)):
    if i%10==0:
        print(i/len(off_arr)*100,"%")
    for j in range(len(len_arr)):
        z_arr = np.copy(z_arr_ext)
        n_arr = np.copy(n_arr_ext)
        
        tpl_offset = off_arr[i]
        tpl_l = len_arr[j]
        
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset)
        
        argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = max(n_arr))
        argon_params['Z'] = z_arr[-1]*1e6; argon_params['Nz']= len(z_arr)
        argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
        argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6, tpl_n, tpl_l, debug)
        
        bmag = PProp.Calc_Bmag(beam_params,n_arr[:endindex], z_arr[:endindex])
        bmag_image[i][j]=bmag
        
minloc = PProp.PlotContour(bmag_image, off_arr*100, len_arr, r'$z_{TPL}$ [cm]', r'$L_{TPL}\,\mathrm{[\mu m]}$')