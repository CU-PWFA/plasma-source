#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 13:11:56 2018

Using gas cell data from Ken and beam parameters for matching from Keenan

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as interp

path = '/home/chris/Desktop/BeamProp/GasCellTest'
excel_filename = '/home/chris/Desktop/Density for tall gas cell lower pinhole .8 density 5.25 us.xlsx'

debug = 1

xl = pd.ExcelFile(excel_filename)
sheet = xl.parse('Density Data')
z_arr = np.array(sheet['Length [m]'])
n_arr = np.array(sheet['Density [cm-3]'])/1e17

extend = 0.0005; num_iter = 10
n_arr_add = np.copy(n_arr); z_arr_add = np.copy(z_arr)
for x in range(num_iter):
    n_arr_add = np.concatenate(([.5*n_arr_add[0]],n_arr_add,[.5*n_arr_add[-1]]))
    z_arr_add = np.concatenate(([z_arr_add[0]-extend],z_arr_add,[z_arr_add[-1]+extend]))

spl = interp.InterpolatedUnivariateSpline(z_arr_add, n_arr_add)
z_arr_ext = np.linspace(z_arr_add[0], z_arr_add[-1], 8000)
n_arr_ext = spl(z_arr_ext)
z_offset = -z_arr_ext[0]
z_arr_ext = z_arr_ext + z_offset
if debug == 1:
    plt.plot(z_arr_ext,n_arr_ext)
    plt.grid(); plt.show()

dump = 10
cores = 4

ramp_end = z_arr[118] + z_offset
endindex = np.nonzero(z_arr_ext>ramp_end)[0][0]

betastar = 0.00213065326633
waist_loc = -0.000206030150754

z_arr = z_arr_ext
n_arr = n_arr_ext

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                               beta_offset=waist_loc, plasma_start=z_offset)
beam = PProp.GaussianBeam(beam_params, debug)

argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                               nset = max(n_arr))
argon_params['Z'] = z_arr[-1]*1e6
argon_params['Nz']= len(z_arr)
argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
argon = PProp.CustomPlasma(argon_params, n_arr, debug)

z_fine = np.copy(z_arr_ext)*1e6#np.linspace(z_arr[0], z_arr[-1], 3000)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)

PProp.PlotEmittance(beam,z_fine,m)
PProp.PlotGamma(beam,z_fine,m)

print("Bmag BP: ",PProp.GetBmag(beam,m))
bmagc = PProp.Calc_Bmag(beam_params,n_arr[:endindex], z_arr[:endindex])
print("Bmag CS: ",bmagc)

PProp.Plot_CSEvo(beam_params, n_arr, z_arr, z_offset)
"""
bmag2 = PProp.Calc_Bmag(beam_params,n_arr, z_arr)
print("Bmag2: ",bmag2)
"""
"""
z_finem = z_fine/1e6
spl = interp.InterpolatedUnivariateSpline(z_arr, n_arr)
n_fine = spl(z_finem)
bmag3 = PProp.Calc_Bmag(beam_params,n_fine, z_finem)
print("Bmag3: ",bmag3)
"""