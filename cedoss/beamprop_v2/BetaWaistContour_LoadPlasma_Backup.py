#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:27:21 2018

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit
import sys
import pandas as pd
import scipy.interpolate as interp

path = '/home/chris/Desktop/BeamProp/GasCellTest'
excel_filename = '/home/chris/Desktop/Density for tall gas cell lower pinhole .8 density 5.25 us.xlsx'

debug = 0

#Load density and distances from excel file
xl = pd.ExcelFile(excel_filename)
print(xl.sheet_names)
sheet = xl.parse('Density Data')
z_arr = np.array(sheet['Length [m]'])
n_arr = np.array(sheet['Density [cm-3]'])/1e17

#Recursively extend the data by setting density to half every extension step
extend = 0.0005; num_iter = 10
n_arr_add = np.copy(n_arr); z_arr_add = np.copy(z_arr)

for x in range(num_iter):
    n_arr_add = np.concatenate(([.5*n_arr_add[0]],n_arr_add,[.5*n_arr_add[-1]]))
    z_arr_add = np.concatenate(([z_arr_add[0]-extend],z_arr_add,[z_arr_add[-1]+extend]))
if debug == 1:
    plt.plot(z_arr_add[:25],n_arr_add[:25])
    plt.grid(); plt.show()
#Spline our extension into a new density and position array
spl = interp.InterpolatedUnivariateSpline(z_arr_add, n_arr_add)
z_arr_ext = np.linspace(z_arr_add[0], z_arr_add[-1], 5000)
n_arr_ext = spl(z_arr_ext)
z_arr_ext = z_arr_ext-z_arr_ext[0]
if debug == 1:
    plt.plot(z_arr_ext,n_arr_ext)
    plt.grid(); plt.show()

dump = 10
cores = 4

doss_addition = 12
rampend_i = 96 + doss_addition

start_loc = z_arr[rampend_i]

beta_arr = np.linspace(0.001, 0.01, 50)
waist_arr = np.linspace(-0.001, 0.001, 50)

bmag_image = np.zeros((len(beta_arr),len(waist_arr)))
for i in range(len(beta_arr)):
    #print(i/len(beta_arr),"%")
    for j in range(len(waist_arr)):
        betastar = beta_arr[i]
        waist_loc = waist_arr[j] - doss_addition*(z_arr[1]-z_arr[0])
        z = np.copy(z_arr_ext)
        n = np.copy(n_arr_ext)
        
        #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=start_loc)
        endindex = np.nonzero(z>start_loc)[0][0]
        bmag = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #bmag = PProp.Calc_Bmag(beam_params,n[:-1000], z[:-1000])
        #bmag = PProp.Calc_Bmag(beam_params,n, z)
        #print("Bmag CS: ",bmag)
        bmag_image[i][j] = bmag

minloc = PProp.PlotContour(bmag_image, beta_arr, waist_arr, r'$\beta^*$ [m]', r'$z_{v}$ [m]')