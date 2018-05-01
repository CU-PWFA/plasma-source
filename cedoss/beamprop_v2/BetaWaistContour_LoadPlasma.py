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
sheet = xl.parse('Density Data')
z_arr = np.array(sheet['Length [m]'])
n_arr = np.array(sheet['Density [cm-3]'])/1e17

#Recursively extend the data by setting density to half every extension step
extend = 0.0010; num_iter = 20 # DEFAULT:  extend = 0.0005; num_iter = 20
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
z_offset = -z_arr_ext[0]
z_arr_ext = z_arr_ext + z_offset
if debug == 1:
    plt.plot(z_arr_ext,n_arr_ext)
    plt.grid(); plt.show()

dump = 10
cores = 4

ramp_end = z_arr[118] + z_offset
z = z_arr_ext
n = n_arr_ext
endindex = np.nonzero(z>ramp_end)[0][0]

beta_arr = np.linspace(0.003, 0.007, 200)
waist_arr = np.linspace(-0.005,-0.001, 200)
#beta_arr = np.linspace(0.001, 0.004, 50) #DEFAULT
#waist_arr = np.linspace(-0.001,0.001, 50)
bmag_image = np.zeros((len(beta_arr),len(waist_arr)))
for i in range(len(beta_arr)):
    if i%10 == 0: print(i/len(beta_arr)*100,"%");
    for j in range(len(waist_arr)):
        betastar = beta_arr[i]
        waist_loc = waist_arr[j]
        
        #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset)
        
        bmag = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #bmag = PProp.Calc_Bmag(beam_params,n, z)
        bmag_image[i][j] = bmag

minloc = PProp.PlotContour(bmag_image, beta_arr, waist_arr, r'$\beta^*$ [m]', r'$z_{v}$ [m]')