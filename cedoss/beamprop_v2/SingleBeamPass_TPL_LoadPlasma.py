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

debug = 1

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
z_arr_ext = np.linspace(z_arr_add[0], z_arr_add[-1], 80000)
n_arr_ext = spl(z_arr_ext)
z_offset = -z_arr_ext[0]
z_arr_ext = z_arr_ext + z_offset
#print(n_arr_ext[740])
if debug == 1:
    plt.scatter(z_arr_add[40:60]*100,n_arr_add[40:60],s=10,c='r',label="Data")
    plt.scatter(z_arr_add[:40]*100,n_arr_add[:40],s=10,c='b',label="Artificial Ramp")
    plt.plot((z_arr_ext-z_offset)*100,n_arr_ext,label="Spline Fit")
    plt.ylim([-1,25])
    plt.xlim([-0.006*100,0.001*100])
    plt.title("Spline fit to artificial density ramp")
    plt.xlabel("z [cm]")
    plt.ylabel(r'$n_p\,\mathrm{[10^{17}cm^{-3}]}$')
    plt.legend(); plt.grid(); plt.show()
    
    plt.scatter(z_arr_add[130:170]*100,n_arr_add[130:170],s=10,c='r',label="Data")
    plt.scatter(z_arr_add[158]*100,n_arr_add[158],s=40,c='g',label="Ramp Top")
    plt.plot((z_arr_ext-z_offset)*100,n_arr_ext,label="Spline Fit")
    plt.ylim([1050,1250])
    plt.xlim([0.0045*100,0.0065*100])
    plt.title("Spline fit at ramp peak")
    plt.xlabel("z [cm]")
    plt.ylabel(r'$n_p\,\mathrm{[10^{17}cm^{-3}]}$')
    plt.legend(); plt.grid(); plt.show()

dump = 10
cores = 4

ramp_end = z_arr[118] + z_offset
endindex = np.nonzero(z_arr_ext>ramp_end)[0][0]

tpl_n = 10.
tpl_l = 74.91
tpl_f = 0.01475

betastar = .10 #0.00213065326633
waist_loc = -0.000206030150754 - tpl_f

tpl_offset = waist_loc

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
argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6, tpl_n, tpl_l, debug)

z_fine = np.copy(z_arr_ext)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)

PProp.PlotEmittance(beam,z_fine,m)
PProp.PlotGamma(beam,z_fine,m)

print("Bmag BP: ",PProp.GetBmag(beam,m))
bmagc = PProp.Calc_Bmag(beam_params,n_arr[:endindex], z_arr[:endindex])
print("Bmag CS: ",bmagc)

PProp.Plot_CSEvo(beam_params, n_arr, z_arr, z_offset)