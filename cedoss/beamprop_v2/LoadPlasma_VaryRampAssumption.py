#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 14:06:12 2018

Assuming the ramp is matched for a ramp half-length of .5mm, vary this length
and measure resultant B-mag.  Also look at optimizing z_tpl, and later
z_tpl = z_waist.  Should give framwork for experimental tolerance of
unknown low density regions

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
z_arr_orig = np.array(sheet['Length [m]'])
n_arr_orig = np.array(sheet['Density [cm-3]'])/1e17

maxdist = 0.035

halflen_arr = np.linspace(0.0001, 0.0010, 50)
tplloc_arr = np.zeros(len(halflen_arr))
bmag_arr = np.zeros(len(halflen_arr))

for i in range(len(halflen_arr)):
    if i%5==0:
        print(i/len(halflen_arr)*100,"%")
        
    extend = halflen_arr[i]#0.0005
    num_iter = int(maxdist/extend)
    n_arr_add = np.copy(n_arr_orig); z_arr_add = np.copy(z_arr_orig)
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
    
    ramp_end = z_arr_orig[118] + z_offset
    endindex = np.nonzero(z_arr_ext>ramp_end)[0][0]
    
    tpl_n = 10.
    tpl_l = 74.91
    tpl_f = 0.01475
    
    betastar = .10 #0.00213065326633
    #waist_loc = -0.000206030150754 - tpl_f
    
    #tpl_offset = waist_loc
    #off_arr = np.linspace(-maxdist, 0, 101)
    off_arr = np.linspace(-0.018, -0.014, 101)
    bmag_inner = np.zeros(len(off_arr))
    for j in range(len(off_arr)):
        z_arr = np.copy(z_arr_ext)
        n_arr = np.copy(n_arr_ext)
                
        tpl_offset = off_arr[j]
        waist_loc = tpl_offset
                
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                        beta_offset=waist_loc, plasma_start=z_offset)
                
        argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = max(n_arr))
        argon_params['Z'] = z_arr[-1]*1e6; argon_params['Nz']= len(z_arr)
        argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
        argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6, tpl_n, tpl_l, debug)
                
        bmag = PProp.Calc_Bmag(beam_params,n_arr[:endindex], z_arr[:endindex])
        bmag_inner[j] = bmag
    
    minloc = np.argmin(bmag_inner)
    bmag_arr[i]=bmag_inner[minloc]
    tplloc_arr[i]=off_arr[minloc]

plt.plot(halflen_arr*1000, bmag_arr, label = "B-mag")
plt.plot([halflen_arr[0]*1000,halflen_arr[-1]*1000],[1.01,1.01], label="1.01")
plt.title("B-mag vs ramp half-length")
plt.xlabel("Half-length [mm]")
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()

plt.plot(halflen_arr*1000, tplloc_arr*100, label = "B-mag")
plt.title("Optimal lens location vs ramp half-length")
plt.xlabel("Half-length [mm]")
plt.ylabel(r'$z_{TPL}\,\mathrm{[cm]}$')
plt.grid(); plt.show()