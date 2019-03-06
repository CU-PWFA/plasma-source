#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 14:16:21 2018

First a copy of Keenan's FieldIonize_v02.ipynb, then I want to basically plot
beam density vs beam ionization to find the upper limit of safe nonlinear PWFA

@author: chris
"""

import sys
sys.path.insert(0,"../")
import Constants.SI as SI
import eBeam_v02 as eb
import numpy as np
from scipy.integrate import simps
import plotty as pl
from importlib import reload as rel
c = SI.lightSpeed;
import matplotlib.pyplot as plt

nlen = 50
q_arr = np.linspace(1e-9, 2e-9, nlen)
np_arr = np.zeros(nlen)

for i in range(len(q_arr)):
    beamParams = eb.get_beam()
    beamParams['charge'] = q_arr[i]
    beamParams['sigma_z'] = 36e-6
    pos = eb.get_pos(beamParams, nxi = 5)
    
    Er, rPeak, EPeak = eb.rad_E_field(pos, beamParams)
    
    W_Ar1 = eb.ionization_rate(Er, 'Ar+')
    W_Ar2 = eb.ionization_rate(Er,  'Ar2+')
    W_He1 = eb.ionization_rate(Er, 'He+')
    W_He2 = eb.ionization_rate(Er,  'He2+')
    
    ion_frac_Ar1, max_ion_Ar1 = eb.ionization_frac(W_Ar1, pos, beamParams)
    ion_frac_Ar2, max_ion_Ar2 = eb.ionization_frac(W_Ar2, pos, beamParams)
    ion_frac_He1, max_ion_He1 = eb.ionization_frac(W_He1, pos, beamParams)
    ion_frac_He2, max_ion_He2 = eb.ionization_frac(W_He2, pos, beamParams)
    
    width_Ar1 = eb.neutral_ring_width(ion_frac_Ar1, pos)
    width_Ar2 = eb.neutral_ring_width(ion_frac_Ar2, pos)
    width_He1 = eb.neutral_ring_width(ion_frac_He1, pos)
    width_He2 = eb.neutral_ring_width(ion_frac_He2, pos)
    widths = [width_Ar1, width_Ar2, width_He1, width_He2]
    
    max_frac = [max_ion_Ar1, max_ion_Ar2, max_ion_He1, max_ion_He2]
    plasmaNames = ['Ar+', 'Ar2+', 'He+', 'He2+']
    
    rel(pl)
    
    #pl.plot_field(Er, pos, beamParams, '$|E_r|$ [GV/m]', ind = 0)
    Ermax = np.amax(Er)
    
    #pl.plot_field(W_Ar1, pos, beamParams, 'W [fs$^{-1}$]', \
    #      ind = 0, gas = True, gasName = 'Argon')
    Wmax = np.amax(W_Ar1)
    
    rel(pl)
    nrz = pl.plot_2D_plasma(W_Ar1, pos, beamParams, 'Singly ionized argon', ind = 0)
    nrzmax = np.amax(nrz)
    np_arr[i] = nrzmax
    
plt.plot(q_arr*1e9, np_arr)
plt.xlabel(r'Charge [$nC$]')
plt.ylabel("Peak Ionization Frac.")
plt.grid(); plt.show()

den_arr = np.zeros(nlen)
e = 1.6022e-19
sigz = beamParams['sigma_z']
sigr = beamParams['sigma_r']
den_arr = q_arr/e/np.power(2*np.pi,3/2)/(sigz*sigr*sigr)/np.power(100,3)

plt.plot(den_arr, np_arr)
plt.xlabel(r'Density [$cm^{-3}$]')
plt.ylabel("Peak Ionization Frac.")
plt.grid(); plt.show()