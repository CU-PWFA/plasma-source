#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:23:56 2018

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np

path = '/home/chris/Desktop/BeamProp/test'

debug = 2 #0 for nothing, 1 for everything, 2 for minimal updates

dump = 1000
cores = 4
tpl_n = 0.5

offset_arr = np.linspace(-0.60e6, -0.20e6, num = 20)
length_arr = np.linspace(0, 300, num = 20)
bmag_image = np.zeros((len(offset_arr),len(length_arr)))

for i in range(len(offset_arr)):
    for j in range(len(length_arr)):
        tpl_offset = offset_arr[i]
        tpl_l = int(length_arr[j])
        
        beam_params = PProp.ReturnDefaultElectronParams(path)
        beam = PProp.GaussianBeam(beam_params, debug)
        
        argon_params = PProp.ReturnDefaultPlasmaParams(path)
        argon = PProp.GaussianRampPlasma_ThinPlasmaLens(argon_params, tpl_offset, tpl_n, tpl_l, debug)
        
        z_orig = np.linspace(0,argon_params['Z'],int(((argon_params['Nz']-1)/5)+1))
        z_fine = PProp.FineSpacingLens(z_orig, argon_params['z0'] + tpl_offset, tpl_l)
        PProp.PropBeamPlasma(beam,argon,z_fine,dump,cores, debug)
        
        m = int(len(z_fine)/dump)
        
        bmag_image[i][j] = PProp.GetBmag(beam,m)
        if debug == 2:
            print(offset_arr[i]," -- ", length_arr[j], " -- ", bmag_image[i][j])
        
minloc = PProp.PlotContour(bmag_image, offset_arr/1e6, length_arr, r'$z_{TPL}$ [m]', r'$TPL_{\rm L}$ [um]')