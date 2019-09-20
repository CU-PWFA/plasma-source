#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 11:38:05 2018

Propagate a beam through free space after loading from a h5 file

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit

path = '/home/chris/Desktop/BeamProp/vorpaltest'
filename = '/home/chris/Desktop/FACETII_NERSC_Run4/PTPLDoubleTanh_WitnessBeam_10.h5'#z_arr to .03, threshold .00001
#filename = '/home/chris/Desktop/SFQED_NERSC_VaryDen/5e15/PTPLDoubleTanh_ElectronBeam_10.h5' #z_arr to 0.7, ind 41, thresh 0.001

debug = 1

dump = 20
cores = 4
threshold = 0.00001#0.00001

#beam = PProp.VorpalBeam(path, filename, threshold, minz=-8e-6, maxz = -7e-6, debug=debug)
beam = PProp.VorpalBeam(path, filename, threshold, debug=debug)
print("N: ",beam.N)

z_arr = np.linspace(0, 0.03, 1000 + 1)
n_arr = np.zeros(len(z_arr))#+1e-9

argon_params = PProp.ReturnDefaultPlasmaParams(path)
argon_params['Z']=z_arr[-1]*1e6
argon_params['z0']=z_arr[0]*1e6
argon_params['Nz']=len(z_arr)
argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, 0,0,0,0)

z_fine = np.copy(z_arr)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)
#PProp.PlotEmittance(beam,z_fine,m)
sig_arr, s_arr = PProp.PlotSigmar(beam,z_fine,m)
#print("Bmag: ",PProp.GetBmag(beam,m))

minloc = np.argmin(sig_arr)
print("index: ", minloc, "|  size [um]: ", sig_arr[minloc])
#beam.plot_phase_at(minloc)

#minloc=41#Used for sfqed 5e15 with linspace 0,0.7,1001 for smallest inner waist
print(beam.get_sigmar_frac(minloc, 58.82e-6))#90% of beam at ind 41 for 5e15

beam.plot_hist_at(minloc)
#beam.plot_phase_hist_at(minloc,fitted=True)
