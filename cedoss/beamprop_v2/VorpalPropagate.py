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
filename = '/media/chris/New Volume/VSimRuns/NERSCResults/3D/TwoBunch/1e18/200um/ThinPlasmaLens3D_WitnessBeam_10.h5'

debug = 1

dump = 20
cores = 4
threshold = 0.001

#beam = PProp.VorpalBeam(path, filename, threshold, minz=-8e-6, maxz = -7e-6, debug=debug)
beam = PProp.VorpalBeam(path, filename, threshold, debug=debug)
print("N: ",beam.N)

z_arr = np.linspace(0, 0.006, 1000 + 1)
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

beam.plot_hist_at(minloc)
