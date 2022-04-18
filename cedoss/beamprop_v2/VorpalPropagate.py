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

import scipy.integrate as Int
def Gauss(p, x):
    return p[1]*np.exp(-.5*np.square(x)/np.square(p[0]))
def GaussPlusGauss_Percent(p):
    pi = [np.abs(p[0]), np.abs(p[2])]
    po = [7.54e-6, np.abs(p[3])]
    inner = Int.quad(lambda x: Gauss(pi, x), -5*pi[0], 5*pi[0])[0]
    outer = Int.quad(lambda x: Gauss(po, x), -5*po[0], 5*po[0])[0]
    print(inner, outer)
    print("Inner: ",inner/(inner+outer)*100,"%")
    print("Outer: ",outer/(inner+outer)*100,"%")
    return  


path = '/home/chris/Desktop/BeamProp/vorpaltest'
#filename = '/home/chris/Desktop/FACETII_NERSC_Run4/PTPLDoubleTanh_WitnessBeam_10.h5'#z_arr to .03, threshold .00001
#filename = '/home/chris/Desktop/SFQED_NERSC_VaryDen/5e15/PTPLDoubleTanh_ElectronBeam_10.h5' #z_arr to 0.7, ind 41, thresh 0.001
#filename = '/home/chris/Desktop/thick4/PTPLDoubleTanh_WitnessBeam_10.h5'

filename = '/home/chris/Desktop/100e16/PTPLDoubleTanh_WitnessBeam_10.h5'

filename = '/home/chris/Desktop/NERSC_Deflection_July/PTPL_Gradient_WitnessBeam_10.h5'

filename = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_Deflection_190um/PTPL_Gradient_WitnessBeam_10.h5'
#filename = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/Feb_190um_separation/PTPL_Gradient_WitnessBeam_10.h5'
#filename = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_Control_114um/PTPL_Gradient_WitnessBeam_10.h5'

filename = '/home/chris/Desktop/singlebunch/PTPL_Gradient_ElectronBeam_20.h5' #z_arr to 0.08, ind 41, thresh 0.001

debug = 1

dump = 20
cores = 4
threshold = 0.001#0.00001

#beam = PProp.VorpalBeam(path, filename, threshold, minz=-8e-6, maxz = -7e-6, debug=debug)
beam = PProp.VorpalBeam(path, filename, threshold, debug=debug)
print("N: ",beam.N)

#z_arr = np.linspace(0, 0.03, 1000 + 1)
#z_arr = np.linspace(0, 0.15, 1000 + 1)
z_arr = np.linspace(0, 0.08, 1000 + 1)
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
print("index: ", minloc, "|  position [cm]: ", s_arr[minloc],"|  size [um]: ", sig_arr[minloc])
#beam.plot_phase_at(minloc)

#minloc=41#Used for sfqed 5e15 with linspace 0,0.7,1001 for smallest inner waist
print(beam.get_sigmar_frac(minloc, 58.82e-6))#90% of beam at ind 41 for 5e15

beam.plot_hist_at(minloc)
GaussPlusGauss_Percent([6.54905155e-07,   6.65681022e-06,   8.06943422e+04,   3.06225426e+03])
#beam.plot_phase_hist_at(minloc,fitted=True)
