#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 15:03:12 2018

Reads in a beam and calculates the emittance under different conditions

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit

path = '/home/chris/Desktop/BeamProp/vorpaltest'
filename = '/home/chris/Desktop/FACETII_NERSC_Run4/PTPLDoubleTanh_WitnessBeam_10.h5'
#filename = '/home/chris/Desktop/SFQED_NERSC_Run1/PTPLDoubleTanh_ElectronBeam_10.h5'
#filename = '/home/chris/Desktop/PTPL_Cylinder/PTPLCylinder_electrons_10.h5'
filename = '/home/chris/Desktop/thick4/PTPLDoubleTanh_WitnessBeam_10.h5'

debug = 0

dump = 20
cores = 4
threshold = 0.00001

beam = PProp.VorpalBeam(path, filename, threshold, debug=debug)
print("N: ",beam.N)

beam.plot_phase_hist_at(0, fitted=True)
print(beam.get_emittance_n(0)[0],beam.get_emittance_n(0)[1])
print("Full emittance: ",np.average(beam.get_emittance_n(0)))

#print("Front: ",np.average(beam.get_emittance_n_zcond(0, 0, 10000)))
#print("Back: ",np.average(beam.get_emittance_n_zcond(0, -10000, 0)))