#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:45:17 2019

@author: mike
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.constants as const
me = const.physical_constants['electron mass energy equivalent in MeV'][0]

# WARGSim Cython code
sys.path.insert(0, "/home/mike/work/github_repos/CU-PWFA/plasma-source/python")
from beam.beams import electronbeam
from beam.elements import plasma_1d as plasma
from beam import interactions
import beam.calc.electron as ecalc


# define electron beam
path = "./"
N      = 10000 # number of particles
gamma  = (10e9)/me # reltivistic Lorentz factor (centroid value)
emit   = 5.0e-6  # m-rad, normalized emittance
betax  = 0.10 # m, x beta function at z=0
betay  = 0.10 # m, y beta function at z=0
alphax = 0.00 # x alpha function at z=0
alphay = 0.00 # y alpha function at z=0
sigmaz = 5e-6 # m, rms bunch length
dE     = 1.0e-2 # relative energy spread

electronParams = {
        'name'      : 'TestBeam',
        'path'      : path,
        'load'      : False,
        'N'         : N,         #10000 for normal, 1000000 for production
        'gamma'     : gamma,
        'emittance' : emit,
        'betax'     : betax,
        'betay'     : betay,
        'alphax'    : alphax,
        'alphay'    : alphay,
        'sigmaz'    : sigmaz,
        'dE'        : dE
    }    
beam = electronbeam.GaussianElectronBeam(electronParams)



# propagate the beam through the plasma
""" Propagate an electron beam through an ion column. 

Parameters
----------
electron : ElectronBeam class
    The electron beam to propagate through the plasma.
plasma : Plasma class
    The plasma that the electron beam will be propagating through.
z : array-like
    The spatial grid used to set the step size for the electron beam.
dumpPeriod : int
    How frequently to save the electron beam to disk.
n : int
    Number of threads to run on.
"""
#interactions.electron_plasma(beam, plasma, z_arr, dumpPeriod, cores)