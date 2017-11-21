#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:06:07 2017

New and improved looper which utilizes the new and improved module version
of FourierRefraction.  Set up the params file to include pulseParams
and then pass to FourierRefraction.py in modules

@author: chris
"""
import sys
import numpy as np
sys.path.insert(0, "../")
from modules import FourierRefraction as frefract

#Script Selection
test_density = 1
free_space = 0
double_x = 1

#Set to some integer to lower grid resolution by a power
reducer = 0

load_density = 1
constant_density = .1

#DONT FORGET TO CHECK THE DIRECTORY AND ITS STUFF!
# initDensity only knows its shape, not dimensions
denfile = '/home/chris/Desktop/DataLoads/DensityFilesNp/initDensity_Ar_6x15x150.npy'
pulsefile = '/home/chris/Desktop/DataLoads/PulseFilesNp/pulseParams_s4_center.npy'

path = '/home/chris/Desktop/FourierPlots/ArJets_DivideDensity/'

params_orig = {'Nx' : 2**(9-reducer)  * (1+double_x),
               'Ny' : 2**(9-reducer),
               'Nz' : 2**(8-reducer),
               'Nt' : 2**7,#7
               'X' : 3e2 * (1+double_x),
               'Y' : 15e2,
               'Z' : 1.5e4,
               'T' : 160,
               'n0' : 1,
               'alpha' : 1.63,
               'EI' : 15.8,
               'E0' : 1,
               'lam' : 0.7853,
               'n' : 1.0,
               'nFinal' : True,
              }

var_loop = [1]#, 2, 3, 5, 8, 10, 15, 20]
for case in var_loop:
    
    """####Set Looped Var####"""
    division = case
    """######################"""
    constant_density=0.1
    
    params = params_orig
    params = frefract.LoadPulseParams(pulsefile,params)
    directory = 'case_' + str(case)
    params['path'] = path+directory+'/'
    params['d_denfile'] = denfile
    params['d_pulsefile'] = pulsefile
    
    density=None
    if load_density == 1:
        density = np.load(denfile)
    elif constant_density != 0:
        params['n0'] = constant_density
    density=density/division
    
    if test_density == 1:
        frefract.TestDensity(density,params)
    elif free_space == 1:
        frefract.RunFreeSpace(params)
    else:
        frefract.RunRefraction(density,params)