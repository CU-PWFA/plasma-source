#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:43:44 2023

New and improved looper which utilizes the new and improved module version
of FourierRefraction.  Set up the params file to include pulseParams
and then pass to FourierRefraction.py in modules

Cleaned up and commented

@author: chris
"""
import sys
import numpy as np
sys.path.insert(0, "../")
from modules import FourierRefraction as frefract

#Script Selection
test_density = False    #Plot the density file to check it is correct
free_space = False      #Propagate pulse through free space to check it is correct
double_x = 0            #Set to 0 for normal, set to 1 for Nx to be 2*Ny

#Set to some integer to lower grid resolution by a power
reducer = 0             #Set between 0-5 or so.  Set to 0 for high resolution production runs.

load_density = False    #True for loading 3D numpy density from file.  False for using a constant density
constant_density = 0.01 #If load_density = False, then the gas is uniform at this density (units of 10^17 cm^-3)

#DONT FORGET TO CHECK THE DIRECTORY
# initDensity only knows its shape, not dimensions
denfile = '/home/chris/Desktop/DataLoads/DensityFilesNp/gasjet_He1e20_20x10x1200/gasdensity.npy' #Location of the 3D numpy density profile, if load_density = True
pulsefile = '/home/chris/Desktop/DataLoads/PulseFilesNp/pulseParams_cylindrical_ex_thesis2.npy'  #Location of the Gaussian pulse file
outputpath = '/home/chris/Desktop/FourierPlots/outputfolder/'                                    #Location of the output of Fourier simulation

#Generate a params dictionary with the default values.
# double_x can make Nx twice as large as Ny if double_x = 1
# reducer is an easy way to lower the computational time.  Set to 0 for production runs.
"""
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
"""
params_orig = frefract.GetDefaultParams(reducer,double_x)

#Make any necessary adjustments to the default parameters.
#Typically, I change 'X' and 'Y' to be as large as they need to be to fully contain the electric field
#           'Z' needs to be changed to match 2*zi from QParamProp (there is a factor of 1e6 here as well)
#           and 'EI' should be set to the ionization energy for whatever gas we are using
params_orig['X']=20e2#50e2
params_orig['Y']=10e2#50e2
params_orig['Z']=12.0e4#7.0e4
params_orig['EI']=15.4
params_orig['lam']=0.800
params_orig['T']=230

#Can make a loop of many cases, or just run a single instance
var_loop = [1]#, 2, 3, 5, 8, 10, 15, 20]
for case in var_loop:
    
    """####Set Looped Var####"""
    #division = case
    division = 1
    """######################"""
    
    params = params_orig
    params = frefract.LoadPulseParams(pulsefile,params)
    directory = 'case_' + str(case)
    params['path'] = outputpath+directory+'/'
    params['d_denfile'] = denfile
    params['d_pulsefile'] = pulsefile
    
    density=None
    if load_density:
        density = np.load(denfile)
        density=density/division
    elif constant_density != 0:
        params['n0'] = constant_density/division
    
    
    if test_density:
        frefract.TestDensity(density,params)
    elif free_space:
        frefract.RunFreeSpace(params)
    else:
        frefract.RunRefraction(density,params)