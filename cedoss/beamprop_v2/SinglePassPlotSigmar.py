#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:23:56 2018

Single pass through plasma with optional TPL using Robert's code


@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit

path = '/home/chris/Desktop/BeamProp/test'

debug = 1
vac = 1

dump = 100
cores = 4

Z_set = 0.25*1e6
#Z_set = PProp.def_Z

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path)
beam = PProp.GaussianBeam(beam_params, debug)

argon_params = PProp.ReturnDefaultPlasmaParams(path, Z_change = Z_set)
argon = PProp.GaussianRampPlasma(argon_params, debug)

z_arr = np.linspace(0,argon_params['Z'],int(((argon_params['Nz']-1)/5)+1))

start = timeit.default_timer()
PProp.PropBeamPlasma(beam, argon, z_arr, dump, cores, debug)
print("prop time: ",str(timeit.default_timer()-start))

m = int(len(z_arr)/dump)
PProp.PlotEmittance(beam, z_arr, m)
PProp.PlotSigmar(beam, z_arr, m)

print("Bmag: ",PProp.GetBmag(beam,m))

if vac ==1:
    beam_params2 = beam_params
    beam_params2['path'] = path+'/vac'
    beam2 = PProp.GaussianBeam(beam_params2, debug)
    
    argon_params2 = PProp.ReturnDefaultPlasmaParams(path+'/vac', Z_change = Z_set)
    argon_params2['n0']=0.0; argon_params2['dgammadz'] = PProp.dgammadz_basic
    argon2 = PProp.GaussianRampPlasma(argon_params2, debug)
    
    z_arr2 = np.linspace(0,argon_params2['Z'],int(((argon_params2['Nz']-1)/5)+1))
    
    start = timeit.default_timer()
    PProp.PropBeamPlasma(beam2, argon2, z_arr2, dump, cores, debug)
    print("prop time: ",str(timeit.default_timer()-start))
    
    m2 = int(len(z_arr2)/dump)
    PProp.PlotEmittance(beam2, z_arr2, m2)
    PProp.PlotSigmar(beam2, z_arr2, m2)
    
    print("Bmag: ",PProp.GetBmag(beam2,m2))
    PProp.PlotEmittance_Compare(beam, beam2,z_arr,m,"With Plasma","Without Plasma")
    PProp.PlotSigmar_Compare(beam, beam2,z_arr,m,"With Plasma","Without Plasma")