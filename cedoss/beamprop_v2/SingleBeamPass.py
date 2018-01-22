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

tpl_offset = -0.37e6
tpl_n = 0.5
#tpl_l = 267
tpl_l = 0

dump = 1000
cores = 4

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path)
beam = PProp.GaussianBeam(beam_params, debug)

argon_params = PProp.ReturnDefaultPlasmaParams(path)
argon = PProp.GaussianRampPlasma_ThinPlasmaLens(argon_params, tpl_offset, tpl_n, tpl_l, debug)

z_orig = np.linspace(0,argon_params['Z'],int(((argon_params['Nz']-1)/5)+1))
z_fine = PProp.FineSpacingLens(z_orig, argon_params['z0'] + tpl_offset, tpl_l)

start = timeit.default_timer()
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)
print("prop time: ",str(timeit.default_timer()-start))

m = int(len(z_fine)/dump)
PProp.PlotEmittance(beam,m)

print("Bmag: ",PProp.GetBmag(beam,m))