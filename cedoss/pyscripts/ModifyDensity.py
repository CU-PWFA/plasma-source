#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 11:18:05 2017

@author: chris
"""

import numpy as np

folder = '/home/chris/Desktop/DataLoads/DensityFilesNp/'
infile = 'initDensity_Ar_3x15x150.npy'
outfile = 'initDensity_Ar_6x15x150.npy'

script = 1

if script == 1:#Double x domain (narrow waist)
    background = 0.1
    den_old = np.load(folder+infile)
    shp = den_old.shape
    den_new = np.full((shp[0]*2, shp[1]*1, shp[2]*1), background)
    den_new[int(shp[0]*1/2):int(shp[0]*3/2),:,:]=den_old[:,:,:]
    np.save(folder+outfile,den_new)