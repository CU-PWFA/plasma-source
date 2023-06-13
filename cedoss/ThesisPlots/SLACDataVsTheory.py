#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:06:53 2023

Quick Plots of Theory vs Experiment for the SLAC Data.

Its gonna look bad, but should be in.

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

#Full Projection Fits, as the e slices dont provide too much
fit_beta_arr =  np.array([151.5,104.5,84.19,80.12,55.56,63.11]) #cm
fit_z0_arr =    np.array([-0.3,-.2,0.077,0.086,-0.034,0.13])#m
fit_emitn_arr = np.array([42.81,71.02,128.63,106.79,190.83,190.46])#mm-mrad

#The datasets and their densities
psi_arr = np.array([0,1,6,25,57.8,115.8])
den_arr = np.array([0,0.27,1.62,6.48,15.3,31.5])#e16 cm-3

#Initial Beam Parameters



#Approximate Plasma Lens Profile


