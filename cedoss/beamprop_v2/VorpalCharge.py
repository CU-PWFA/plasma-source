#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:25:30 2019

The goal is to make a script which reads a rho file and outputs a total charge in nC

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit

import sys
sys.path.insert(0, "../../python")

from vsim import load

#filename = '/home/chris/Desktop/thick4/PTPLDoubleTanh_rhoDrive_10.h5'
filename = '/home/chris/Desktop/200e16/PTPLDoubleTanh_rhoPlasma_5.h5'
#filename = '/home/chris/Desktop/afterglow_run1/PTPLDoubleTanh_rhoDrive_10.h5'

#data = load.get_field_data(filename,'rhoWitness')
#data = load.get_field_data(filename,'rhoDrive')
data = load.get_field_data(filename,'rhoPlasma')
data = np.array(data)

total_rho = np.sum(data)

#Manually input for now, add other stuff later to auto this if neded

#thick4 AND #cylinder
#nowafterglow
dx = 5.629629629629629e-07; dy = 5.639097744360901e-07; dz = 5.639097744360901e-07;

total_charge = total_rho * dx * dy * dz
print("Charge = ",total_charge*1e9,"nC")