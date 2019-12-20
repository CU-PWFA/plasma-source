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
filename = '/home/chris/Desktop/thick4/PTPLDoubleTanh_rhoPlasma_3.h5'

#data = load.get_field_data(filename,'rhoWitness')
#data = load.get_field_data(filename,'rhoDrive')
data = load.get_field_data(filename,'rhoPlasma')
data = np.array(data)

total_rho = np.sum(data)

#Manually input for now, add other stuff later to auto this if neded

#thick4 AND #cylinder
dx = 5.590551181102362e-07; dy = 5.597014925373134e-07; dz = 5.597014925373134e-07;

total_charge = total_rho * dx * dy * dz
print("Charge = ",total_charge*1e9,"nC")