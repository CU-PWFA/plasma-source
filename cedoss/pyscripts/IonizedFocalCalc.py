#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 14:27:17 2017

Given some narrow and wide waist, and some initial power/laser size, 
calculate the focal length of an ionized plasma lens assuming no
refraction and all Gaussian/thin assumptions are held

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../../python")
from ionization import ionization
from ionization import adk

sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

#pow0 = 50e9
pow0 = 0.5e12
w0 = 5e-3
lmda = 785.3e-9

#wx0 = 25.0801912297e-6
wx0 = 25.0801912297e-5
wy = 226.156532851e-6

pulsetime = 35
ionization_chi = 15.8

back_den = 0.1

#def wx(z):
#    return wx0*np.sqrt(1+np.square((z*lmda)/(np.pi*wx0**2)))

def PowerToIntensity(pow_init):
    return 2*pow_init/np.pi/np.power(w0,2)

def GetDomain(w_long):
    return np.linspace(-4*w_long, 4*w_long, 1000)

def Intensity_1D(pow_init, w_init, w_long, w_trans):
    xrange = GetDomain(w_long)
    intensity_init = PowerToIntensity(pow_init)
    return intensity_init*w_init**2/(w_long*w_trans)*np.exp(-2*np.square(xrange/w_long))

def Ionization_1D(intensity_arr, chi_energy, delt_t):
    E = ionization.field_from_intensity(intensity_arr)
    return adk.gaussian_frac(chi_energy,E,delt_t,1)

x = GetDomain(wx0) * 1e6
I = Intensity_1D(pow0, w0, wx0, wy) / 1e14 / 1e4
H = Ionization_1D(I, ionization_chi, pulsetime)
den = H*back_den

plt.plot(x, I)
plt.show()
plt.plot(x, den)
plt.show()
Foc.Calc_Focus(den,x)