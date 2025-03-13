#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:27:48 2018

Quick script to take a double tanh profile and return the focal length

@author: chris
"""

import sys
sys.path.insert(0, "../")
import numpy as np
import matplotlib.pyplot as plt

from modules import ThreeDimensionAnalysis as ThrDim
from modules import TPLFocalLength as Foc

#from scipy.signal import find_peaks, savgol_filter

gam = Foc.gam_def

"""
a = 86.8450398321
b = 19.8872039701
n_0 = 0.502043587936
"""
a = 35.9713173965
b = 2.67654175646*4
n_0 = 5.09212414537# / 1e17
"""
a = 150
b = 20
n_0 = 0.2# / 1e17
"""

a = 40
b = 2.67654175646*3
n_0 = 1# / 1e17

fit = [a,b,n_0]
edge = 1/2*(3*a+6*b)
z_arr = np.linspace(-edge, edge, 201)
n_arr = ThrDim.DoubleTanh(fit,z_arr)

focal = Foc.Calc_Focus(n_arr, z_arr, gam)
thick = Foc.Calc_Square_Lens(n_0*1e17, focal, gam)

plt.plot(z_arr, n_arr,c='b')
plt.plot([-edge,-0.5*thick],[0,0],c='r',ls='--')
plt.plot([-0.5*thick,-0.5*thick],[0,n_0],c='r',ls='--')
plt.plot([-0.5*thick,0.5*thick],[n_0,n_0],c='r',ls='--')
plt.plot([0.5*thick,0.5*thick],[n_0,0],c='r',ls='--')
plt.plot([0.5*thick,edge],[0,0],c='r',ls='--')
#plt.title("Density along beam axis")
plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
plt.grid(); plt.show()

print(thick,"um equivalent thickness")