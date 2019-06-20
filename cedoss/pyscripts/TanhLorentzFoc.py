#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:47:31 2019

Quick script to take a double tanh-Lorentzian profile and return the focal length

@author: chris
"""

import sys
sys.path.insert(0, "../")
import numpy as np
import matplotlib.pyplot as plt

from modules import ThreeDimensionAnalysis as ThrDim
from modules import TPLFocalLength as Foc

gam = 30 # Foc.gam_def

"""
a = 86.8450398321
b = 19.8872039701
n_0 = 0.502043587936
"""
"""
a = 10 * 1000
b = 4 * 1000
n_0 = 4.9e12 / 1e17
"""
a = 298.465129594
b = 78.1121221798
n_0 = 1.00251813855

A = 7.56e3
gamma = 4812.9
x_0 = 0

fit_tanh = [a,b,n_0]
fit_lorentz = [A, gamma, x_0]
edge = 1/2*(3*a+6*b)
z_arr = np.linspace(-edge, edge, 201)
n_arr = ThrDim.DoubleTanh_Lorentzian(fit_tanh, fit_lorentz, z_arr)

plt.plot(z_arr, n_arr)
plt.title("Density along beam axis")
plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
plt.grid(); plt.show()

focal = Foc.Calc_Focus(n_arr, z_arr, gam)
thick = Foc.Calc_Square_Lens(n_0*1e17, focal, gam)
print(thick,"um equivalent thickness")

#Below is to compare with TanhFoc.py
"""
plt.plot(z_arr, n_arr-n_arrb)
plt.title("Density difference after Lorentzian")
plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
plt.ylabel(r'Density Change$\mathrm{\ [10^{17}cm^{-3}]}$')
plt.grid(); plt.show()
"""