#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:52:57 2018
Quick script to take a Gaussian profile and return the focal length
@author: chris
"""
import sys
sys.path.insert(0, "../")
import numpy as np
import matplotlib.pyplot as plt
from modules import ThreeDimensionAnalysis as ThrDim
from modules import TPLFocalLength as Foc

gam = Foc.gam_def

n_0 = 3.82965364484e+16
sig = 1111.05229343
x_0 = -1.09058401518e-09

n_0=n_0/1e17
fit = [n_0, sig, x_0]
edge = 6*sig

z_arr = np.linspace(-edge, edge, 201)
n_arr = ThrDim.Gaussian(fit,z_arr)

plt.plot(z_arr, n_arr)
plt.title("Density along beam axis")
plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
plt.grid(); plt.show()

focal = Foc.Calc_Focus(n_arr, z_arr, gam)
thick = Foc.Calc_Square_Lens(n_0*1e17, focal, gam)
print(thick,"um equivalent thickness") 