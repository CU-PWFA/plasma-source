#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:03:22 2018

Quick script to take a gaussian profile and return the focal length

@author: chris
"""

import sys
sys.path.insert(0, "../")
import numpy as np
import matplotlib.pyplot as plt

from modules import ThreeDimensionAnalysis as ThrDim
from modules import TPLFocalLength as Foc

gam = Foc.gam_def

sig = 19.95*4
n_0 = 1e18 / 1e17

fit = [n_0,sig,0]
edge = 5*sig
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