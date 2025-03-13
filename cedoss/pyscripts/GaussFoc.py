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
from modules import TPLFocalLength as Foc

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def SuperGauss(x, *p):
    A, mu, sigma, n = p
    return A*np.exp(-np.power( (x-mu)**2/(2.*sigma**2),n ))

gam = Foc.gam_def
"""
n_0 = 3.82965364484e+16
sig = 1111.05229343
x_0 = -1.09058401518e-09

n_0=n_0/1e17
edge = 6*sig
gcoeff=[n_0, x_0, sig]

z_arr = np.linspace(-edge, edge, 201)
n_arr = Gauss(z_arr,*gcoeff)
"""

n_0 = 31.5e+16
n_0=n_0/1e17

scoeff=[  8.65825767e+02 ,  5.02935833e-02 ,  6.76584409e-01 ,  2.43666488e+00]
scoeff[0]=n_0
scoeff[2]=scoeff[2]*1e4
edge = 6*scoeff[2]

z_arr = np.linspace(-edge, edge, 201)
n_arr = SuperGauss(z_arr,*scoeff)

plt.plot(z_arr, n_arr)
plt.title("Density along beam axis")
plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
plt.grid(); plt.show()

focal = Foc.Calc_Focus(n_arr, z_arr, gam)
thick = Foc.Calc_Square_Lens(n_0*1e17, focal, gam)
print(thick,"um equivalent thickness") 