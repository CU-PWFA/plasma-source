#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:49:38 2019

Plotting up some of the results in analyzing the SFQED simulations for
single bunch focusing

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
from modules import TPLFocalLength as Foc

density = np.array([1.00E+14, 2.00E+14,4.00E+14,6.00E+14,8.00E+14,1.00E+15,2.00E+15,5.00E+15,1.00E+16])
inrloca = np.array([180.3,440,405,312,256,208,112,56.1,47.7])
inrperc = np.array([100,46,46,46,45,43,36,33,26])
inrsigr = np.array([36.28,27.9,19.31,14.56,11.71,9.89,5.87,3.96,3.89])
inrener = np.array([1.0400851396,0.8090081063,1.6888719872,2.9705613754,4.4926457129,6.0183760496,14.3030858841,28.8089225589,23.5221813231])

tpl_l = 737e-6 * 100#cm
beta_i = 10.0 *100 #cm
gamma = 25440
emit = 3.5e-6 * 100 #cm

density_theory = np.linspace(density[0],density[-1],1000)

inrloca_theory = np.zeros(len(density_theory))
inrsigr_theory = np.zeros(len(density_theory))
for i in range(len(density_theory)):
    tpl_k = Foc.Calc_K(density_theory[i], gamma)
    inrloca_theory[i] = Foc.Calc_ThickWaistPos_DeltaOff_UnNormalized(tpl_k, tpl_l, beta_i, 0)
    inrsigr_theory[i] = np.sqrt(Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(tpl_k, tpl_l, beta_i, 0) * emit / gamma)*1e4
    
plt.semilogx(density, inrloca)
plt.scatter(density, inrloca, label = "Simulations")
plt.semilogx(density_theory, inrloca_theory, label = "Ideal Theory")
plt.title("Location of Inner Waist")
plt.xlabel(r'$\mathrm{Density \ }[cm^{-3}]$')
plt.ylabel(r'$\mathrm{Waist \ Location \ }[cm]$')
plt.grid(); plt.legend(); plt.show();

plt.semilogx(density, inrperc)
plt.scatter(density, inrperc, label = "Simulations")
plt.title("Percentage of Beam in Inner Waist")
plt.xlabel(r'$\mathrm{Density \ }[cm^{-3}]$')
plt.ylabel(r'$\mathrm{Inner \ Percent \ }[\%]$')
plt.grid(); plt.show();

plt.semilogx(density, inrsigr)
plt.scatter(density, inrsigr, label = "Simulations")
plt.semilogx(density_theory, inrsigr_theory, label = "Ideal Theory")
plt.title("Size of Inner Waist")
plt.xlabel(r'$\mathrm{Density \ }[cm^{-3}]$')
plt.ylabel(r'$\mathrm{Waist \ Size \ }[\mu m]$')
plt.grid(); plt.legend(); plt.show();

plt.semilogx(density, inrener)
plt.scatter(density, inrener, label = "Simulations")
plt.title("Inner Waist Density Factor Increase")
plt.xlabel(r'$\mathrm{Density \ }[cm^{-3}]$')
plt.ylabel(r'$\mathrm{Factor \ Increase}$')
plt.grid(); plt.show();