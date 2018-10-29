#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 10:49:30 2018

Simple plot of beam density vs waist beta

@author: chris
"""
import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
import numpy as np
import matplotlib.pyplot as plt

beta_i_arr = np.linspace(10, 20, 100) #cm
nmax_arr = np.zeros(len(beta_i_arr))

#emitx = 5.3e-6 *100 #cm-rad
#emitx = 3e-6*100 #cm-rad
emitx = 100e-9*100 #cm-rad
gam = Foc.gam_def
nbeam = 6e9
sigz = 5.2e-6*100 #cm
nmax_arr = nbeam/(2*np.pi)**(3/2)/sigz/(beta_i_arr*emitx/gam)

plt.plot(beta_i_arr,nmax_arr)
plt.title("Beam Density vs Initial Beta")
plt.ylabel(r'$n_0\mathrm{\ [cm^{-3}]}$')
plt.xlabel(r'$\beta_i$ [cm]')
plt.grid(); plt.show()