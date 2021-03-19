#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:08:16 2020

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

gasden_arr = np.array([1e16,2e16,4e16,6e16,8e16,1e17])
hwhmz0_arr = np.array([34.7,34.9,34.9,34.3,33.8,33.6])
hwhmzc_arr = np.array([34.7,34.9,35.7,36.2,36.4,36.5])
zcent_arr  = np.array([-0.7,-1.5,-2.3,-3.1,-4.7,-5.5])*-1
nmaxz0_arr = np.array([1e16,2e16,4e16,5.8e16,6.9e16,7.4e16])
nmaxzc_arr = np.array([1e16,2e16,4e16,5.95e16,7.7e16,9.3e16])
charge_arr = np.array([1.213e-7,2.314e-7,4.208e-7,5.755e-7,7.046e-7,8.155e-7])
idealc_arr = np.array([1.213e-7,2.426e-7,4.852e-7,7.278e-7,9.704e-7,1.213e-6])

plt.title("Peak Plasma Density Along Beam Axis")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("Peak Plasma Density " + r'$(cm^{-3})$')
plt.plot(gasden_arr,gasden_arr, c='black', ls = 'dotted', label = "Ideal")
plt.plot(gasden_arr,nmaxzc_arr, c='blue', ls = 'solid', label = r'$x_{off}$')
plt.plot(gasden_arr,nmaxz0_arr, c='red', ls = '--', label = r'$x_{off}=0$')
plt.legend(); plt.grid(); plt.show()

plt.title("Plasma HWHM Along Beam Axis")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("HWHM " + r'$(\mu m)$')
plt.plot(gasden_arr,np.array([34.5,34.5,34.5,34.5,34.5,34.5]), c='black', ls = 'dotted', label = "Ideal")
plt.plot(gasden_arr,hwhmzc_arr, c='blue', ls = 'solid', label = r'$x_{off}$')
plt.plot(gasden_arr,hwhmz0_arr, c='red', ls = '--', label = r'$x_{off}=0$')
plt.legend(); plt.grid(); plt.show()

plt.title("Horizontal Offset towards Peak Density")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("Upstream Horizontal Offset " + r'$(mm)$')
plt.plot(gasden_arr,zcent_arr, c='blue', ls = 'solid', label = r'$x_{off}$')
plt.legend(); plt.grid(); plt.show()

plt.title("Total Ionized Charge vs Uniform Gas Density")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("Total Charge " + r'$(\mu C)$')
plt.plot(gasden_arr,charge_arr*1e6,c='blue', ls = 'solid',label="Simulation")
plt.plot(gasden_arr,idealc_arr*1e6,c='black', ls = 'dotted',label="Ideal Guess")
plt.legend(); plt.grid(); plt.show()
