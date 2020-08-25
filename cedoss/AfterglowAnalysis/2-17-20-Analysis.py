#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:45:28 2020

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

gasden_arr = np.array([1e16,2e17,5e17,7e17,1e18,2e18,3e18])
hwhmz0_arr = np.array([16.2,19.2,25.7,27.6,26.8,11.5,9.4])
hwhmzc_arr = np.array([32.9,34.9,38.1,40.0,42.8,48.6,49.6])
charge_arr = np.array([3.090e-8,6.932e-7,1.732e-6,2.218e-6,2.825e-6,4.956e-6,7.133e-6])
idealc_arr = np.array([3.09e-8,6.18e-7,1.545e-6,2.163e-6,3.09e-6,6.18e-6,9.27e-6])

plt.title("Plasma HWHM Along Beam Axis")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("HWHM " + r'$(\mu m)$')
#plt.plot(gasden_arr,np.array([34.5,34.5,34.5,34.5,34.5,34.5]), c='black', ls = 'dotted', label = "Ideal")
plt.plot(gasden_arr,hwhmzc_arr, c='blue', ls = 'solid', label = r'$x_{off}=-1.7$ mm')
plt.plot(gasden_arr,hwhmz0_arr, c='red', ls = '--', label = r'$x_{off}=0$')
plt.legend(); plt.grid(); plt.show()

plt.title("Total Ionized Charge vs Uniform Gas Density")
plt.xlabel("Gas Density " + r'$(cm^{-3})$')
plt.ylabel("Total Charge " + r'$(\mu C)$')
plt.plot(gasden_arr,charge_arr*1e6,c='blue', ls = 'solid',label="Simulation")
plt.plot(gasden_arr,idealc_arr*1e6,c='black', ls = 'dotted',label="Ideal Guess")
plt.legend(); plt.grid(); plt.show()