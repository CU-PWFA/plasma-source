#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 11:45:56 2023

With the betatron fits I got by manually recording all the good ones, make some presentation plots

This version is using the set energy slices for my thesis

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#beta_0   = np.array([1.903, 1.83, 1.913, 1.737,1.318, 2.312, 1.82e-3])*100
beta_0   = np.array([1.903, 1.83, 1.913, 1.737,1.318, 2.312])*100
beta_1   = np.array([1.28, 1.17, .991, 1.43, 1.92, 2.51, 1.71])*100
#beta_6   = np.array([.934, .0077, .813, .484, .485, .933, 1.524])*100
beta_6   = np.array([.934, .813, .484, .485, .933, 1.524])*100
beta_24  = np.array([1.28, 1.15, .645, .429, .457, 1.21, 1.44])*100
#beta_57  = np.array([1.39e-4, 7.32e-4, .827, .473, .321, .805, .890])*100
beta_57  = np.array([.827, .473, .321, .805, .890])*100
#beta_115 = np.array([1.96, 6.09e-3, .942, .75, .608, .683, 1.04])*100
beta_115 = np.array([1.96, .942, .75, .608, .683, 1.04])*100

#emit_0   = np.array([2.21e-5, 2.24e-5, 1.94e-5, 2.42e-5, 3.94e-5, 3.77e-5, 1.80e-8])*1e6
emit_0   = np.array([2.21e-5, 2.24e-5, 1.94e-5, 2.42e-5, 3.94e-5, 3.77e-5])*1e6
emit_1   = np.array([5.12e-5, 4.30e-5, 4.38e-5, 3.00e-5, 2.43e-5, 3.65e-5, 1.17e-4])*1e6
#emit_6   = np.array([8.63e-5, 2.62e-8, 9.62e-5, 8.01e-5, 9.33e-5, 9.76e-5, 6.94e-5])*1e6
emit_6   = np.array([8.63e-5, 9.62e-5, 8.01e-5, 9.33e-5, 9.76e-5, 6.94e-5])*1e6
emit_24  = np.array([1.27e-4, 1.06e-4, 8.95e-5, 7.22e-5, 7.08e-5, 5.32e-5, 5.79e-5])*1e6
#emit_57  = np.array([2.03e-8, 2.38e-8, 1.36e-4, 1.04e-4, 9.04e-5, 1.07e-4, 1.57e-4])*1e6
emit_57  = np.array([1.36e-4, 1.04e-4, 9.04e-5, 1.07e-4, 1.57e-4])*1e6
#emit_115 = np.array([3.70e-5, 4.13e-8, 1.18e-4, 1.01e-4, 1.19e-4, 1.37e-4, 1.41e-4])*1e6
emit_115 = np.array([3.70e-5, 1.18e-4, 1.01e-4, 1.19e-4, 1.37e-4, 1.41e-4])*1e6

slice_arr = np.array([9.85, 9.9, 9.95, 10.0, 10.05, 10.1, 10.15])

slice_arr_p0 = np.array([9.85, 9.9, 9.95, 10.0, 10.05, 10.1])
slice_arr_p6 = np.array([9.85, 9.95, 10.0, 10.05, 10.1, 10.15])
slice_arr_p57 = np.array([9.95, 10.0, 10.05, 10.1, 10.15])
slice_arr_p115 = np.array([9.85, 9.95, 10.0, 10.05, 10.1, 10.15])

pressure_arr = np.array([0,1,6,24,57.8,115.8])
density_arr = pressure_arr * 2.7e15
for i in range(len(density_arr)):
    density_arr[i] = float('%.2g' % density_arr[i])

colors = plt.cm.brg(np.linspace(0, 1, len(pressure_arr)))

plt.figure(figsize=(7,5))
plt.plot(slice_arr_p0,beta_0,c=colors[0],label="No Plasma")
plt.scatter(slice_arr_p0,beta_0,c=colors[0])
plt.plot(slice_arr,beta_1,c=colors[1],label=r'$\mathrm{n_p \approx }$'+str(density_arr[1])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_1,c=colors[1])
plt.plot(slice_arr_p6,beta_6,c=colors[2],label=r'$\mathrm{n_p \approx }$'+str(density_arr[2])+r'$\ cm^{-3}$')
plt.scatter(slice_arr_p6,beta_6,c=colors[2])
plt.plot(slice_arr,beta_24,c=colors[3],label=r'$\mathrm{n_p \approx }$'+str(density_arr[3])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_24,c=colors[3])
plt.plot(slice_arr_p57,beta_57,c=colors[4],label=r'$\mathrm{n_p \approx }$'+str(density_arr[4])+r'$\ cm^{-3}$')
plt.scatter(slice_arr_p57,beta_57,c=colors[4])
plt.plot(slice_arr_p115,beta_115,c=colors[5],label=r'$\mathrm{n_p \approx }$'+str(density_arr[5])+r'$\ cm^{-3}$')
plt.scatter(slice_arr_p115,beta_115,c=colors[5])
plt.xlabel("Projection Slice (GeV)")
plt.ylabel(r'$\beta^*$'+" from Betafunction Fit (cm)")
plt.legend();plt.show()

plt.figure(figsize=(7,5))
plt.plot(slice_arr_p0,emit_0,c=colors[0],label="No Plasma")
plt.scatter(slice_arr_p0,emit_0,c=colors[0])
plt.plot(slice_arr,emit_1,c=colors[1],label=r'$\mathrm{n_p \approx }$'+str(density_arr[1])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_1,c=colors[1])
plt.plot(slice_arr_p6,emit_6,c=colors[2],label=r'$\mathrm{n_p \approx }$'+str(density_arr[2])+r'$\ cm^{-3}$')
plt.scatter(slice_arr_p6,emit_6,c=colors[2])
plt.plot(slice_arr,emit_24,c=colors[3],label=r'$\mathrm{n_p \approx }$'+str(density_arr[3])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_24,c=colors[3])
plt.plot(slice_arr_p57,emit_57,c=colors[4],label=r'$\mathrm{n_p \approx }$'+str(density_arr[4])+r'$\ cm^{-3}$')
plt.scatter(slice_arr_p57,emit_57,c=colors[4])
plt.plot(slice_arr_p115,emit_115,c=colors[5],label=r'$\mathrm{n_p \approx }$'+str(density_arr[5])+r'$\ cm^{-3}$')
plt.scatter(slice_arr_p115,emit_115,c=colors[5])
plt.xlabel("Projection Slice (GeV)")
plt.ylabel(r'$\epsilon_N$'+" from Betafunction Fit (mm-mrad)")
plt.legend();plt.show()






