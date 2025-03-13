#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 15:58:56 2022

With the betatron fits I got by manually recording all the good ones, make some presentation plots

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

beta_0   = np.array([178.7,135,126.9,155.8,212.2,254.2,212.1])
beta_1   = np.array([115.4,120.9,150.2,183.8,178.3,125.3,116.6])
beta_6   = np.array([70.01,63.44,59.45,59.02,67.66,83.94,153.1])
beta_24  = np.array([88.22,75.74,66.85,63.94,60.36,63.03,64.37])
beta_57  = np.array([27.4,25.08,30.89,47.81,65.62,79.32,92.18])
beta_115 = np.array([39.77,47.52,59.25,56.53,54.59,59.92,69.79])

emit_0   = np.array([37.42,38.55,35.92,26.6,19.16,16.64,21.41])
emit_1   = np.array([92.81,49.99,31.41,23.54,25.21,39.91,45.35])
emit_6   = np.array([123.8,93.35,77.87,74.67,76.15,107.5,139.57])
emit_24  = np.array([84.94,59.09,56.02,60.14,77.5,113.9,178.4])
emit_57  = np.array([149.7,99.71,83.52,97.8,139.4,219.63,344.3])
emit_115 = np.array([193.2,142.2,110.6,109.8,129,178.6,232.7])

slice_arr = np.array([20,30,40,50,60,70,80])
pressure_arr = np.array([0,1,6,24,57.8,115.8])
density_arr = pressure_arr * 3.0e15
for i in range(len(density_arr)):
    density_arr[i] = float('%.2g' % density_arr[i])

colors = plt.cm.brg(np.linspace(0, 1, len(pressure_arr)))

plt.figure(figsize=(7,5))
plt.plot(slice_arr,beta_0,c=colors[0],label="No Plasma")
plt.scatter(slice_arr,beta_0,c=colors[0])
plt.plot(slice_arr,beta_1,c=colors[1],label=r'$\mathrm{n_p \approx }$'+str(density_arr[1])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_1,c=colors[1])
plt.plot(slice_arr,beta_6,c=colors[2],label=r'$\mathrm{n_p \approx }$'+str(density_arr[2])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_6,c=colors[2])
plt.plot(slice_arr,beta_24,c=colors[3],label=r'$\mathrm{n_p \approx }$'+str(density_arr[3])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_24,c=colors[3])
plt.plot(slice_arr,beta_57,c=colors[4],label=r'$\mathrm{n_p \approx }$'+str(density_arr[4])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_57,c=colors[4])
plt.plot(slice_arr,beta_115,c=colors[5],label=r'$\mathrm{n_p \approx }$'+str(density_arr[5])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,beta_115,c=colors[5])
plt.xlabel("Projection Slice as Percent of Charge (%)")
plt.ylabel(r'$\beta^*$'+" from Betafunction Fit (cm)")
plt.legend();plt.show()

plt.figure(figsize=(7,5))
plt.plot(slice_arr,emit_0,c=colors[0],label="No Plasma")
plt.scatter(slice_arr,emit_0,c=colors[0])
plt.plot(slice_arr,emit_1,c=colors[1],label=r'$\mathrm{n_p \approx }$'+str(density_arr[1])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_1,c=colors[1])
plt.plot(slice_arr,emit_6,c=colors[2],label=r'$\mathrm{n_p \approx }$'+str(density_arr[2])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_6,c=colors[2])
plt.plot(slice_arr,emit_24,c=colors[3],label=r'$\mathrm{n_p \approx }$'+str(density_arr[3])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_24,c=colors[3])
plt.plot(slice_arr,emit_57,c=colors[4],label=r'$\mathrm{n_p \approx }$'+str(density_arr[4])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_57,c=colors[4])
plt.plot(slice_arr,emit_115,c=colors[5],label=r'$\mathrm{n_p \approx }$'+str(density_arr[5])+r'$\ cm^{-3}$')
plt.scatter(slice_arr,emit_115,c=colors[5])
plt.xlabel("Projection Slice as Percent of Charge (%)")
plt.ylabel(r'$\epsilon_N$'+" from Betafunction Fit (mm-mrad)")
plt.legend();plt.show()






