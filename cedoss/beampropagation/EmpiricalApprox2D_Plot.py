#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:46:21 2018

plots the contour over beta and sigma

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp
import matplotlib.pyplot as plt

betalabel = 10

path = '/home/chris/Desktop/DataLoads/ContourEmpiricalApprox1/'
print(path)
bmag_image = np.load(path+'bmagarr.npy')
kb_arr = np.load(path+'kb.npy')
beta_arr = np.load(path+'beta.npy')

minloc = PProp.PlotContour(bmag_image, kb_arr, beta_arr, r'$k_b \,\mathrm{[m^{-1}]}$', r'$\beta^*\,\mathrm{[m]}$')

count = 0
jlen = len(beta_arr)
nlen = len(kb_arr)*len(beta_arr)
norm_arr = np.zeros(nlen)
bmag_norm = np.zeros(nlen)
for i in range(len(kb_arr)):
    for j in range(len(beta_arr)):
        beta_n = kb_arr[i]*beta_arr[j]
        norm_arr[i*jlen+j] = beta_n
        bmag = bmag_image[i][j]
        if bmag < 1.01:
            count = count + 1
#        else:
#            bmag = 1.1
        bmag_norm[i*jlen+j] = bmag

plt.scatter(norm_arr, bmag_norm, s=1)
plt.plot([0,150],[1.01,1.01],ls='--' , color='orange')
plt.title('B-mag vs Normalized ' + r'$\beta^*$')
plt.ylabel("B-mag")
plt.xlabel(r'$\widetilde{\beta^*}$')
plt.ylim(0.98,1.1)
plt.xlim(0, 150)
plt.grid(); plt.show()
tol_arr = np.zeros(count)
count = 0
for k in range(len(bmag_norm)):
    if bmag_norm[k] < 1.01:
        tol_arr[count] = norm_arr[k]
        count=count+1
        
print("Tolerances for B-mag < 1.01:")
print(" Min: ", min(tol_arr))
print(" Max: ", max(tol_arr))