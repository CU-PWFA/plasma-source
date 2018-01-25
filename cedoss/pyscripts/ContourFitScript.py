#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:36:14 2018

Assuming we have a 2d array, bmag_image, this script attempts to fit the
data to a square root scaling.  For use with gamma vs npl0 in
'DensityEnergyContour.py' and keep npl0_arr and gamma_arr as well

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim
from beampropagation import PlasmaPropagation as PProp

#npl0_arr = np.linspace(3e16, 1e17, num = 20)
#gamma_arr = np.linspace(1e4, 8e4, num = 20)
gamma_mat = np.zeros(len(npl0_arr))
#bmag_image = np.zeros((len(npl0_arr),len(gamma_arr)))

x_label = r'$ npl0 $ [cm^-3]'
y_label = r'$ \gamma $'

minloc = PProp.PlotContour(bmag_image, npl0_arr, gamma_arr, x_label, y_label)

for i in range(len(npl0_arr)):
    minlocy = np.argmin(bmag_image[i,:])
    gamma_mat[i] = gamma_arr[minlocy]

X = np.tile(npl0_arr.reshape(-1,1),(1,len(gamma_arr)))
Y = np.tile(gamma_arr.T,(len(npl0_arr),1))

levels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
labels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
fig, axes = plt.subplots(1,1, sharey=True)
CS = plt.contour(X,Y,bmag_image,levels,cmap=plt.get_cmap('Vega20b'))
plt.clabel(CS,labels,fontsize=9,inline=1,fmt='%1.2f')
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$B_m$')
cbar.set_ticks(levels)
plt.ylabel(y_label)
plt.xlabel(x_label)
plt.plot(npl0_arr, gamma_mat)
plt.show()

x_arr = np.linspace(min(gamma_mat), max(gamma_mat), 200)

lfit = np.polyfit(gamma_mat, npl0_arr, 1)
llabel = '{0:.2e}'.format(lfit[0]) + r'$\gamma$'
llabel = llabel + " + " + '{0:.2e}'.format(lfit[1])

plt.plot(gamma_mat, npl0_arr, label='Contour Data')
plt.plot(x_arr, lfit[0] * x_arr + lfit[1], label = llabel)
plt.title("Linear Fit")
plt.ylabel(x_label)
plt.xlabel(y_label)
plt.grid(); plt.legend(); plt.show()

y_label2 = r'$ \sqrt{\gamma} $'

rfit = np.polyfit(np.sqrt(gamma_mat), npl0_arr, 1)
rlabel = '{0:.2e}'.format(rfit[0]) + r'$\sqrt{\gamma}$'
rlabel = rlabel + " + " + '{0:.2e}'.format(rfit[1])

plt.plot(np.sqrt(gamma_mat), npl0_arr, label='Contour Data')
plt.plot(np.sqrt(x_arr), rfit[0] * np.sqrt(x_arr) + rfit[1], label=rlabel)
plt.title("Square Root Fit")
plt.ylabel(x_label)
plt.xlabel(y_label2)
plt.grid(); plt.legend(); plt.show()

plt.plot(gamma_mat, npl0_arr, label='Contour Data')
plt.plot(x_arr, lfit[0] * x_arr + lfit[1], label = llabel)
plt.plot(x_arr, rfit[0] * np.sqrt(x_arr) + rfit[1], label=rlabel)
plt.title("Linear Fit and Square Root Fit")
plt.ylabel(x_label)
plt.xlabel(y_label)
plt.grid(); plt.legend(); plt.show()