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
from beampropagation import PlasmaPropagation as PProp

betalabel = 10

path = '/home/chris/Desktop/DataLoads/ContourBetaWaistSigma_10GeV_HighRes_PostFix/'+str(betalabel)+'cm/'
print(path)
bmag_image = np.load(path+'bmagarr.npy')
sigma_arr = np.load(path+'sig.npy')
zbeta_arr = np.load(path+'beta.npy')

#npl0_arr = np.linspace(3e16, 1e17, num = 20)
#gamma_arr = np.linspace(1e4, 8e4, num = 20)
sigma_mat = np.zeros(len(sigma_arr))
bmag_mat = np.zeros(len(sigma_arr))
#bmag_image = np.zeros((len(npl0_arr),len(gamma_arr)))

x_label = r'$ z_{\beta^*} $ [m]'
y_label = r'$ \sigma_{hw} $ [m]'

minloc = PProp.PlotContour(bmag_image, zbeta_arr, sigma_arr, x_label, y_label)

print("dx: ",zbeta_arr[1]-zbeta_arr[0])
print("dy: ",sigma_arr[1]-sigma_arr[0])

for i in range(len(zbeta_arr)):
    minlocy = np.argmin(bmag_image[i,:])
    sigma_mat[i] = sigma_arr[minlocy]
    bmag_mat[i] = bmag_image[i,minlocy]

X = np.tile(zbeta_arr.reshape(-1,1),(1,len(sigma_arr)))
Y = np.tile(sigma_arr.T,(len(zbeta_arr),1))
titleadd = '; '+r'$\beta^* = $'+str(betalabel)+' cm'

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
plt.title("Minimum Valley Contour"+titleadd)
plt.plot(zbeta_arr, sigma_mat)
plt.show()

x_arr = np.linspace(min(zbeta_arr), max(zbeta_arr), 200)
lfit = np.polyfit(zbeta_arr[:-1], sigma_mat[:-1], 1)

llabel = '{0:.2e}'.format(lfit[0]) + r'$ z_{\beta^*} $'
llabel_simp = llabel
llabel = llabel + " + " + '{0:.2e}'.format(lfit[1])

plt.plot(zbeta_arr, sigma_mat, label='Contour Data')
plt.plot(x_arr, lfit[0] * x_arr + lfit[1], label = llabel)
plt.plot(x_arr, lfit[0] * x_arr, label = llabel_simp)
plt.title("Linear Fit"+titleadd)
plt.ylabel(y_label)
plt.xlabel(x_label)
plt.grid(); plt.legend(); plt.show()

plt.plot(zbeta_arr, bmag_mat)
plt.title("B-mag along matched valley"+titleadd)
plt.ylabel("B-mag")
plt.xlabel(x_label)
plt.grid(); plt.show();

thresh_arr = [1.01]
for thresh in thresh_arr:
    print("Threshold = "+str(thresh))
    for i in range(1,len(bmag_mat)-1):
        if (bmag_mat[i] < thresh) & (bmag_mat[i-1] > thresh):
            zmin = zbeta_arr[i]
            print(" Minimum z_beta for B-mag < "+str(thresh)+": ",zmin," [m]")
            sigmax = lfit[0]*zmin + lfit[1]
            print(" Corresponding maximum sigma_hw: ",sigmax, "[m]")
        if (bmag_mat[i] < thresh) & (bmag_mat[i+1] > thresh):
            zmax = zbeta_arr[i]
            print(" Maximum z_beta for B-mag < "+str(thresh)+": ",zmax," [m]")
            sigmin = lfit[0]*zmax + lfit[1]
            print(" Corresponding minimum sigma_hw: ",sigmin, "[m]")