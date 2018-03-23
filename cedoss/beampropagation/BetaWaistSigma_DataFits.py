#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:26:53 2018

@author: chris
"""

import matplotlib.pyplot as plt
import numpy as np

"""#for k_b=210.37 
beta_arr = [0.01,0.025,0.05,0.10,0.20,0.30,0.50]
slop_arr = [-4.16e-1, -3.58e-1, -3.26e-1, -3.03e-1, -2.78e-1, -2.68e-1, -2.56e-1]
intr_arr = [0.41e-2, 0.81e-2, 1.35e-2, 2.30e-2, 4.29e-2, 5.89e-2, 9.16e-2]
waist_arr= [-0.0098, -0.051, -0.152, -0.393, -0.959, -1.53, -2.98]
sigma_arr= [0.00816, 0.0265, 0.0629, 0.142, 0.310, 0.471, 0.857]
"""
"""#for k_b=210.37, without the endpoints
beta_arr = [0.025,0.05,0.10,0.20,0.30]
slop_arr = [-3.58e-1, -3.26e-1, -3.03e-1, -2.78e-1, -2.68e-1]
intr_arr = [0.81e-2, 1.35e-2, 2.30e-2, 4.29e-2, 5.89e-2]
waist_arr= [-0.051, -0.152, -0.393, -0.959, -1.53]
sigma_arr= [0.0265, 0.0629, 0.142, 0.310, 0.471]
"""
#for k_b=212.5, PRL values
beta_arr = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
waist_arr= [-0.1490, -0.3816, -0.6469, -0.9449, -1.245, -1.543]
sigma_arr= [0.06204, 0.1384, 0.2194, 0.3051, 0.3898, 0.4724]
"""
#for k_b=297.51
beta_arr = [0.025,0.05,0.10,0.20,0.30]
sigma_arr= [0.0293,0.0652,0.145,0.321,0.495]
waist_arr= [-0.0648,-0.169,-0.424,-1.041,-1.683]
"""
"""
#for k_b = 94.08
beta_arr = [0.025,0.05,0.10,0.20,0.30]
sigma_arr= [0.02069,0.05241,0.1224,0.2759,0.4259]
waist_arr= [-0.02655,-0.09931,-0.2862,-0.7414,-1.2276]
"""
#First is sigma vs beta
x_arr=np.linspace(beta_arr[0], beta_arr[-1], 200)
lfit1 = np.polyfit(beta_arr, sigma_arr, 1)
pfit1 = np.polyfit(beta_arr, sigma_arr, 2)
#lfit1=[1.65,-0.0170]

llabel1 = '{0:.2e}'.format(lfit1[0]) + r'$ \beta^* $'
llabel_simp1 = llabel1
llabel1 = llabel1 + " + " + '{0:.2e}'.format(lfit1[1])

plabel1 = '{0:.2e}'.format(pfit1[0]) + r'$ (\beta^*)^2 $'
plabel1 = plabel1 + " + " + '{0:.2e}'.format(pfit1[1]) + r'$ \beta^* $'
plabel1 = plabel1 + " + " + '{0:.2e}'.format(pfit1[2])

plt.scatter(beta_arr, sigma_arr, label="Contour Data")
plt.plot(x_arr, lfit1[0] * x_arr + lfit1[1], label=llabel1)
plt.plot(x_arr, pfit1[0] * np.square(x_arr) + pfit1[1] * x_arr + pfit1[2], label=plabel1)
plt.title("Matched ramp half-length vs beta")
plt.ylabel(r'$\sigma_{hw,m}$ [m]')
plt.xlabel(r'$\beta^*$[m]')
plt.grid(); plt.legend(); plt.show()

pfit = np.polyfit(beta_arr, waist_arr, 2)
#pfit=[-3.92,-4.67,0.067]
llabel2 = '{0:.2e}'.format(pfit[0]) + r'$ (\beta^*)^2 $'
llabel_simp2 = llabel2
llabel2 = llabel2 + " + " + '{0:.2e}'.format(pfit[1]) + r'$ \beta^* $'
llabel_simp3 = llabel2
llabel2 = llabel2 + " + " + '{0:.2e}'.format(pfit[2])

plt.scatter(beta_arr, waist_arr, label="Contour Data")
#plt.plot(x_arr, pfit[0] * np.square(x_arr) + pfit[1] * x_arr + pfit[2], label=llabel2)
plt.plot(x_arr, pfit[0] * np.square(x_arr) + pfit[1] * x_arr + pfit[2], label="all 3 terms")
plt.plot(x_arr, pfit[1] * x_arr + pfit[2], label="no quad")
plt.plot(x_arr, pfit[0] * np.square(x_arr), label="just quad")
plt.plot(x_arr, pfit[0] * np.square(x_arr) + pfit[1] * x_arr, label="no offset")
plt.title("Matched waist position vs beta")
plt.ylabel(r'$z_{\beta*,m}$ [m]')
plt.xlabel(r'$\beta^*$[m]')
plt.grid(); plt.legend(); plt.show()