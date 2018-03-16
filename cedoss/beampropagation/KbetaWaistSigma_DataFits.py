#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:26:53 2018

@author: chris
"""

import matplotlib.pyplot as plt
import numpy as np

#for beta=10cm, no energy gain
kbeta_arr = [100., 150., 200., 250., 300.]
waist_arr= [-0.291, -0.342, -0.383, -0.408, -0.423]
sigma_arr= [0.123, 0.132, 0.140, 0.143, 0.145]
"""
#for beta=10cm, with energy gain
kbeta_arr = [100., 150., 200., 250., 300.]
waist_arr= [-0.296, -0.342, -0.383, -0.408, -0.423]
sigma_arr= [0.125, 0.132, 0.140, 0.143, 0.145]
"""
kbeta_arr_inv = 1/np.array(kbeta_arr)

#First is sigma vs beta
x_arr=np.linspace(kbeta_arr[0], kbeta_arr[-1], 200)
lfit1 = np.polyfit(kbeta_arr_inv, sigma_arr, 1)
#compfit1  = [-5.2272, 0.165]
#lfit1=[1.71,-0.0192]

llabel1 = '{0:.2e}'.format(lfit1[0]) + r'$ /k_{\beta^*} $'
llabel1 = llabel1 + " + " + '{0:.2e}'.format(lfit1[1])

plt.scatter(kbeta_arr, sigma_arr, label="Contour Data")
plt.plot(x_arr, lfit1[0] * 1/x_arr + lfit1[1], label=llabel1)
#plt.plot(x_arr, compfit1[0] * 1/x_arr + compfit1[1], label="kb scaling model")
plt.plot(x_arr, 1.546e-5*x_arr + 0.1537 - 3.608/x_arr, label="kb scaling model")
plt.title("Matched ramp half-length vs k_b")
plt.ylabel(r'$\sigma_{hw,m}$ [m]')
plt.xlabel(r'$k_{\beta^*}\,\mathrm{[m^{-1}]}$')
plt.grid(); plt.legend(); plt.show()

#################################################################

#Second is sigma vs beta
lfit2 = np.polyfit(kbeta_arr_inv, waist_arr, 1)
#lfit1=[1.71,-0.0192]

llabel2 = '{0:.2e}'.format(lfit2[0]) + r'$ /k_{\beta^*} $'
llabel_simp2 = llabel2
llabel2 = llabel2 + " + " + '{0:.2e}'.format(lfit2[1])

plt.scatter(kbeta_arr, waist_arr, label="Contour Data")
plt.plot(x_arr, lfit2[0] * 1/x_arr + lfit2[1], label=llabel2)
plt.plot(x_arr, -1.567e-4*x_arr - 0.447 + 19.06/x_arr, label="kb scaling model")
plt.title("Matched beam waist vs k_b")
plt.ylabel(r'$z_{\beta^*,m}$ [m]')
plt.xlabel(r'$k_{\beta^*}\,\mathrm{[m^{-1}]}$')
plt.grid(); plt.legend(); plt.show()