#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:42:08 2019

I felt bad for not including sig varying with axial distance.  To rectify this
I simply measured three values of sigma in the domains of 5mm +/- 1.2mm and
will do a quad fit.

This version uses the Lorentzian fit for the Radial profiles

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

y_arr = [-1200, 0, 1200] #um
gam_arr = [3111.51164927, 4020.58922632, 5058.88702969] #um
n0p_arr = [3.0919794422e+19, 2.56260482473e+19, 2.2206935743e+19]
x_0_arr = [-0.189293914047, -0.200114558806, -0.211465405831]

pfit = np.polyfit(y_arr, gam_arr, 2)
print('gam(y)=',pfit[2],'+',pfit[1],'y','+',pfit[0],'y**2')

print("gfit: ",pfit)
plt.plot(y_arr, gam_arr)
plt.title("Gamma Array")
plt.show()

pfit = np.polyfit(y_arr, n0p_arr, 2)
print('n0(y)=',pfit[2],'+',pfit[1],'y','+',pfit[0],'y**2')

print("nfit: ",pfit)
plt.plot(y_arr, n0p_arr)
plt.title("N0 Array")
plt.show()