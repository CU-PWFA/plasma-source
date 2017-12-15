#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:13:31 2017

@author: chris
"""
import matplotlib.pyplot as plt
import numpy as np

den = [0.1,    0.78,  1.04,  1.57,  1.96,  3.13,  5.22, 7.83, 15.67]
fcl = [480.19, 61.71, 45.68, 29.67, 23.31, 13.92, 7.95, 5.41, 3.14]

den_np = np.linspace(.1,16)
fcl_np = 1129.6*1/den_np/24

plt.loglog(den,fcl,label='Simulations')
plt.plot(den_np,fcl_np,label='Analytic')
plt.axis([0.1,16,3,500])
plt.title("Logscale Focal Length vs Central Density from Gas Jet")
plt.xlabel("Gas Density (10^17 cm^-3)")
plt.ylabel("Integrated Focal Length (cm)")
plt.axvspan(3,16,facecolor = 'r', alpha = 0.3)
plt.text(3.1,110,"Significant Refraction")
plt.scatter(den,fcl)
plt.grid(); plt.legend(); plt.show()