#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 12:39:26 2019

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt

#%%
Delay= np.array([0, 30, 35, 40, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 10, 5, 15, 20, 25])
Energy= np.array([480, 388, 330, 216, 140.3, 120.7, 97.3, 69.8, 61.3, 46.4, 36.3, 28.1, 23.2, 18.6, 16.5, 469, 479, 463, 449, 424])

#%%
fitfn= lambda p, x: -p[0]/((np.exp(p[1]*(((1536-x)/p[2])-1))+1))+p[3]
errfunc = lambda p, x, y: fitfn(p, x) - y;

p0= [25000, 20, 300, 2000];

p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y));

#%%
plt.plot(Delay, Energy, '.')