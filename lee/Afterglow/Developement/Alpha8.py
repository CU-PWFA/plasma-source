#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 23:43:49 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
Te= 16000
#Te= 116045
ne= 1e17
ne= np.linspace(1e14, 1e18, 1000)
alpha= 1.46e-8*ne*Te**(-9/2)

#%%
print(alpha, '$cm^-3/s$')
#%%
plt.title('Te= 10eV= 116045K')
plt.plot(ne, alpha)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('$n_e ({cm}^{-3}$)')
plt.ylabel('$\\alpha_8 ({cm}^{-3}/s)$')