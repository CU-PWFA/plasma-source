#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 11:55:57 2020

@author: valentinalee
"""

#%% Import
import numpy as np;
import matplotlib.pyplot as plt;

#%%
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
#%%
#Den= np.load('LineOutArray_vac.npy')
Den_NGFirst= np.load('LineOutArray.npy')
#Den_NGLast= np.load('LineOutArray_LastFrame.npy')
den= Den_NGFirst
#LineOut= np.load('TestLineout.npy')
x= np.linspace(-1500, 1500, 200)
#%%
def FWHM(X,Y):
    half_max = max(Y) / 2.
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    return X[right_idx] - X[left_idx] 
#%%
FWHMarray= np.zeros(31)
for n in range (0, 31):
    y= den[:, int(n)]
    FWHMarray[int(n)]= FWHM(x, y)
#%%
time= np.linspace(0, 100, 31)
#%%
#FWHM_Vac= FWHMarray
FWHM_DF= FWHMarray
#%%
plt.plot(time, FWHM_Vac, label= 'FWHM_Diffuse to Vac')  
plt.plot(time, FWHM_DF, label= 'FWHM_Diffuse to Neutral')  
#plt.plot(time, FWHM_NGLast,  label= 'FWHM_NGLast')
plt.xlabel('time (ns)')
plt.ylabel('FWHM ($\mu$m)')
plt.legend()