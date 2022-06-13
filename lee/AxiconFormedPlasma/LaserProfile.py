#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 13:42:32 2020

@author: valentinalee
"""


"""
This is a file to establish the initial laser pulse parameters for ionization simulation.
*All the units are in SI unit*
*only plotting unit are as labels for comprehensing*

The spatial profile fitting equation is: 28.447e18*exp(-(r^2/(2*(29.204*1e-6)^2))^1.1799)

The temporal profile fitting equation is: 28.447e18*exp((-1/2)*(t/12.8e-15)^2)
which is a gaussian of a FWHM=30fs

Both spatial and temporal are plotted in this document for visualization
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
I0= 28.447e18
#%%
r= np.linspace(-200e-6, 200e-6, 1000)
Ir= np.exp(-(r**2/(2*(29.204*1e-6)**2))**1.1799)
#%%
plt.plot(r*1e6, Ir*I0*1e-4)
plt.xlabel('r ($\mu$m)')
plt.ylabel('I (W/cm^2)')
#%%
t= np.linspace(-100e-15, 100e-15, 1000)
It= np.exp((-1/2)*(t/12.8e-15)**2)
#%%
plt.plot(t*1e15, It*I0*1e-4)
plt.xlabel('Time (fs)')
plt.ylabel('I (W/cm^2)')
#%%
IR, IT= np.meshgrid(Ir, It)
Irt= IR*IT
plt.pcolormesh(t*1e15, r*1e6, Irt*I0*1e-4)
plt.xlabel('Time (fs)')
plt.ylabel('r ($\mu$m)')
plt.colorbar(label='I (W/cm^2)')
#%%
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
def FWHM(X,Y):
    half_max = max(Y) / 2.
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    return X[right_idx] - X[left_idx] 
#%%
FWHM(t, It)