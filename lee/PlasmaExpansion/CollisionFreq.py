#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 17:26:13 2020

@author: valentinalee
"""

#%% Import
import numpy as np;
import matplotlib.pyplot as plt;

#%%
n=1e23
me=9.1e-31
e=1.6e-19
epslon0=8.85e-12
kb=1.38e-23
Te=1*11604.45
eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epslon0)**2*np.sqrt(kb*Te)**3)*np.log(10)
nu=n*e**2*eta/me
nu1=nu*1e-9

print(nu1)
#%%
v=np.sqrt(kb*Te/me)
nu_0=n*e**4*8*np.pi/((4*np.pi*epslon0)**2*me**2*v**3)*np.log(10)
print(nu_0)