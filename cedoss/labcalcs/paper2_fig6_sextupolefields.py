#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 18:13:52 2021

Just plotting some Sextupole Fields

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

reso = 100
y=np.linspace(-1,1,reso)
x=np.linspace(-1,1,reso)
z_sex=np.zeros((reso,reso))
z_tpl=np.zeros((reso,reso))
for i in range(len(x)):
    for j in range(len(y)):
        z_sex[j,i]=y[j]**3-3*x[i]**2*y[j]
        z_tpl[j,i]=y[j]**3+x[i]**2*y[j]
        
fig, (ax0,ax1) = plt.subplots(nrows=2)
fig.set_size_inches(5,8.5)
        
cf = ax0.contourf(x,y,z_sex)
fig.colorbar(cf, ax=ax0)
ax0.set_title(r'$\psi \sim y^3-3x^2y$')
ax0.set_ylabel("y")
ax0.text(-.95,0.87,"(a)",color='white')

cf = ax1.contourf(x,y,z_tpl)
fig.colorbar(cf, ax=ax1)
ax1.set_title(r'$\psi \sim y^3+x^2y$')
plt.xlabel("x")
plt.ylabel("y")
ax1.text(-.95,0.87,"(b)",color='black')

plt.show()