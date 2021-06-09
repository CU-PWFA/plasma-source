#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 19:16:03 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as const
#%%
dt= 200e-12
dx= dt*const.c
theta_i= np.deg2rad(20)
m= 1
dlam= 100e-9
lam0= 800e-9
lamb= lam0-dlam/2
lamr= lam0+dlam/2
sp= 1e-3/1200

theta_i_n= 1000
L_n= 2000
Diffmisdx= np.empty((theta_i_n, L_n))
theta_i_a= np.linspace(0, np.pi/2, theta_i_n)
L_a= np.linspace(0, 1, L_n)

#theta_i= np.deg2rad(20)
#L= 0.5
for t_i in range (len(theta_i_a)):
    for l_i in range (len(L_a)):
        
        theta_i= theta_i_a[t_i]
        L= L_a[l_i]
        theta_mb= np.arcsin(np.sin(theta_i)-m*lamb/sp)
        theta_mr= np.arcsin(np.sin(theta_i)-m*lamr/sp)
        
        #print('b', theta_mb)
        #print('r', theta_mr)
        
        diff1= L/np.cos(theta_mb)-L/np.cos(theta_mr)
        
        if theta_mb*theta_mb>0:
            #Same sign
            diff2= abs(L*np.tan(theta_mr)-L*np.tan(theta_mb))*np.cos(theta_i)
        else:
            diff2= 0-abs(L*np.tan(theta_mr)+L*np.tan(theta_mb))*np.cos(theta_i)
        
        totdiff= diff1+diff2
        #print(diff1, diff2)
        #print(totdiff)
        
        Diffmisdx[t_i, l_i]= abs(totdiff-dx)

#%%
plt.pcolormesh(L_a, theta_i_a*180/np.pi, Diffmisdx, cmap= 'nipy_spectral')
plt.colorbar()
plt.xlabel('L (m)')
plt.ylabel('Incident Angle (°)')
#%%
theta_i= np.deg2rad(17)
L= 0.72
theta_mb= np.arcsin(np.sin(theta_i)-m*lamb/sp)
theta_mr= np.arcsin(np.sin(theta_i)-m*lamr/sp)

print('b', np.rad2deg(theta_mb))
print('r', np.rad2deg(theta_mr))

diff1= L/np.cos(theta_mb)-L/np.cos(theta_mr)

if theta_mb*theta_mb>0:
    #Same sign
    diff2= abs(L*np.tan(theta_mr)-L*np.tan(theta_mb))*np.cos(theta_i)
else:
    diff2= 0-abs(L*np.tan(theta_mr)+L*np.tan(theta_mb))*np.cos(theta_i)

totdiff= diff1+diff2
print(diff1, diff2)
print(abs(totdiff-dx))

