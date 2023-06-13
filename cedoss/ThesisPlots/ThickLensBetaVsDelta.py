#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:42:01 2023

Plot Beta vs delta for thick lenses

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

delta_arr = np.linspace(-0.11,0.11,100)

re = 2.8179e-13 #cm
n0 = 3e16 #cm-3
gam0 = 19569.5
gam_arr = gam0*(1+delta_arr)

k0 = 2*np.pi*re*n0/gam0
#k_arr = k0/(1+delta_arr)
k_arr = 2*np.pi*re*n0/gam_arr

#betastar = 0.05 *100#cm
betastar = 0.2 *100#cm
emit0 = 3e-6 *100#cm-rad
d = 0

beta0 = betastar**2 +d**2/betastar
alph0 = -d/betastar
gamm0 = 1/betastar

sqrtkl_arr = np.array([0.1,0.3,0.6,0.9])
colors = plt.cm.brg(np.linspace(0, 1, 3*len(sqrtkl_arr)))

for i in range(len(sqrtkl_arr)):
    sqrtkl = sqrtkl_arr[i]
    lp = sqrtkl/np.sqrt(k0)
    
    betaf_arr = 1/(k_arr*beta0*np.square(np.sin(np.sqrt(k_arr)*lp))
            + gamm0*np.square(np.cos(np.sqrt(k_arr)*lp))  
            + np.sqrt(k_arr)*np.sin(2*np.sqrt(k_arr)*lp))
    
    fit = np.polyfit(delta_arr,betaf_arr,1)
    
    
    plt.semilogy(delta_arr,betaf_arr*np.sqrt(k0),c=colors[3*i],label=r'$\sqrt{K}L=$'+str(sqrtkl_arr[i]))
    plt.semilogy(delta_arr,(delta_arr*fit[0]+fit[1])*np.sqrt(k0),ls='--',c=colors[3*i+1])
plt.xlabel(r'$\delta$')
plt.ylabel(r'$\tilde\beta^*_{L,thick}(\delta)$')
plt.legend();plt.show()