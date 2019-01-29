#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:18:26 2018

Given some 3D Gaussian distribution, what % of particles are in  the beam for
various cutoff weights?

@author: chris
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def Weights(z,y,x,sigx,sigy,sigz,cutoff):
    xdist = np.exp(-np.square(x)/(2*np.square(sigx)))
    ydist = np.exp(-np.square(y)/(2*np.square(sigy)))
    zdist = np.exp(-np.square(z)/(2*np.square(sigz)))
    weight = xdist * ydist * zdist
    if weight < cutoff:
        weight = 0
    return weight

sigx = 3.9e-6
sigy = sigx
sigz = 5.2e-6

mode = 0 #0 for single, 1 for array, 2 for just plotting 

dist = 7
ans0 = 1.2456679616112006e-15

if mode == 0:
    cutoff = 0.0034
    
    f = lambda x,y,z:Weights(z,y,x,sigx,sigy,sigz,cutoff)
    
    ans = integrate.tplquad(f, -dist*sigx, dist*sigx, lambda x: -dist*sigy, 
                            lambda x: dist*sigy, lambda x,y: -dist*sigz, lambda x,y: dist*sigz)[0]
    print(ans)
    print("Cutoff "+str(cutoff) + ": "+str(ans/ans0*100)+"%")
    
if mode == 1:
    cutoff_arr = np.linspace(0.00, 0.01, 100)
    ans_arr = np.zeros(len(cutoff_arr))
    for i in range(len(cutoff_arr)):
        cutoff = cutoff_arr[i]
        f = lambda x,y,z:Weights(z,y,x,sigx,sigy,sigz,cutoff)
        
        ans = integrate.tplquad(f, -dist*sigx, dist*sigx, lambda x: -dist*sigy, 
                                lambda x: dist*sigy, lambda x,y: -dist*sigz, lambda x,y: dist*sigz)[0]
        ans_arr[i] = ans/ans0*100
    mode = 2

if mode == 2:
    plt.plot(cutoff_arr,ans_arr)
    plt.xlabel("Weight Cutoff")
    plt.ylabel("% of charge")
    plt.grid(); plt.show();