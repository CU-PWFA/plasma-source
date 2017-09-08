#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:28:08 2017

@author: chris
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim
from modules import ProfileAnalyzer as Prof

#Gauss
n_0 = 1.29763441697e+19
sig = 2099.27300216
x_0 = -7.06194282574
p_gauss = [n_0,sig,x_0]
"""
#Tanh Tester
n_0 = 1
sig = 300
x_0 = 0
p_gauss = [a,b,n_0]
"""
p1 = p_gauss

#Build our 3D approximate density
xwindow = 7500; xstep = .5 #microns
x = np.arange(-xwindow, xwindow, xstep)
arr = ThrDim.Gaussian(p1,x)

fwhm = Prof.CalcFWHM(x,arr)
flattop = Prof.CalcTopWidth(x,arr,.95)

plt.plot(x,arr, label = 'Gauss Profile')
plt.title("Radial Distribution")
plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
plt.axvspan(-fwhm/2, -flattop/2, facecolor = 'b', alpha = 0.3, label='FWHM: '+str(fwhm))
plt.axvspan(-flattop/2, flattop/2, facecolor = 'g', alpha = 0.3, label='95% Len: '+str(flattop))
plt.axvspan(flattop/2, fwhm/2, facecolor = 'b', alpha = 0.3)
plt.grid(); plt.legend(); plt.show()

print("Parameters for distribution: "+str(p1))
print("Maximum Value              : "+str(max(arr)))
print("Full Width Half Max        : "+str(fwhm))
print("Flat Top Length (>95%)     : "+str(flattop))
print("Up Ramp Length (10%-95%)   : "+str(Prof.CalcUpRamp(x,arr)))
#print("Down Ramp Length (10%-90%) : "+str(Prof.CalcDownRamp(x,arr)))