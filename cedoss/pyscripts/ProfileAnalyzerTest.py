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

#Tanh
a = 391.484902775
b = 121.168309664
n_0 = 1.34929464489e+19
p_tanh = [a,b,n_0]
"""
#Gauss
n_0 = 1.46608792576e+19
sig = 314.194613916
x_0 = -14.4412644071
p_gauss = [n_0,sig,x_0]
"""
"""
#Tanh Tester
a = 9
b = 2
n_0 = 1
p_tanh = [a,b,n_0]
"""
p1 = p_tanh
#p1 = p_gauss

#Build our 3D approximate density
xwindow = 1500; xstep = .5 #microns
x = np.arange(-xwindow, xwindow, xstep)
arr = ThrDim.DoubleTanh(p1,x)

fwhm = Prof.CalcFWHM(x,arr)
flattop = Prof.CalcTopWidth(x,arr,.95)

plt.plot(x,arr, label = 'Tanh Profile')
plt.title("Distribution")
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