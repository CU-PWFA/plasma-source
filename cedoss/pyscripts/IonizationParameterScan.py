#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 09:17:07 2022

Script to plot a 2D contour of ionization thresholds for H2

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

#Set up constants and test that a single calculation works

chi = 15.4 # GeV, ionization energy of H2

"""
w0 = 12.096e-6 # um, spot size at focus
tau = 40e-15 # s, pulse width
El = 20e-3 #J, laser energy

P = El/tau
I = np.array([2*P/np.pi/np.power(w0,2)])/1e4
H = ThrDim.IonFracFromIntensity_1D(I,chi,tau)

print(I[0], "peak intensity")
print(H[0],"% ionization")
"""
#Set up parameter scan
#I'll run  this script at fixed values of tau, and find a 2d array of H over
# a scan of w0 and El.  The ThrDim fucntion is built for 1D arrays, so i can
# just do a single for loop

tau = 70e-15
w0_len = 300
El_len = 300
w0_array = np.linspace(10e-6,300e-6,w0_len) #spot sizes from 10um to 100um
El_array = np.linspace(0.5e-3, 20e-3, El_len) #laser energy from 0.5 mJ to 20 mJ
H_cont = np.zeros([len(w0_array),len(El_array)])


#Loop through the arrays and calculate ionziation

for i in range(len(w0_array)):
    P_array = El_array/tau
    w0 = w0_array[i]
    I_array = 2*P_array/np.pi/np.power(w0,2)/1e4
    H = ThrDim.IonFracFromIntensity_1D(I_array,chi,tau)
    H_cont[i]=H


#2D plot or something

plt.imshow(H_cont, interpolation="none", origin="lower",
    extent=[El_array[0]*1e3,El_array[-1]*1e3,w0_array[0]*1e6,w0_array[-1]*1e6],
    aspect='auto')
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$Ionization \ Fraction$')
plt.contour(El_array*1e3,w0_array*1e6,H_cont,levels = [0.99],colors='r')
plt.title("Ionization Rate with "+str(tau*1e15)+" fs Pulse Width")
plt.xlabel(r'$Laser \ Energy \ \mathrm{(mJ)}$')
plt.ylabel(r'$Spot \ Size \ \mathrm{(\mu m)}$')