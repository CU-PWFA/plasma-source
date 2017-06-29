# -*- coding: utf-8 -*-
"""
Tests the GaussianBeam.py for creating gaussian beam pulses and analyzes
them in 2D

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../")

from modules import Doss_Ionization as adk
from modules import GaussianBeam

#Set up constants and size of our domain in z and r
size = 500
domain = np.arange(-size,size)
counter = domain+size
r=domain/size*1e-4
z=(domain/size)*2e-1
z0 = 0
I0 = 1

wavelength = 100e-9
k = np.pi*2/(wavelength)
w0 = np.sqrt(0.159*wavelength/2)*0.4
zR = np.power(w0,2)*np.pi/wavelength

#Print some diagnostic values
print("Waist:",w0)
print("Rayleigh:",GaussianBeam.Gauss_zR(w0,wavelength))
print("Wavelength:",wavelength)

#Gathers the itensity of our gaussian beam in both z and r
#Our domain is -size to size, so our 2D array is of 2*size
I=np.empty([size*2,size*2])
for i in counter:
    for j in counter:
        I[i][j]=GaussianBeam.GaussianBeamIntensity(I0,k,w0,z0,z[i],r[j])

#Takes our intensity and uses ADK to calculate ionization fraction
#Uses some values for ADK:  I_0=2.0e14, delt_t=100e-15, E_i=13.6, Z=1
H=np.empty([size*2,size*2])
for i in domain:
    H[i] = adk.TotalProb(adk.Elefield(I[i]*2.0e14),100e-15,13.6,1)

#Plot intensity
plt.imshow(I, interpolation="none", origin="lower",
           extent=[r[0],r[2*size-1],z[0],z[2*size-1]], aspect='.001')
CB=plt.colorbar()
CB.set_label('I/I0')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('Relative Intensity of a Gaussian Beam')
plt.show()

#Contour plot of intensity
plt.figure()
levels = np.arange(0.05,0.95,0.1)
CS=plt.contour(r*10000,z,I,levels)
CB=plt.colorbar(CS,extend='both')
CB.set_label('I/I0')
plt.ylabel('z (m)')
plt.xlabel('r (10^-4 m)')
plt.title('Relative Intensity of a Gaussian Beam')
plt.show()

#Plot ionization fraction
plt.imshow(H, interpolation="none", origin="lower",
           extent=[r[0],r[2*size-1],z[0],z[2*size-1]], aspect='.001')
plt.colorbar()
CB.set_label('n_i/n_0')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('Ionization Fraction for I_0 = 2.0e14 and E_ionization = 13.6')
plt.show()
