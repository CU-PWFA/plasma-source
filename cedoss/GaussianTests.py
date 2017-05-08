# -*- coding: utf-8 -*-
"""
@author: chris
"""

import numpy as np
import GaussianBeam
import matplotlib.pyplot as plt

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

print("Waist:",w0)
print("Rayleigh:",GaussianBeam.Gauss_zR(w0,wavelength))
print("Wavelength:",wavelength)

I=np.empty([size*2,size*2])
for i in counter:
    for j in counter:
        I[i][j]=GaussianBeam.GaussianBeamIntensity(I0,k,w0,z0,z[i],r[j])
"""
H=np.empty([size,size])
for i in domain:
    H[i] = adk.f(adk.Elefield(I[i]*1.4e14),13.6,1)
"""
plt.imshow(I, interpolation="none", origin="lower",
           extent=[r[0],r[2*size-1],z[0],z[2*size-1]], aspect='.001')
CB=plt.colorbar()
CB.set_label('I/I0')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('Relative Intensity of a Gaussian Beam')
plt.show()

plt.figure()
levels = np.arange(0.05,0.95,0.1)
CS=plt.contour(r*10000,z,I,levels)
CB=plt.colorbar(CS,extend='both')
CB.set_label('I/I0')
plt.ylabel('z (m)')
plt.xlabel('r (10^-4 m)')
plt.title('Relative Intensity of a Gaussian Beam')
plt.show()
"""
plt.imshow(H, interpolation="none", origin="lower",
           extent=[r[0],r[size-1],z[0],z[size-1]], aspect='.001')
plt.colorbar()
plt.show()
"""