#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 17:28:47 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
from lens import phaselens

#%%
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
def FWHM(X,Y):
    half_max = max(Y) / 2.
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    return X[right_idx] - X[left_idx] 
#%%
lam= 800e-9
k0= 2*np.pi/lam

X= 160e-3
Y= 160e-3
xgrid= np.linspace(-X/2, X/2, 1000)
ygrid= np.linspace(-Y/2, Y/2, 1000)
Xgrid, Ygrid= np.meshgrid(xgrid, ygrid)
r= np.sqrt(Xgrid**2+Ygrid**2)
NF= np.exp(-r**2/(2*(8e-3)**2))

kx=fftfreq(NF.shape[1], X/(NF.shape[1]))*2*np.pi
ky=fftfreq(NF.shape[0], Y/(NF.shape[0]))*2*np.pi
kX, kY= np.meshgrid(kx, ky);
kz= np.sqrt(k0**2-kX**2-kY**2);

f1= 0.75
lens1= np.exp(1j*k0*r**2/f1/2)
f2= 0.5
lens2= np.exp(1j*k0*r**2/f2/2)

Zout1= 0.75

d= 1.276
#%%
lam= 800e-9
k0= 2*np.pi/lam

X= 160e-3
xgrid= np.linspace(-X/2, X/2, 10000)
NF= np.exp(-xgrid**2/(2*(8e-3)**2))

f1= 0.75
f2= 0.5

Zout1= 0.75

d= 1.276
Eout1= phaselens.phase_lens_1d(X, NF, f1, 0, d, lam)
Eout2= phaselens.phase_lens_1d(X, Eout1, f2, 0, 9.7, lam)

#%%
plt.plot(abs(Eout2)**2)
#%%
EFWHM= np.zeros((10776))
for dz in range (0, 1277, 1):
    Eout1= phaselens.phase_lens_1d(X, NF, f1, 0, dz/1000, lam)
    EFWHM[dz]= FWHM(xgrid, Eout1)

for dz2 in range(0, 9500, 1):
    Eout2= phaselens.phase_lens_1d(X, Eout1, f2, 0, dz2/1000, lam)
    EFWHM[dz+dz2]= FWHM(xgrid, Eout2)

#%%
distance= np.linspace(0, 1.276, 1277)
np.linspace(0, 9.5, 9499)
#%%
Earray= np.zeros((2500, 1640))
for dz in range (0, 1277, 2):
    Eout1= phaselens.phase_lens_1d(X, NF, f1, 0, dz/1000, lam)
    Earray[:, int(dz/2)]= abs(Eout1[int(Eout1.shape[0]*3/8):int(Eout1.shape[0]*5/8)])

for dz2 in range(0, 1001, 1):
    Eout2= phaselens.phase_lens_1d(X, Eout1, f2, 0, dz2/1e2, lam)
    Earray[:, int(dz/2+dz2+1)]= abs(Eout2[int(Eout2.shape[0]*3/8):int(Eout2.shape[0]*5/8)])
#%%
dist= np.linspace(0, 1.276, 639)
dist= np.append(dist, np.linspace(0, 10, 1001)+1.276)
#dist= np.append(dist, np.zeros(1001))
raxis= xgrid[int(Eout1.shape[0]*3/8):int(Eout1.shape[0]*5/8)]
#%%
#plt.pcolormesh(Earray)
plt.pcolormesh(dist, raxis, Earray, cmap='gist_stern')
plt.colorbar()
plt.xlabel('Distance (m)')
plt.ylabel('r (m)')
#%%
#%%
distFW= np.linspace(0, 1.276, 1276)
distFW= np.append(distFW, np.linspace(0, 9.5, 9500)+1.276)

plt.plot(distFW, EFWHM)
plt.xlabel('Distance (m)')
plt.ylabel('FWHM (m)')
#%%
plt.plot(xgrid, abs(NF)**2)
#plt.plot(xgrid, abs(Eout2)**2)