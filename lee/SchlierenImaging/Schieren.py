#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 19:35:39 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
import scipy
import cmath
#%%
def phase_lens_2d(X, Y, Ein, F, Zin, Zout, lam):
    """ Generates the intensity profile that is an initial E field goes through a phase only lens at distance Z after the lens

    Parameters
    ----------
    X : double
        Physical width of the grid in the x direction, the grid goes from [-X/2, X/2].
    Y : double
        Physical width of the grid in the y direction, the grid goes from [-Y/2, Y/2].
    Ein : array-like, 2d [x][y]
        E field at the lens
    F : double
        focal lenth of the length
    Zin : double
        The distance from the input E field to the phase mask
    Zout : double
        Intensity at Zout where the intensity output will be. The distance from phase mask to output screen.
    lam : double
        The vacuum wavelength of the laser radiation. in micron

    Returns
    -------
    Eout : array-like, 2d
        Out put Efield 
    """
    k0= 2*np.pi/lam
    
    kx=fftfreq(2*Ein.shape[1], X/(2*Ein.shape[1]))*2*np.pi
    ky=fftfreq(2*Ein.shape[0], Y/(2*Ein.shape[0]))*2*np.pi
    kX, kY= np.meshgrid(kx, ky);
    kz= np.sqrt(k0**2-kX**2-kY**2);
    
    E= np.zeros((Ein.shape[0]*2, Ein.shape[1]*2), dtype=complex)
    E[int(E.shape[0]/4):int(E.shape[0]*3/4), int(E.shape[1]/4):int(E.shape[1]*3/4)]= Ein
    
    E0=ifft2(fft2(E)*np.exp(-1j*kz*Zin))

    xreal=np.linspace(-X, X, Ein.shape[1]*2)
    yreal=np.linspace(-Y, Y, Ein.shape[0]*2)
    Xreal, Yreal= np.meshgrid(xreal, yreal)                 
    phase_mask= np.exp(1j*k0*(Xreal**2+Yreal**2)/F/2/4)
    
    E_pm= E0*phase_mask
    Eout= ifft2(fft2(E_pm)*np.exp(-1j*kz*Zout))
    
    Eout= Eout[int(E.shape[0]/4):int(E.shape[0]*3/4), int(E.shape[1]/4):int(E.shape[1]*3/4)]
#    Iout= abs(Eout)**2
    return Eout

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
Den= 1e23
W= 50e-6
X= 10e-3
Z= 0.001
Y= 10e-3
xgrid= np.linspace(-X/2, X/2, 1000)
zgrid= np.linspace(-Z/2, Z/2, 500)
ygrid= np.linspace(-Y/2, Y/2, 1000)
Zg, Yg= np.meshgrid(zgrid, ygrid)
shift= 0
r= np.sqrt((Zg+shift)**2+Yg**2)
den= Den*np.exp(-(r**5)/(2*W**5))
c= 3e8
lam=800e-9
k0= 2*np.pi/lam
w= c/lam*2*np.pi
e= 1.6e-19
eps0= 8.85e-12
me= 9.1e-31
wp= np.sqrt(den*e**2/eps0/me)
npls= -(wp**2/(2*w**2))

#%%
#for zidx in range (50, 51):
for zidx in range (0, len(zgrid)):
    nx= np.empty(len(xgrid))
    nx.fill(1)
    NY, NX= np.meshgrid(npls[:, zidx], nx)
    np.save('n_z='+str(zidx), NY)
    print(zidx)
#%%
XG, YG= np.meshgrid(xgrid, ygrid)
rxy= np.sqrt(XG**2+YG**2)
E= np.exp(-rxy**2/(2*1.5e-3**2))
#%%
idx= 0
Phase= np.zeros(1)
kx=fftfreq(E.shape[1], X/(E.shape[1]))*2*np.pi
ky=fftfreq(E.shape[0], Y/(E.shape[0]))*2*np.pi
kX, kY= np.meshgrid(kx, ky);
kz= np.sqrt(k0**2-kX**2-kY**2);

#for idx in range(50, 51):
while zgrid[idx+1]< np.amax(zgrid):
    dz= zgrid[idx+1]-zgrid[idx]
    npls= np.load('n_z='+str(zidx)+'.npy')
    E= E*np.exp(1j*k0*dz*npls)
    E=ifft(fft(E)*np.exp(-1j*kz*dz))
    idx= idx+1
    print(idx)

#%%
f1= 0.1
E1= phase_lens_2d(X, Y, E, f1, 0.099, 0.1, lam)
#%%
block= np.zeros((len(xgrid), len(ygrid)))+1
block[0:int(len(ygrid)/2+8), :]= 0
#%%
E1b= E1*block
#%%
f2= 0.05
E2= phase_lens_2d(X, Y, E1b, f2, 0.05, 0.05, lam)

#%%
E2lineout= E2[:, int(E2.shape[1]/2)]
nplsL= npls[:, int(npls.shape[1]/2)]
#%%
E2L= E2[0:int(E2.shape[0]/2), int(E2.shape[1]/2)]

#%%
plt.pcolormesh(ygrid*1e3, xgrid*1e3, abs(E2)**2)
#plt.pcolormesh(ygrid, xgrid, abs(Etest)**2)
plt.axis('scaled')
#plt.title('E0')
plt.title('At Second Focus')
plt.xlabel('mm')
plt.ylabel('mm')
#%%
#plt.pcolormesh(block)
plt.pcolormesh(npls)
plt.colorbar()

#%%
#plt.plot(ygrid[0:int(E2.shape[0]/2)], abs(E2L)**2)
plt.plot(ygrid, abs(E2lineout)**2/np.amax(abs(E2lineout)**2), label= 'Intensity')
plt.plot(ygrid, -nplsL/np.amax(abs(nplsL)), label= 'Plasma')
plt.legend()