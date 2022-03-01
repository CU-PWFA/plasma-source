#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:30:22 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
import scipy.interpolate as spi
#%%
lam0=800e-9 
c=3e8/1e12
k0= 2*np.pi/lam0
w0= c/lam0*2*np.pi
t= np.linspace(-1000, 1000, 5000)
#t= np.linspace(-1000e-12, 1000e-12, 5000)
x= np.linspace(-80e-3, 80e-3, 500)
kx= fftfreq(len(x), np.amax(x)*2/len(x))*2*np.pi
ft= fftfreq(len(t), np.amax(t)*2/len(t))*2*np.pi
#chirpR= 3e-15/100
chirpR= -1e-10
ft= 3.5e-6*t+2.65e-15
#TW= 2
TW= 200
XW= 10e-3
Et= np.exp(-t**2/TW**2)
#Et= np.exp(-t**2/TW**2)*np.exp(-1j*chirpR*t**2)*np.exp(1j*w0*t)
#Et= np.exp(-t**2/TW**2)*np.exp(1j*w0*t)
Ef= np.exp(-x**2/XW**2)
#%%
EF, ET= np.meshgrid(Ef, Et)
E= ET*EF

#%%
FTE= (fft2(E))
#%%
iFTE= ifft2(fftshift(FTE))
#%%
plt.figure(2)
plt.title('Initial E Field')
plt.pcolormesh(np.real(E))
#plt.pcolormesh(np.real((iFTE)))
#plt.xlabel('x (m)')
#plt.ylabel('t (ps)')
plt.colorbar()
#%%
plt.title('Initial E Field in Fourier Plane')
plt.pcolormesh(fftshift(kx), ft, fftshift(np.real(FTE)))
#plt.pcolormesh(np.real((iFTE)))
plt.xlabel('k (1/m)')
plt.ylabel('f (1/ps)')
plt.colorbar()

#%%
def prop(FTEfield, dz):
    
    Ep= np.zeros((E.shape[0], E.shape[1]), complex)
    for slices in range (len(t)):
        k= (ft[slices])*2*np.pi/c
        kz= np.sqrt((k0+k)**2-kx**2)
        fftEs= FTE[slices, :]
        Eps=(fftEs*np.exp(-1j*kz*dz))
        Ep[slices, :]= Eps
    return Ep
#%%
def grating (d, FTE, theta_i=0):
    
    Eshr= np.zeros((E.shape[0], E.shape[1]), complex)
    SFTE= fftshift(FTE)
    for slices in range (len(t)):
        k= (ft[slices])*2*np.pi/c
        theta_d= np.arcsin(2*np.pi/k0/d-np.sin(theta_i))
        fftEs= SFTE[slices, :]
        knew= fftshift(kx)+1/np.cos(theta_d)*k*100
        f= spi.interp1d(knew, fftEs, fill_value="extrapolate")
        Eshr[slices, :]= f(fftshift(kx))
    return (fftshift(Eshr))

#%%
d= 1e-3/1200
l= 0.238
theta_i1= 0.45
#theta_i1= 0.01
#theta_i2= 1
theta_i2= -np.pi/2+1.4756
theta_d1= np.arcsin(lam0/d-np.sin(theta_i1))
theta_d2= np.arcsin(lam0/d-np.sin(theta_i2))
EFTgrd1= grating(d, FTE, theta_i1)
EFTprop= prop(EFTgrd1, l)
EFTgrd2= grating(d, EFTprop, theta_i2)
EFTprop2= prop(EFTgrd2, l)
EFTgrd3= grating(d, EFTprop2, theta_d2)
EFTprop3= prop(EFTgrd3, l)
EFTgrd4= grating(d, EFTprop3, theta_d1)

#%%
shear= grating(d, E)
#%%
#%%
plt.pcolormesh(np.real(ifft2(shear)))
#plt.pcolormesh(np.real(ifft2(EFTgrd4)))
#%%
Efinal= np.zeros(5000, complex)
Efinal[350:2850]= np.real(ifft2(EFTgrd4))[0:2500:, 250]/np.amax(np.real(ifft2(EFTgrd4))[0:2500, 250])
#%%
plt.plot(t, np.real(E)[:, 250], label= 'Initial')
plt.plot(t, Efinal, label= 'After Compressed')
plt.xlabel('t (ps)')
plt.ylabel('E (Normalized)')
plt.legend()
#%%
plt.pcolormesh(fftshift(kx), ft, np.real(EFTgrd))
plt.xlabel('k (1/m)')
plt.ylabel('f (1/ps)')
plt.colorbar()
#%%
Ep= prop(FTE, 0.1)
plt.pcolormesh(np.real(ifft2(Ep)))

#%%
prop(FTE)