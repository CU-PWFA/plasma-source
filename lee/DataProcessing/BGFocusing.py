#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:38:55 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from lens import phaselens
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
#%%
def us(x):
    return (np.sign(x)+1)/2

def rect(x, halfwidth):
    return us(x+halfwidth)-us(x-halfwidth)

#%%
#z0= 1
#z= 0
W0=1e-3
lam= 800e-9
#W0= np.sqrt(lam*z0/np.pi)
z0= W0**2*np.pi/lam
Wz= W0*np.sqrt(1+(z/z0)**2)
Rz= z*(1+(z0/z)**2)
eta= np.arctan(z/z0)
#phase= (2*np.pi/lam)*z-eta+((2*np.pi/lam)*r**2/(2*Rz))
#%%
f= 3.39
Wo= Wz
W0p= lam/np.pi/Wo*f
#%%
z= np.linspace(-20, 20, 1000)
Wz= W0*np.sqrt(1+(z/z0)**2)
#%%
plt.plot(z, Wz, 'b')
plt.plot(z, -Wz, 'b')
plt.xlabel('z (m)')
plt.ylabel('W(z) (m)')
plt.title('z0='+str(round(z0, 2)))

#%%
BG= np.sqrt(np.array(Image.open('2006260039/19423601_2006260039_0001.tiff')))
#%%
#%%Imaging
xpixel= 0.8e-5
ypixel= 0.8e-5
Physical_X= BG.shape[1]*xpixel
Physical_Y= BG.shape[0]*ypixel
xreal=np.linspace(-Physical_X/2, Physical_X/2, BG.shape[1])
yreal=np.linspace(-Physical_Y/2, Physical_Y/2, BG.shape[0])
imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, BG, 10, 0, 10, 800e-9)
#%%
plt.figure(2)
plt.pcolormesh(xreal, yreal, abs(imaging))
#%%
xr= np.zeros(BG.shape[1])+1
Xr, Yr= np.meshgrid(xr, yreal)
wire= (-rect(Yr, 100e-6)+1)*Xr
#%%
newIm= imaging*wire
#%%
k0= 2*np.pi/lam

kx=fftfreq(BG.shape[1], Physical_X/(BG.shape[1]))*2*np.pi
ky=fftfreq(BG.shape[0], Physical_Y/(BG.shape[0]))*2*np.pi
kX, kY= np.meshgrid(kx, ky);
kz= np.sqrt(k0**2-kX**2-kY**2);

prop= ifft2(fft2(newIm)*np.exp(-1j*kz*1))

#%%
imaging2=phaselens.phase_lens_2d(Physical_X, Physical_Y, newIm, 0.5, 1.8, 0.69, 800e-9)
