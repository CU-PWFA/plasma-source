#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:42:41 2020

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2

#%%
pixel= 2e-5 #m
lambda0=800e-9;
c= 3e8;
f0=c/lambda0;
k0=2*np.pi/lambda0; 
w0= 2*np.pi*f0
ne= 5*10*22 #(m^-3)
e= 1.6e-19;
me=9.11e-31;
epsilon0= 8.854e-12;
wp=np.sqrt(ne*e**2/(epsilon0*me))
epsilon= 1-(wp**2/w0**2)
n= np.sqrt(epsilon)
#%%
#Unit step fn
def us(x):
    return (np.sign(x)+1)/2;

#rect fn
def rect(x, halfwidth):
    return us(x+halfwidth)-us(x-halfwidth);

#%%
x= np.linspace(-600, 600, 1800)
y= np.linspace(-600, 600, 1800)
X, Y= np.meshgrid(x, y)
R= np.sqrt(X**2+Y**2)

sigma= 100
x0=0
Probe= np.exp(-((R-x0)**2/(2*sigma**2))**4)
#%%
W= abs(np.exp(-((Y-x0)**2/(2*2**2))**2)-1)
#%%
ProbeI= Probe*W
#%%
kx=fftfreq(1800, pixel)*2*np.pi;
ky=fftfreq(1800, pixel)*2*np.pi;
kX, kY= np.meshgrid(kx, ky);
kz= np.sqrt(k0**2-kX**2-kY**2);
L=4
FFL= (ifft2(fft2(ProbeI)*np.exp(-1j*kz*L)));

#%%
Winst= np.zeros((1800, 1800))
PFnext= Probe
for l in range (0, 50):
    Winst= Winst+1
    Winst[890: 910, (650+l*10): (650+(l+1)*10)]= W[890: 910, (650+l*10): (650+(l+1)*10)]
    ProbeI= PFnext*Winst
    PFnext= (ifft2(fft2(ProbeI)*np.exp(-1j*kz*0.02)));
    print(l)
   

#%%
FFL= (ifft2(fft2(PFnext)*np.exp(-1j*kz*2.5)));
    
    
#%%
plt.figure(3)
plt.pcolormesh(abs(FFL))
plt.axis('scaled')