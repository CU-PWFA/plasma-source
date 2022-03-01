#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:38:16 2018

@author: valentina
"""
#%%
import numpy as np;
import matplotlib.pyplot as plt;
import math;
import cmath;
from scipy.fftpack import fft2, ifft2, fftshift, ifftshift, fftfreq;
#%%
def FermiDirac (ToverTF, EoverEF):
    TFoverT= 1/ToverTF;
    muoverEF= 0.001-(math.pi**2/12)*((ToverTF)**2);
    kB= 1.38e-23;
    return 1/(np.exp(TFoverT*((EoverEF)-muoverEF))+1);
#%%
#field distribution
length= 0.01;
npoint=2000;
TTFratio=0.0001;
x= np.linspace(-length, length, npoint);
y= np.linspace(-length, length, npoint);
X, Y = np.meshgrid(x, y);
R= np.sqrt(X**2+Y**2);
I=FermiDirac (TTFratio, R);

It=FermiDirac (TTFratio, R);
plt.figure(1);
plt.plot(R, It);

#%%
#L is the screen location, f is the focal length
L=0.1;
Rr=0.1;

#some fixed parameters
lambda0= 780*math.pow(10, -9);
k0= 2*math.pi/ lambda0;
c=3*math.pow(10, 8);
f0= c/lambda0;
omega0=2*math.pi*f0;

#get e field from intensity
E= np.sqrt(I);

plt.figure(2);
plt.pcolormesh(X, Y, E);
plt.gca().set_aspect('equal');
plt.xlabel('x(m)');
plt.ylabel('y(m)');
plt.title('I @ z=0');

#parabolic phase change, (thin lens)
ps= np.exp(1j*k0*(X**2+Y**2)/Rr/2);
#E field after thin lens 
Em= E*ps;

#take a fft of E field
Ef= ((fft2(Em)));

kx= fftfreq(npoint, 2*length/npoint)*2*math.pi;
ky= kx;
kX, kY = np.meshgrid(kx, ky);

kxy=np.sqrt(k0**2-kX**2-kY**2);
H= np.exp(-1j*kxy*L);
G= Ef*H;
absEf=abs(Ef);
absG=abs(G);

plt.figure(3);
plt.pcolormesh(kx, ky, np.abs(Ef));
plt.gca().set_aspect('equal');
plt.xlabel('x(m)');
plt.ylabel('y(m)');
plt.title('Ef');

plt.figure(4);
plt.pcolormesh(kx, ky, np.abs(G));
plt.gca().set_aspect('equal');
plt.xlabel('x(m)');
plt.ylabel('y(m)');
plt.title('G');


#take inver fft get back to space domain
g= abs((ifft2((G))));

#everything makes sense.....
plt.figure(5);
plt.pcolormesh(X, Y, g);
plt.gca().set_aspect('equal');
plt.xlabel('x(m)');
plt.ylabel('y(m)');
plt.title('I @ focal waist');

#%%
