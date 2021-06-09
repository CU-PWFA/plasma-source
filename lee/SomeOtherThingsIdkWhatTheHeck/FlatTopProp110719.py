#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:50:55 2019

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image
from scipy import optimize

#%%
im01 = Image.open('17583372_1911050044_0007.tiff')
FP= np.array(im01)
#%%
pixel=3.45e-6
xaxis=np.linspace(-1023, 1024, 2048)*pixel*10**3
yaxis=np.linspace(-767, 768, 1536)*pixel*10**3

#%%
LineOut= FP[:,1024]

#%%
x=np.linspace(0, 1536, 1536)
y= -1/((np.exp(2*(((1536-x)/300)-1))+1))
plt.figure(2)
plt.plot(x, y)

#%%
x=np.linspace(0, 1536, 1536)
y= LineOut
plt.figure(3)
plt.plot(x, y, 'r.')
#%%
fitfn= lambda p, x: -p[0]/((np.exp(((1536-x)/p[2])**(4))+1))+p[3]
errfunc = lambda p, x, y: fitfn(p, x) - y;

p0= [25000, 1, 300, 2000];

p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y));

plt.figure(4);
plt.plot(x, y, ".", x, fitfn(p1, x), "r-");
plt.ylabel('Intensity (counts)')
plt.xlabel('X')
plt.title('order=4')
#%%
xlength=2048;
ylength=1536;
x=np.linspace(-xlength/2, xlength/2, xlength);
y=np.linspace(-ylength/2, ylength/2, ylength);
X, Y=np.meshgrid(x, y);
R=np.sqrt(X**2+Y**2);
FlatTop=fitfn(p0, R)
#%%
plt.figure(2)
plt.pcolormesh(xaxis, yaxis, FlatTop)
plt.axis('scaled')
#plt.colorbar(label='TW/cm^2')
#plt.xlabel('mm')
#plt.ylabel('mm')

#%%
plt.figure(2)
plt.pcolormesh(xaxis, yaxis, FP)
plt.axis('scaled')
plt.colorbar(label='TW/cm^2')
plt.xlabel('mm')
plt.ylabel('mm')
