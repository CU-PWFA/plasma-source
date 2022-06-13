#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 19:25:29 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
#%%
#ne= np.load('ne.npy')
n_lineout= np.load('Lineout.npy')
I_lineout= np.load('I.npy')
ext= np.load('Grid.npy')
#Int= np.load('I.npy')
#%%
x_grid= np.linspace(ext[2], ext[3], n_lineout.shape[0])
#%%
n_lineout_new= n_lineout[502:522]
x_grid_new= x_grid[502:522]
x_grid_new= x_grid_new-np.argmax(n_lineout_new)
#%%
x= x_grid_new
y= n_lineout_new/np.amax(n_lineout)
#%%
fitfn= lambda p, x: p[0]*np.exp(-((x-p[2])**2/(2*p[1]**2))**2)
#fitfn= lambda p, x: p[0]*np.exp(-1/2*((x-p[4])/p[1])**2)+p[2]*np.exp(-((x-p[4])**2/(2*p[3]**2))**2)

errfunc = lambda p, x, y: fitfn(p, x) - y;
p0= [5e16, 20, 9];
#p0= [1e17, 20, 1e16, 200, 9];
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)
#%%
xNew= np.linspace(np.amin(x_grid_new), np.amax(x_grid_new), 1000)
plt.plot(x, y, '.', label= 'plasma data')
plt.plot(xNew, fitfn(p1, xNew), '-', label='fitting')
plt.legend()

#%%
#plt.plot(x, y*1e6, '.', label= 'Split Step')
plt.plot(x*1e-6, y*1e6, '.', label= 'Split Step')
plt.legend()

#%%
#%%
#%%
I= np.zeros(49)
#I[46:53]= Int[99, 509:516]
I[22:29]= I_lineout[509:516]
#%%
I_All= I_lineout[500:525]
#%%
Igridsize= 2*ext[3]/I_lineout.shape[0]
#%%
x= np.linspace(-12*Igridsize, 12*Igridsize, 25)
y= I_All/np.amax(I_All)
#%%
fitfn= lambda p, x: p[0]*np.sinc(p[1]*x)+p[2]+p[3]*np.sinc(p[4]*x)
#fitfn= lambda p, x: p[0]*np.exp(-1/2*((x-p[4])/p[1])**2)+p[2]*np.exp(-((x-p[4])**2/(2*p[3]**2))**2)

errfunc = lambda p, x, y: fitfn(p, x) - y;
p0= [17, 0.04, 1, 1, 0.001];
#p0= [1e17, 20, 1e16, 200, 9];
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y), epsfcn= 2.5e-34);
print(p1)
#%%
xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.figure(2)
plt.plot(x, y, '.', label= 'plasma data_Intensity')
plt.plot(xNew, fitfn(p1, xNew), '-', label='fitting')
plt.legend()
#%%
#%%
x_grid= np.linspace(ext[2], ext[3], ne.shape[1])
z_grid= np.linspace(0, ext[1], ne.shape[0])
#%%
plt.plot(x_grid, n_lineout)
#%%
f= interp2d(x_grid, z_grid, ne, kind='cubic')
x_new= np.linspace(ext[2], ext[3], 50000)
z_new= np.linspace(0, ext[1], ne.shape[0]*2)
n= f(x_new, z_new)
#%%
n_lineout_new= n_lineout[502:522]
x_grid_new= x_grid[502:522]
xMax= np.amax(x_grid_new)

#%%
plt.pcolormesh(x_new, z_new, n)#%%
plt.plot(x_grid_new, n_lineout_new)
#%%
tck = interpolate.splrep(x_grid_new, n_lineout_new)
newX= np.linspace(-xMax, xMax, 2**10)

ne = interpolate.splev(newX, tck, ext=0) # Return 0 if out of range
#%%
plt.plot(newX, ne)
#%%
spl = UnivariateSpline(x_grid_new, n_lineout_new)
newX= np.linspace(-xMax, xMax, 2**10)
plt.plot(newX, spl(newX), 'g', lw=1)
#%%
f= interp2d(x_grid_new, n_lineout_new, kind='cubic')
newX= np.linspace(-xMax, xMax, 2**10)
n= f(newX)
