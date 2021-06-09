#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 23:08:18 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
from lens import phaselens
from scipy import optimize

#%%
lam= 800e-9
k0= 2*np.pi/lam
points= 7000
X= 70e-3
xgrid= np.linspace(-X/2, X/2, points)
NF= np.exp(-(xgrid**2/(2*(8e-3)**2))**1)

f1= 0.75
f2= 0.5


d= 1.276
Eout1= phaselens.phase_lens_1d(X, NF, f1, 0, d, lam)
Eout2= phaselens.phase_lens_1d(X, Eout1, f2, 0, 8, lam)
#8.8
#9.9
#%%
Etest= phaselens.phase_lens_1d(X, NF, f1, 0, 0.5, lam)
FFphase= np.angle(Etest)

#%%
plt.plot(xgrid, abs(Etest))

#%%
plt.plot(xgrid, np.unwrap(FFphase))
#%%
plt.figure(1)
plt.plot(xgrid, abs(Eout2))
#%%
FFphase= np.angle(Eout2)
#%%
plt.figure(1)
plt.plot(xgrid, np.unwrap(FFphase))
#%%
radVals= np.unwrap(FFphase)[int(len(FFphase)/2):len(FFphase)]
azm= np.linspace(0, 2*np.pi, radVals.size)
th, r= np.meshgrid(radVals, azm)

radius= np.linspace(0, 1, int(len(xgrid)/2))
angle= np.linspace(0, 2*np.pi, radius.size)
r_grid, a_grid = np.meshgrid(radius, angle)
#%%
def polar_to_cartesian(data):
    new = np.zeros_like(data)+0
#    new = np.zeros_like(data) * np.nan
    x = np.linspace(-1, 1, new.shape[1])
    y = np.linspace(-1, 1, new.shape[0])
    for i in range(new.shape[0]):
        for j in range(new.shape[1]):
            x0, y0 = x[j], y[i]
            r, a = np.sqrt(x0**2 + y0**2), np.arctan2(y0, x0)
            data_i = np.argmin(np.abs(a_grid[:, 0] - a))
            data_j = np.argmin(np.abs(r_grid[0, :] - r))
            val = data[data_i, data_j]

            if r <= 1:
                new[i, j] = val
            print(i, j)

    return new

#%%
new = polar_to_cartesian(th)
newgrid= np.linspace(-X/2, X/2, int(points/2))
#%%
plt.figure(1)
plt.pcolormesh(newgrid, newgrid, new)
plt.axis('scaled')
plt.colorbar()

#%%
y= abs(Eout2)
x= xgrid
fitfn= lambda p, x: p[0]*np.exp(-((x)**2/(p[1]**2))**p[2])

errfunc = lambda p, x, y: fitfn(p, x) - y
p0= [3, 0.02, 2]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)
#%%
xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.plot(x, y, '.', label= 'plasma data')
plt.plot(xNew, fitfn(p1, xNew), '-', label='fitting')
plt.legend()

#%%
np.save('laser_para_88', p1)
np.save('FFphase_88', new)
np.save('grid_88', newgrid)

#%%
plt.figure(2)
plt.pcolormesh(th)
plt.colorbar()

#%%
laserPara= np.load('laser_para_88.npy')
FFphase=np.load('FFphase_88.npy')
PhaseGrid=np.load('grid_88.npy')

xg= PhaseGrid*1e6
yg= PhaseGrid*1e6
XG, YG= np.meshgrid(xg, yg)
R= np.sqrt(XG**2+YG**2)
beam= laserPara[0]*np.exp(-((R)**2/((laserPara[1]*1e6)**2))**laserPara[2])

#%%
beamPhase= beam*np.exp(-1j*FFphase)
propE= laser.fourier_prop2(beamPhase, xg, yg, 0e6, 0.8)

#%%
plt.figure(0)
plt.pcolormesh(xg, yg, (abs(propE[0, :, :])**2))