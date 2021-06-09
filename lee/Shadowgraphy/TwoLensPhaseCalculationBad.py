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
from beam.beams import laserbeam
from beam import interactions
from scipy import optimize
#%%
basePath = 'lee/'
path = basePath + 'LaserRefraction/'
beamParams = {
    'Nx' : 2**12,
    'Ny' : 2**12,
    'X' : 0.04e6,
    'Y' : 0.04e6,
    'lam' : 0.8, #wavelength
    'path' : path,
    'name' : 'ProbePulse',
    'load' : False,
    'threads' : 8,
    'cyl' : False,
    'E0' : 5.473494, # Gives peak input intensity of 3.98e10, that is, E=0.2mJ, beam dia= 8mm, pulse width= 100ps
    'waist' :  5e3,#The spot size of the flattop region.
    'z0' : 0.0,
    'order' : 3,
    'theta' : 0, #The angle of propagation in degrees, positive is angled upward.
    'dx' : 0
}

beam = laserbeam.GeneralGaussianLaser(beamParams)
#beam = laserbeam.GeneralSuperGaussianLaser(beamParams)

#%%lens phase
xgrid= np.linspace(-beamParams['X']/2, beamParams['X']/2, beamParams['Nx'], dtype='double')
ygrid= np.linspace(-beamParams['Y']/2, beamParams['Y']/2, beamParams['Ny'], dtype='double')
Xgrid, Ygrid= np.meshgrid(xgrid, ygrid)
rho= np.sqrt((Xgrid-beamParams['dx'])**2+Ygrid**2)
f1= 0.75e6
f2= 0.5e6
k0= 2*np.pi/(beam.lam)
phase1= np.exp(-1j*k0*rho**2/f1/2)
phase2= np.exp(-1j*k0*rho**2/f2/2)

interactions.beam_phase(beam, phase1)
#%%
beam.propagate(1.276e6, 1) #distance between two lenses --put in experimental distance

#%%
plt.figure(1)
#plt.pcolormesh(np.real(phase1))
#plt.colorbar()
plt.pcolormesh(abs(beam.e)**2)


#%%
interactions.beam_phase(beam, phase2)
beam.propagate(7.25e6, 1) #distance from lens 2 to 'before plasma interaction code' --put in experimental distance

#%%
beam.e

#%%
FFphase= np.angle(beam.e)

#%%
plt.pcolormesh((FFphase))

#%%
plt.figure(2)
plt.plot(y)
#%%
y= FFphase[:, int(beamParams['Nx']/2)]
#%%
x= xgrid
fitfn= lambda p, x: p[0]*np.exp(-((x)**2/(2*p[1]**2))**p[2])

errfunc = lambda p, x, y: fitfn(p, x) - y
p0= [5e16, 5e3, 2]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)
#%%
xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.plot(x, y, '.', label= 'plasma data')
plt.plot(xNew, fitfn(p1, xNew), '-', label='fitting')
plt.legend()

#%%
np.save('FFphase_para', p1)

#%%


