#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 13:27:08 2020

@author: valentina_lee
"""

#%%
import sys
import os
sys.path.insert(0, "python")
import pyximport; pyximport.install()
import numpy as np
from beam.beams import laserpulse
from beam.beams import laserbeam
from beam.elements import plasma
from beam import interactions
from ionization import ionization
from lens import profile
from lens import bessel
from lens import ray
import matplotlib.pyplot as plt
from unwrap import unwrap

basePath = 'Lee/'
#%%
path = basePath + 'LaserRefraction/'
#%%
# Build the plasma
plasmaParams = {
    'Nx' : 2**12,
    'Ny' : 2**10,
    'Nz' : 2**8,
    'X' : 120e3,
    'Y' : 60e3,
    'Z' : 2.0e6, 
    'n0' : 1.0, #Nominal gas number density in 10^17 cm^-3.
    'atom' : ionization.Ar,
    'path' : path,
    'name' : 'ArgonPlasma',
    'load' : False,
    'cyl' : True
}
#%%
w = 200
z0 = 0.75e6 #The distance at which the uniform fully ionized plasma starts.
zf = plasmaParams['Z']
dz = 0.5e6 #The length of the fully ionized plasma.
sigmaIn = 15e4
sigmaOut = 15e4
z, ne = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
#%%
argon = plasma.Plasma(plasmaParams)
#%%
n = np.full((argon.Nx, argon.Ny, argon.Nz), plasmaParams['n0'], dtype='double')
r2 = argon.x[:, None, None]**2 + argon.y[None, :, None]**2
#%%
ne = plasmaParams['n0'] * ne[None, None, :] * np.exp(-r2/w**2)
argon.initialize_plasma(n, ne)
del n
del ne
#argon.plot_long_density_center()
#%%
# Build the beam
beamParams = {
    'Nx' : plasmaParams['Nx'],
    'Ny' : plasmaParams['Ny'],
    'X' : plasmaParams['X'],
    'Y' : plasmaParams['Y'],
    'lam' : 0.8, #wavelength
    'path' : path,
    'name' : 'ProbePulse',
    'load' : False,
    'threads' : 4,
    'cyl' : False,
    'E0' : 5.473494, # Gives peak input intensity of 3.98e10, that is, E=0.2mJ, beam dia= 8mm, pulse width= 100ps
    'waist' : 8000, #The spot size of the flattop region.
    #'z0' : 0.0,
    'order' : 4,
    'theta' : -2, #The angle of propagation in degrees, positive is angled upward.
    'dx' : -36e3
}

beam = laserbeam.GeneralSuperGaussianLaser(beamParams)
#%%
beam.plot_current_intensity()
#%%
interactions.beam_plasma(beam, argon)
#%%
beam.plot_intensity_at(50)
#%%
e = np.zeros((argon.Nz, beam.Nx), dtype='complex128')
Z = argon.Z
X = plasmaParams['X']
#%%
for i in range(argon.Nz):
    e[i, :] = beam.load_field(i)[0][:, int(beam.Ny/2)]
#%%
I = beam.intensity_from_field(e)
plt.figure(figsize=(14, 5))
plt.imshow(beam.prep_data(I), aspect='auto', extent=[0, Z, -X/2, X/2])
cb = plt.colorbar()
cb.set_label(r'Intensity')
plt.set_cmap('viridis')
plt.title('Longitudinal beam intensity')
plt.xlabel(r'z')
plt.ylabel(r'x')
plt.show()
#%%
e, z = beam.load_field(plasmaParams['Nz']-1)
I = beam.intensity_from_field(e)
phi = unwrap(np.angle(e)) - beam.k*np.radians(beam.theta)*beam.x[:, None]
phi -= np.amin(phi)
X = plasmaParams['X']
Y = plasmaParams['Y']
xlim = [0, 30e3]
ylim = [-15e3, 15e3]

plt.figure(figsize=(14, 6))
plt.subplot(121)
plt.imshow(beam.prep_data(I), aspect='auto', extent=[-X/2, X/2, -Y/2, Y/2])
cb = plt.colorbar()
cb.set_label(r'Intensity')
plt.set_cmap('viridis')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.xlim(xlim)
plt.ylim(ylim)

plt.subplot(122)
plt.imshow(beam.prep_data(phi), aspect='auto', extent=[-X/2, X/2, -Y/2, Y/2])
cb = plt.colorbar()
cb.set_label(r'Phase (rad)')
plt.set_cmap('cool')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.xlim(xlim)
plt.ylim(ylim)
plt.clim(3,6)
plt.show()

