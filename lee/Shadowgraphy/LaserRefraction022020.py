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
from beam.beams import laserbeam
from beam.elements import plasma
from beam import interactions
from ionization import ionization
from lens import profile
from lens import phaselens
from lens import ray
import matplotlib.pyplot as plt
from unwrap import unwrap
from propagation import laser
from scipy import interpolate; 

'''
If laser plasma interaction fail, check the grid density. If the spacing k*theta is larger than the grid spacig, it will fail.
Try with a smaller crossing angle to confirm the code works. 
'''
#%%
basePath = 'lee/'
path = basePath + 'LaserRefraction/'

# Build the plasma
plasmaParams = {
    'Nx' : 2**12,
    'Ny' : 2**8,
    'Nz' : 2**8,
    'X' : 0.100e6,
    'Y' : 0.060e6,
    'Z' : 2.0e6, 
    'n0' : 1, #Nominal gas number density in 10^17 cm^-3.
    'atom' : ionization.He,
    'path' : path,
    'name' : 'HePlasma',
    'load' : False,
    'cyl' : True
}


w = 200
z0 = 0.95e6 #The distance at which the uniform fully ionized plasma starts.
zf = plasmaParams['Z']
dz = 0.1e6 #The length of the fully ionized plasma.
sigmaIn = 0.15e6
sigmaOut = 0.15e6
z, ne = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
He = plasma.Plasma(plasmaParams)

n = np.full((He.Nx, He.Ny, He.Nz), plasmaParams['n0'], dtype='double')
r2 = He.x[:, None, None]**2 + He.y[None, :, None]**2



ne = plasmaParams['n0'] * ne[None, None, :] * np.exp(-r2/w**2)
He.initialize_plasma(n, ne)
del n
del ne

He.plot_long_density_center()
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
    'threads' : 8,
    'cyl' : False,
    'E0' : 5.473494, # Gives peak input intensity of 3.98e10, that is, E=0.2mJ, beam dia= 8mm, pulse width= 100ps
    'waist' : 8000, #The spot size of the flattop region.
    'z0' : 0.0,
    'order' : 3,
    'theta' : 0.573,#0.859,#2.06175769, #The angle of propagation in degrees, positive is angled upward.
    'dx' : -0.01e6#-0.015e6
}

beam = laserbeam.GeneralSuperGaussianLaser(beamParams)
#%%
#beam.plot_current_field()
plt.figure(2)
beam.plot_current_intensity()
#%%
interactions.beam_plasma(beam, He)
#%%
plt.figure(4)
beam.plot_intensity_at(255)
#%%
#fix propagation angle
NewBeam = beam.e*np.exp(-1j*beam.k*np.radians(beam.theta)*beam.x[:, None])

#Centeralized the beam in XY plane
crossZ=abs(beamParams['dx'])/np.tan(np.radians(beamParams['theta']))
dZ_2nd=beam.z[int(plasmaParams['Nz']-1)]-crossZ
x_shift= abs(beamParams['dx'])*dZ_2nd/crossZ

xpixel= beamParams['X']/beamParams['Nx']
ypixel= beamParams['Y']/beamParams['Ny']

nth_pixel= round((x_shift/xpixel)+beamParams['Nx']/2)
#indexSY=int(beamParams['Ny']/2-1.5*(round(beamParams['waist']/ypixel)))
#indexFY=int(beamParams['Ny']/2+1.5*(round(beamParams['waist']/ypixel)))
indexSY=0
indexFY=beamParams['Ny']
indexSX=int(nth_pixel-2*(round(beamParams['waist']/xpixel)))
indexFX=int(nth_pixel+2*(round(beamParams['waist']/xpixel)))
E_afterPlasma= np.swapaxes(np.asarray(NewBeam[indexSX:indexFX, indexSY:indexFY]), 0, 1)

del NewBeam

#np.save('E_afterPlasma', E_afterPlasma)
#%%
'''
#Interpolate
Int_afterPlasma= abs(E_afterPlasma)
Phase_afterPlasma= np.angle(E_afterPlasma)

#del E_afterPlasma

xI=np.linspace(-Int_afterPlasma.shape[1]/2*xpixel, Int_afterPlasma.shape[1]/2*xpixel, int(Int_afterPlasma.shape[1]))
yI=np.linspace(-Int_afterPlasma.shape[0]/2*ypixel, Int_afterPlasma.shape[0]/2*ypixel, int(Int_afterPlasma.shape[0]))
xF=np.linspace(-15e3, 15e3, 2**10)
yF=np.linspace(-15e3, 15e3, 2**16)
Int_f=interpolate.interp2d(xI, yI, Int_afterPlasma, kind='cubic');
Phase_f=interpolate.interp2d(xI, yI, Phase_afterPlasma, kind='cubic');

Int=Int_f(xF, yF)
Phase=Phase_f(xF, yF)

E_initial=Int*np.exp(1j*Phase)

del Int_afterPlasma
del Phase_afterPlasma
'''
#%%
#np.save('E_initial_30_1016', E_initial)
#%%
Xvalue= E_afterPlasma.shape[1]*xpixel
Yvalue= E_afterPlasma.shape[0]*ypixel
imaging=phaselens.phase_lens_2d(Xvalue, Yvalue, E_afterPlasma, 1e6, 0.8e6, 2.25e6, 0.8)

#%%
xF=np.linspace(-Xvalue/2, Xvalue/2, 4916)
yF=np.linspace(-Yvalue/2, Yvalue/2, 256)
#plt.pcolormesh(xF, yF, (abs(imaging)**2))
plt.pcolormesh((abs(imaging)**2))

plt.axis('scaled')
plt.title('E field @ 2.25 (Imagining at the center of plasma by f=1m lens)')
#%%
plt.pcolormesh(abs(E_afterPlasma))
