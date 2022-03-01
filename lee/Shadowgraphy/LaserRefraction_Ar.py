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
#import pyximport; pyximport.install()
import numpy as np
from beam.beams import laserbeam
from beam.elements import plasma
from beam import interactions
from ionization import ionization
from lens import profile
from lens import phaselens
from lens import ray
import matplotlib.pyplot as plt
from propagation import laser
from scipy import interpolate

'''
If laser plasma interaction fail, check the grid density. If the spacing k*theta is larger than the grid spacig, it will fail.
Try with a smaller crossing angle to confirm the code works. 
'''
#%%
basePath = 'lee/'
path = basePath + 'LaserRefraction/'

# Build the plasma
plasmaParams = {
    'Nx' : 2**14,
    'Ny' : 2**8,
    'Nz' : 2**9,
    'X' : 0.120e6,
    'Y' : 0.030e6,
    'Z' : 2.0e6, 
    'n0' : 0.25, #Nominal gas number density in 10^17 cm^-3.
    'atom' : ionization.Ar,
    'path' : path,
    'name' : 'Ar',
    'load' : False,
    'cyl' : True
}


w = 100
z0 = 0.5e6 #The distance at which the uniform fully ionized plasma starts.
zf = plasmaParams['Z']
dz = 1e6 #The length of the fully ionized plasma.
sigmaIn = 0.02e6
sigmaOut = 0.02e6
z, ne = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
Ar = plasma.Plasma(plasmaParams)

n = np.full((Ar.Nx, Ar.Ny, Ar.Nz), plasmaParams['n0'], dtype='double')
r2 = Ar.x[:, None, None]**2 + Ar.y[None, :, None]**2



ne = plasmaParams['n0'] * ne[None, None, :] * np.exp(-(r2/(2*w**2)))
Ar.initialize_plasma(n, ne)
del n
del ne

Ar.plot_long_density_center()
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
    'theta' : 2.2906,#The angle of propagation in degrees, positive is angled upward.
    'dx' : -0.04e6
}

beam = laserbeam.GeneralSuperGaussianLaser(beamParams)
#%%
beam.plot_current_field()
beam.plot_current_intensity()
#%%
c= 3e8
epsilon= 8.85e-12
I= abs(beam.e)**2*c*epsilon/2
#%%
interactions.beam_plasma(beam, Ar)
#%%
NewBeam = beam.e*np.exp(-1j*beam.k*np.radians(beam.theta)*beam.x[:, None])
#%%
Lineout= abs(NewBeam[13653, :])**2
#%%
plt.figure()
beam.plot_intensity_at(511)
#%%
plt.figure()
plt.plot(abs(NewBeam[13653, :])**2)
#%%
#fix propagation angle
NewBeam = beam.e*np.exp(-1j*beam.k*np.radians(beam.theta)*beam.x[:, None])

#Centeralized the beam in XY plane
x_final= plasmaParams['Z']*np.tan(np.radians(beamParams['theta']))+beamParams['dx']
xpixel= beamParams['X']/beamParams['Nx']
ypixel= beamParams['Y']/beamParams['Ny']

nth_pixel= round((x_final/xpixel)+beamParams['Nx']/2)

#indexSY=int(beamParams['Ny']/2-1.5*(round(beamParams['waist']/ypixel)))
#indexFY=int(beamParams['Ny']/2+1.5*(round(beamParams['waist']/ypixel)))
indexSY=0
indexFY=beamParams['Ny']
indexSX=int(nth_pixel-1.5*(round(beamParams['waist']/xpixel)))
indexFX=int(nth_pixel+1.5*(round(beamParams['waist']/xpixel)))
E_afterPlasma= np.swapaxes(np.asarray(NewBeam[indexSX:indexFX, indexSY:indexFY]), 0, 1)

del NewBeam

#np.save('E_initial_30_1016', E_initial)
#%%Imaging
Physical_X= E_afterPlasma.shape[1]*xpixel
Physical_Y= E_afterPlasma.shape[0]*ypixel
imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, E_afterPlasma, 0.5e6, 0.75e6, 0.7e6, beamParams['lam'])
#%%
Iimg= abs(imaging)**2*sum(sum(I))/sum(sum(abs(imaging)**2))
#%%

plt.figure(3)
xF=np.linspace(-Physical_X/2, Physical_X/2, imaging.shape[1])
yF=np.linspace(-Physical_Y/2, Physical_Y/2, imaging.shape[0])
plt.pcolormesh(xF, yF, Iimg, cmap='viridis')
plt.axis('scaled')
plt.title('Normalized Intensity (Imagining at the center of plasma by f=1m lens) \n n= 8e17')    
plt.xlabel('x $\mu$m')
plt.ylabel('y $\mu$m')
plt.colorbar(label='10^14 W/cm^2')

#%%
#Lineout1e16= abs(Iimg[:,int(imaging.shape[1]/2)])
#Lineout5e16= abs(Iimg[:,int(imaging.shape[1]/2)])
#Lineout1e17= abs(Iimg[:,int(imaging.shape[1]/2)])
#Lineout5e17= abs(Iimg[:,int(imaging.shape[1]/2)])
Lineout2_5e16= abs(Iimg[:,int(imaging.shape[1]/2)])


#%%

cpixel= round(beamParams['Nx']/2+beamParams['dx']/xpixel)
ILineout= I[cpixel, :]
SB1e16= ILineout- Lineout1e16
SB5e16= ILineout- Lineout5e16
SB1e17= ILineout- Lineout1e17
SB5e17= ILineout- Lineout5e17

#%%
#plt.plot(yF, ILineout, label= 'Initial Probe')
#plt.plot(yF, Lineout1e16, label='1e16')
#plt.plot(yF, Lineout5e16, label='5e16')
#plt.plot(yF, Lineout1e17, label='1e17')
#plt.plot(yF, Lineout5e17, label='5e17')
plt.figure(5)
plt.plot(yF, Lineout2_5e16, label='2.6e16')

plt.legend()
plt.xlabel('Verticle Lineout $\mu$m')
plt.ylabel('Normalized Intensity 10^14 W/cm^2')


#%%

plt.plot(yF, ILineout, label= 'Initial Probe')
plt.plot(yF, SB1e16, label='1e16')
plt.plot(yF, SB5e16, label='5e16')
plt.plot(yF, SB1e17, label='1e17')
#plt.plot(yF, SB5e17, label='5e17')

plt.legend()
plt.xlabel('Verticle Lineout $\mu$m')
plt.ylabel('Normalized Intensity 10^14 W/cm^2')
'''