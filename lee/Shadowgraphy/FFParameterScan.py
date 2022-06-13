#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 22:55:11 2021

@author: valentinalee
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
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

'''
If laser plasma interaction fail, check the grid density. If the spacing k*theta is larger than the grid spacig, it will fail.
Try with a smaller crossing angle to confirm the code works. 
'''
#%%
c= 3e8
epsilon= 8.85e-12

laserPara= np.load('laser_para.npy')
FFphase=np.load('FFphase.npy')
PhaseGrid=np.load('xgrid.npy')
#%%
basePath = 'lee/'
path = basePath + 'LaserRefraction/'

power= 2
PlasmaWidthS= 70
# Build the plasma
plasmaParams = {
    'Nx' : 2**16, #2**17,
    'Ny' : 2**7,
    'Nz' : 2**10,
    'X' : 0.14e6,#0.060e6,
    'Y' : 0.006e6,#0.030e6,
    'Z' : 1.5e6, 
    'n0' : 1, #Nominal gas number density in 10^17 cm^-3.
    'atom' : ionization.He,
    'path' : path,
    'name' : 'HePlasma',
    'load' : False,
    'cyl' : True
}

w = PlasmaWidthS
z0 = 0.6e6 #The distance at which the uniform fully ionized plasma starts.
zf = plasmaParams['Z']
dz = 0.3e6 #The length of the fully ionized plasma.
sigmaIn = 0.15e6
sigmaOut = 0.15e6
z, nez = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
He = plasma.Plasma(plasmaParams)
r = np.sqrt(He.x[:, None]**2 + He.y[None, :]**2)

#%%
os.chdir("lee/LaserRefraction/elements/element_HePlasma")
for z in range(plasmaParams['Nz']):
    print(z)
    ne = plasmaParams['n0'] * nez[z] * np.exp(-r**power/(2*w**power))
    ne= ne[:, int(plasmaParams['Ny']/2)]
    np.save(str(plasmaParams['name']) + '_plasmaDensity_' + str(z) + '.npy', ne)

os.chdir("../../../../")
    

#%%
del nez
del ne

#%%
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
    'waist' : laserPara[1]*1e6,
    'z0' : 0.0,
    'order' : laserPara[2],
    'theta' : 5, #The angle of propagation in degrees, positive is angled upward.
    'dx' : -0.0656e6
}

#beam = laserbeam.GeneralGaussianLaser(beamParams)
beam = laserbeam.GeneralSuperGaussianLaser(beamParams)

#%%
beam.plot_current_intensity()

#%%
xg= PhaseGrid
yg= PhaseGrid
Xg, Yg= np.meshgrid(xg, yg)

xgrid= np.linspace(-beamParams['X']/2, beamParams['X']/2, beamParams['Nx'], dtype='double')
ygrid= np.linspace(-beamParams['Y']/2, beamParams['Y']/2, beamParams['Ny'], dtype='double')
Xgrid, Ygrid= np.meshgrid(xgrid, ygrid)
f= interp2d(xg, yg, FFphase, kind='cubic')
phase2= f(xgrid, ygrid)
#%%
interactions.beam_phase(beam, phase2)
#%%
I= abs(beam.e)**2*c*epsilon/2
interactions.beam_plasma(beam, He)
#%%
beam.plot_intensity_at(1023)

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

#%%
#%%Imaging
Physical_X= E_afterPlasma.shape[1]*xpixel
Physical_Y= E_afterPlasma.shape[0]*ypixel
#imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, E_afterPlasma, 1e6, 2.25e6, 1.5e6, beamParams['lam'])
#imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, E_afterPlasma, 0.5e6, 0.7115e6, 0e6, beamParams['lam'])
imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, E_afterPlasma, 0.5e6, 1.4615e6, 0.76e6, beamParams['lam'])
#%%
Iimg= abs(imaging)**2*sum(sum(I))/sum(sum(abs(imaging)**2))
#%%
plt.figure(3)
xF=np.linspace(-Physical_X/2, Physical_X/2, imaging.shape[1])
yF=np.linspace(-Physical_Y/2, Physical_Y/2, imaging.shape[0])
plt.pcolormesh(xF, yF, Iimg, cmap='viridis')
plt.axis('scaled')
plt.title('Normalized Intensity (Imagining at the end of the plasma by f=1m lens) \n n= 1e17 w=70$\mu$m')    
plt.xlabel('x $\mu$m')
plt.ylabel('y $\mu$m')
plt.colorbar(label='10^14 W/cm^2')

#%%
lineout= abs(Iimg[:,int(imaging.shape[1]/2)])
f= interp1d(yF, lineout, kind='cubic')
Newgrid= 10000
newY= np.linspace(-Physical_Y/2, Physical_Y/2, Newgrid)
LONew= f(newY)

#%%
plt.plot(newY, LONew)
plt.xlabel('$\mu m$')