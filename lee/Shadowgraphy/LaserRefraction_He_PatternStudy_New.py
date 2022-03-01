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
from scipy.interpolate import interp1d

'''
If laser plasma interaction fail, check the grid density. If the spacing k*theta is larger than the grid spacig, it will fail.
Try with a smaller crossing angle to confirm the code works. 
'''
#%%
basePath = 'lee/'
path = basePath + 'LaserRefraction/'

power= 2
PlasmaWidthS= 70
# Build the plasma
plasmaParams = {
    'Nx' : 2**13,
    'Ny' : 2**8,#2**8,
    'Nz' : 2**8,
    'X' : 0.110e6,#0.060e6,
    'Y' : 0.040e6,#0.030e6,
    'Z' : 2.0e6, 
    'n0' : 1, #Nominal gas number density in 10^17 cm^-3.
    'atom' : ionization.He,
    'path' : path,
    'name' : 'HePlasma',
    'load' : False,
    'cyl' : True
}


w = PlasmaWidthS
z0 = 0.85e6 #The distance at which the uniform fully ionized plasma starts.
zf = plasmaParams['Z']
dz = 0.3e6 #The length of the fully ionized plasma.
sigmaIn = 0.10e6
sigmaOut = 0.10e6
z, ne = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
He = plasma.Plasma(plasmaParams)

n = np.full((He.Nx, He.Ny, He.Nz), plasmaParams['n0'], dtype='double')
r = np.sqrt( He.x[:, None, None]**2 + He.y[None, :, None]**2)



ne = plasmaParams['n0'] * ne[None, None, :] * np.exp(-r**power/(2*w**power))
He.initialize_plasma(n, ne)
del n
del ne

#He.plot_long_density_center()
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
    'waist' :  16000, #8000,#The spot size of the flattop region.
    'z0' : 0.0,
    'order' : 4,
    'theta' : 1.6,#2.06175769, #The angle of propagation in degrees, positive is angled upward.
    'dx' : -0.02793e6
}

beam = laserbeam.GeneralSuperGaussianLaser(beamParams)
#%%
#beam.plot_current_field()
#beam.plot_current_intensity()
#%%
c= 3e8
epsilon= 8.85e-12
I= abs(beam.e)**2*c*epsilon/2
#%%
interactions.beam_plasma(beam, He)
#plt.figure(4)
#beam.plot_intensity_at(255)
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
plt.pcolormesh(abs(E_afterPlasma))
#%%
#np.save('E_initial_30_1016', E_initial)
'''
#%%Imaging
Physical_X= E_afterPlasma.shape[1]*xpixel
Physical_Y= E_afterPlasma.shape[0]*ypixel
imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, E_afterPlasma, 0.5e6, 0.8e6, 0.6923e6, beamParams['lam'])
#%%
Iimg= abs(imaging)**2*sum(sum(I))/sum(sum(abs(imaging)**2))
#%%

plt.figure(3)
xF=np.linspace(-Physical_X/2, Physical_X/2, imaging.shape[1])
yF=np.linspace(-Physical_Y/2, Physical_Y/2, imaging.shape[0])
plt.pcolormesh(xF, yF, Iimg, cmap='viridis')
plt.axis('scaled')
plt.title('Normalized Intensity (Imagining at the center of plasma by f=1m lens) \n n= 1e17 w=200$\mu$m')    
plt.xlabel('x $\mu$m')
plt.ylabel('y $\mu$m')
plt.colorbar(label='10^14 W/cm^2')

#%%
lineout= abs(Iimg[:,int(imaging.shape[1]/2)])
f= interp1d(yF, lineout, kind='cubic')
Newgrid= 10000
newY= np.linspace(-Physical_Y/2, Physical_Y/2, Newgrid)
LONew= f(newY)
'''
#%%
'''
CIdx= int(LONew.shape[0]/2)
MaxIdx= int(LONew.shape[0])
for n in range (CIdx, int(MaxIdx-250)):
    if LONew[n+1]-LONew[n] <0:
        if LONew[n+50]-LONew[n] <0:
            if LONew[n+100]-LONew[n] <0:            
                if LONew[n+150]-LONew[n] <0:            
                    if LONew[n+250]-LONew[n] <0:            
                        if LONew[n-1]-LONew[n] <0:            
                            if LONew[n-50]-LONew[n] <0:            
                                if LONew[n-100]-LONew[n] <0:            
                                    if LONew[n-150]-LONew[n] <0:            
                                        if LONew[n-250]-LONew[n] <0:            
                                            print(n)
                                            PeakIdxH= n
                                            break
for n in range (CIdx, 250, -1):
    if LONew[n+1]-LONew[n] <0:
        if LONew[n+50]-LONew[n] <0:
            if LONew[n+100]-LONew[n] <0:            
                if LONew[n+150]-LONew[n] <0:            
                    if LONew[n+250]-LONew[n] <0:            
                        if LONew[n-1]-LONew[n] <0:            
                            if LONew[n-50]-LONew[n] <0:            
                                if LONew[n-100]-LONew[n] <0:            
                                    if LONew[n-150]-LONew[n] <0:            
                                        if LONew[n-250]-LONew[n] <0:            
                                            print(n)
                                            PeakIdxL= n
                                            break
#%%
PatternL= PeakIdxH-PeakIdxL
pixel= Physical_Y/Newgrid
PeakSpacing= PatternL*pixel
PeakHeigh= LONew[PeakIdxH]-LONew[CIdx]
#%%
plt.plot(newY, LONew)
plt.xlabel('Vertical Lineout $\mu$m')
plt.ylabel('Intensity 10^14 W/cm^2')
plt.title('Plasma W=200')
#%%
#W30= PeakSpacing
#H30= PeakHeigh
#W40= PeakSpacing
#H40= PeakHeigh
#W50= PeakSpacing
#H50= PeakHeigh
#W60= PeakSpacing
#H60= PeakHeigh
#W80= PeakSpacing
#H80= PeakHeigh
#W100= PeakSpacing
#H100= PeakHeigh
#W120= PeakSpacing
#H120= PeakHeigh
#W140= PeakSpacing
#H140= PeakHeigh
#W160= PeakSpacing
#H160= PeakHeigh
#%%
np.save('W200', LONew)

#%%
WHMatrix= np.zeros((9, 2))
WHMatrix[0][0]= W30
WHMatrix[0][1]= H30
WHMatrix[1][0]= W40
WHMatrix[1][1]= H40
WHMatrix[2][0]= W50
WHMatrix[2][1]= H50
WHMatrix[3][0]= W60
WHMatrix[3][1]= H60
WHMatrix[4][0]= W80
WHMatrix[4][1]= H80
WHMatrix[5][0]= W100
WHMatrix[5][1]= H100
WHMatrix[6][0]= W120
WHMatrix[6][1]= H120
WHMatrix[7][0]= W140
WHMatrix[7][1]= H140
WHMatrix[8][0]= W160
WHMatrix[8][1]= H160
#%%
np.save('WHMatrix', WHMatrix)
#%%
PlasmaW=np.array([30, 40, 50, 60, 80, 100, 120, 140, 160]) 
plt.plot(PlasmaW, WHMatrix[:, 0], 'o')
plt.xlabel('Plasma Width $\mu$m')
#plt.ylabel('Peak-Valley W/cm^2')
plt.ylabel('FirstPeak Spacing $\mu$m')
'''