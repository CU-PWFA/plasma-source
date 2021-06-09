#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 17:05:18 2020

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
from multiprocessing import Pool
'''
If laser plasma interaction fail, check the grid density. If the spacing k*theta is larger than the grid spacig, it will fail.
Try with a smaller crossing angle to confirm the code works. 
'''
#%%
def ParameterScan(GPList):
    GP_begining= GPList[0]
    GP_end= GPList[1]
    for power in range (int(GP_begining), int(GP_end)):
    
        for PlasmaWidthS in range (20, 40, 10):
    
            basePath = 'lee/'
            path = basePath + 'LaserRefraction/'
            
            plasmaParams = {
                    'Nx' : 2**12,
                    'Ny' : 2**9,#2**8,
                    'Nz' : 2**8,
                    'X' : 0.080e6,#0.060e6,
                    'Y' : 0.060e6,#0.030e6,
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
            sigmaIn = 0.15e6
            sigmaOut = 0.15e6
            z, ne = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
            He = plasma.Plasma(plasmaParams)
            
            n = np.full((He.Nx, He.Ny, He.Nz), plasmaParams['n0'], dtype='double')
            r = np.sqrt( He.x[:, None, None]**2 + He.y[None, :, None]**2)
    
    
    
            ne = plasmaParams['n0'] * ne[None, None, :] * np.exp(-r**power/(2*w**power))
    
            He.initialize_plasma(n, ne)
            del n
            del ne
    
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
                    'order' : 3,
                    'theta' : 0.859,#2.06175769, #The angle of propagation in degrees, positive is angled upward.
                    'dx' : -0.015e6
                    }
    
            beam = laserbeam.GeneralSuperGaussianLaser(beamParams)
    
            c= 3e8
            epsilon= 8.85e-12
            I= abs(beam.e)**2*c*epsilon/2
            interactions.beam_plasma(beam, He)
    
            NewBeam = beam.e*np.exp(-1j*beam.k*np.radians(beam.theta)*beam.x[:, None])
    
            x_final= plasmaParams['Z']*np.tan(np.radians(beamParams['theta']))+beamParams['dx']
            xpixel= beamParams['X']/beamParams['Nx']
            ypixel= beamParams['Y']/beamParams['Ny']
    
            nth_pixel= round((x_final/xpixel)+beamParams['Nx']/2)
    
            indexSY=0
            indexFY=beamParams['Ny']
            indexSX=int(nth_pixel-1.5*(round(beamParams['waist']/xpixel)))
            indexFX=int(nth_pixel+1.5*(round(beamParams['waist']/xpixel)))
            E_afterPlasma= np.swapaxes(np.asarray(NewBeam[indexSX:indexFX, indexSY:indexFY]), 0, 1)
            
            del NewBeam
    
            Physical_X= E_afterPlasma.shape[1]*xpixel
            Physical_Y= E_afterPlasma.shape[0]*ypixel
            imaging=phaselens.phase_lens_2d(Physical_X, Physical_Y, E_afterPlasma, 1e6, 0.8e6, 2.25e6, beamParams['lam'])
            Iimg= abs(imaging)**2*sum(sum(I))/sum(sum(abs(imaging)**2))
            lineout= abs(Iimg[:,int(imaging.shape[1]/2)])
            xF=np.linspace(-Physical_X/2, Physical_X/2, imaging.shape[1])
            yF=np.linspace(-Physical_Y/2, Physical_Y/2, imaging.shape[0])
    
            f= interp1d(yF, lineout, kind='cubic')
            Newgrid= 10000
            newY= np.linspace(-Physical_Y/2, Physical_Y/2, Newgrid)
            LONew= f(newY)
    
    
            np.save('PW'+str(PlasmaWidthS)+'_GP'+str(power), LONew)    
    
if __name__ == '__main__':
    with Pool(4) as p:
        p.map(ParameterScan, [[1, 2], [2, 3]])
