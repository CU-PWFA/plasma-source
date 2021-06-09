#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:23:06 2020

@author: valentina_lee
"""

#setup cython
#import pyximport
#pyximport.install()

import sys
import os
sys.path.insert(0, "python")
import numpy as np
from ionization import ionization
from lens import profile
import matplotlib.pyplot as plt
from beam.beams import laserpulse
from beam.beams import laserbeam
from beam.elements import plasma
from beam.elements import optic
from beam import interactions
from lens import design
from scipy.interpolate import interp1d
#%%
#%
#plt.style.use('notes')
#%load_ext autoreload
#%autoreload 2
Eng= [2.5, 3, 3.5, 4]
for energy in Eng:
    
    basePath = 'lee/'
#    path = basePath + 'AxiconFormedPlasma/'
    path = basePath + 'AxiconFormedPlasma/'+str(energy)+'/'
    
    lam = 0.796
    
    # All length units are in umr
    
    # Create a laser beam, pulse and axicon lens
    Nx = 2**13
    X = 60e3
    tau = 30
    pulseParams = {'Nx' : Nx,
                   'Ny' : Nx,
                   'Nt' : 2**6,
                   'X' : X,
                   'Y' : X,
                   'T' : 3*tau,
                   'lam' : lam,
                   'path' : path,
                   'name' : 'FlattopBeam',
                   'threads' : 20,
                   'cyl' : True,
                   'tau' : tau,
                   'load' : False}
    
    beam = laserbeam.Laser(pulseParams)
    
    alpha = 0.5
    n = 1.458
    beta = np.degrees(np.arcsin(n*np.sin(np.radians(alpha)))) - alpha
    print('Axicon deflection angle : %0.3f deg' % beta)
    lensParams = {'Nx' : Nx,
                  'Ny' : Nx,
                  'X' : X,
                  'Y' : X,
                  'path' : path,
                  'name' : 'Axicon',
                  'lam' : lam,
                  'beta' : beta,
                  'load' : False}
    axicon = optic.AxiconLens(lensParams)
    
    apertureParams = {'Nx' : pulseParams['Nx'],
                      'Ny' : pulseParams['Ny'],
                      'X' : pulseParams['X'],
                      'Y' : pulseParams['Y'],
                      'path' : path,
                      'name' : 'Aperturen',
                      'lam' : lam,
                      'r' : 50.8e3/2,
                      'load' : False}
    
    aperture = optic.Aperture(apertureParams)
    
    #%
    # Create a super gaussian intensity profile and pass it through a lens.
    w0 = 15e3
    E0 = energy
    n = 8
    x2 = np.reshape(beam.x, (beam.Nx, 1))**2
    y2 = np.reshape(beam.y, (1, beam.Ny))**2
    r = np.sqrt(x2 + y2)
    e = E0 * np.exp(-(r/w0)**n)
    
    beam.initialize_field(e)
    #beam.plot_current_intensity()
    print('Power:', beam.total_cyl_power(beam.x[int(beam.Nx/2):],
                                             beam.intensity_from_field(beam.e[int(beam.Nx/2):, int(beam.Ny/2)])))
    interactions.beam_phase(beam, axicon)
    interactions.beam_intensity(beam, aperture)
    #beam.plot_current_intensity()
    r = beam.x[int(Nx/2):]
    E = beam.e[int(Nx/2):, int(Nx/2)]
    del beam
    #%
    # Create the gas density the laser is going into
    ne0 = 1 #*1e17
    z = np.linspace(0, 4e6, 10000)
    sim_start = 0#30e4
    n0= np.zeros(z.shape)+ne0
    
    n= interp1d(z, n0)
    #sim_start, n_plot, n = profile.lithium_oven_profile(z, 40e4+sim_start, ne0)
    sim_length = 380e4#3e6 #80e4
    #np.save(path+'sim_size.npy', [sim_start, sim_length])
    
    
    X = 30e3
    Nx = 2**11
    beam0, pulseParams = design.propagate_to_start(r, E, sim_start, X, Nx, path, lam, tau, 20, [-1, 1])
    
    #Nx = 2**10
    #Nz = 100
    #X =30e3
    #design.domain_test(X, Nx, sim_length, Nz, beam0, pulseParams, z, np.zeros(len(z)), sim_start, [-800, 800]);
    
    Nx = 2**11
    Nz = 100
    ext = [0, sim_length/1e4, -X/2, X/2]
    pulse, I, ne = design.plasma_refraction(X, Nx, sim_length, Nz, beam0, pulseParams, ionization.He, n, sim_start, 1, ne0, Eq=E0)
    
#    design.plot_plasma_density(pulse, ne, ne0, ext, lines=[250, 300, 350])
    
    np.save('ne_'+str(E0), ne)
    np.save('ne0_'+str(E0), ne0)
    np.save('ext_'+str(E0), ext)