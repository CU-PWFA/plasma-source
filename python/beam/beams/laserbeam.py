#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:47:12 2017

@author: robert
"""

import os
import pyfftw
import numpy as np
from beam.beams import beam
from beam.calc import laser


class Laser(beam.Beam):
    """ A laser beam class that stores the field on a two dimensional grid. """
    keys = ['Nx',
            'Ny',
            'X',
            'Y',
            'lam',
            'path',
            'name',
            'threads']
    
    def __init__(self, params):
        super().__init__(params)
        self.k = 2*np.pi / self.params['lam']
        # Create a folder to store the beam data in
        self.dirName = dirName = self.path + 'beams/beam_' + self.name + '/'
        self.filePre = dirName + self.name
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        self.clear_dir()
        # Create internal variables
        self.create_grid()
        self.create_fft()
        self.initialize_field()
        self.save_initial()        
    
    def create_grid(self):
        """ Create an x-y rectangular grid. """
        X = self.X
        Y = self.Y
        self.x = np.linspace(-X/2, X/2, self.Nx, False, dtype='double')
        self.y = np.linspace(-Y/2, Y/2, self.Ny, False, dtype='double')
        self.z = []
    
    def create_fft(self):
        """ Create the fftw plans. """
        threads = self.threads
        # Allocate space to carry out the fft's in
        efft = pyfftw.empty_aligned((self.Nx, self.Ny), dtype='complex128')
        self.fft = pyfftw.builders.fft2(efft, overwrite_input=True,
                                         avoid_copy=True, threads=threads)
        self.ifft = pyfftw.builders.ifft2(efft, overwrite_input=True, 
                                           avoid_copy=True, threads=threads)
        
    def initialize_field(self, e=None):
        """ Create the array to store the electric field values in. """
        if e is None:
            self.e = np.zeros((self.Nx, self.Ny), dtype='complex128')
        else:
            self.e = e
        self.saveInd = 0
        self.z = []
        self.save_field(self.e, 0.0)
        
    def set_field(self, e):
        """ Set the value of the electric field. """
        self.e = np.array(e, dtype='complex128')
        self.save_field(self.e, self.z[-1])
        
    def save_initial(self):
        """ Save the initial params object and the grid. """
        filePre = self.filePre
        np.save(filePre + '_params.npy', self.params)
        np.save(filePre + '_x.npy', self.x)
        np.save(filePre + '_y.npy', self.y)
    
    def save_field(self, e, z):
        """ Save the current electric field to file and adavnce z. """
        np.save(self.filePre + '_field_' + str(self.saveInd) + '.npy', e)
        self.saveInd += 1
        self.z.append(z)
        np.save(self.filePre + '_z.npy', self.z)
    
    def clear_dir(self):
        """ Clear all files from the beam directory. """ 
        # TODO implemet this function
        
    def propagate(self, z, n):
        """ Propagate the field to an array of z distances.
        
        Prameters
        ---------
        z : array-like
            Array of z distances from the current z to calculate the field at. 
            Does not need to be evenly spaced.
        n : double
            Index of refraction of the medium the wave is propagating through.
        """
        z = np.array(z, ndmin=1, dtype='double')
        self.e = laser.fourier_prop(self.e, self.x, self.y, z, self.lam, n, 
                                    self.fft, self.ifft, self.save_field)


class GaussianLaser(Laser):
    """ A laser beam class that creates a Gaussian electric field. """
    
    def __init__(self, params):
        self.keys.extend(
                ['E0',
                 'waist',
                 'z'])
        super().__init__(params)
    
    def initialize_field(self):
        """ Create the array to store the electric field values in. 
        
        Fills the field array with the field of a Gaussian pulse.
        """
        k = self.k
        w0 = self.params['waist']
        z = self.params['z']
        E0 = self.params['E0']
        x2 = np.reshape(self.x, (self.params['Nx'], 1))**2
        y2 = np.reshape(self.y, (1, self.params['Ny']))**2
        # Calculate all the parameters for the Gaussian beam
        r2 = x2 + y2
        zr = np.pi*w0**2 / self.params['lam']
        wz = w0 * np.sqrt(1+(z/zr)**2)
        Rz = z * (1 + (zr/z)**2)
        psi = np.arctan(z/zr)
        # Create the Gaussian field
        e = E0 * w0 / wz * np.exp(-r2/wz**2) \
                 * np.exp(-1j*(k*z + k*r2/(2*Rz) - psi))
        super().initialize_field(e)


class SuperGaussianLaser(Laser):
    # TODO create a super Gaussian beam
    """ A laser beam class that creates a super-Gaussian electric field. """
    
    def __init__(self, params):
        self.keys.extend(
                ['E0',
                 'waist',
                 'order'])
        super().__init__(params)
    
    def initialize_field(self):
        """ Create the array to store the electric field values in. 
        
        Fills the field array with the field of a Gaussian pulse.
        """
        w0 = self.params['waist']
        E0 = self.params['E0']
        n = self.params['order']
        x2 = np.reshape(self.x, (self.params['Nx'], 1))**2
        y2 = np.reshape(self.y, (1, self.params['Ny']))**2
        # Calculate all the parameters for the Gaussian beam
        r = np.sqrt(x2 + y2)
        e = E0 * np.exp(-(r/w0)**n)
        super().initialize_field(e)
