#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:47:12 2017

@author: robert
"""

import pyfftw
import numpy as np
from beam.beams import beam
from beam.calc import laser
import matplotlib.pyplot as plt


class Laser(beam.Beam):
    """ A laser beam class that stores the field on a two dimensional grid. 
    
    This class stores a two dimensional grid upon which a complex scalar field
    is stored. The grid has a point at exactly 0, (Nx/2, Ny/2). The grid is
    meant for use with the discrete fourier transform. A note on units, any
    unit of length may be used as long as it is consistent for all variables.
    Generally, units of microns are appropriate.
    
    Parameters
    ----------
    Nx : int
        Number of grid points in the x direction, a power of 2 is recommended.
    Ny : int
        Number of grid points in the y direction, a power of 2 is recommended.
    X : double
        Width of the grid in the x direction, the grid goes from [-X/2, X/2).
    Y : double
        Width of the grid in the y direction, the grid goes from [-Y/2, Y/2).
    lam : double
        The vacuum wavelength of the laser radiation.
    path : string
        The path for the calculation. This class will create a folder inside
        the path to store all output data in.
    name : string
        The name of the beam, used for naming files and folders.
    threads : int
        The number of processors to parallelize the fft over.
    cyl : bool
        Whether the beam is cylindrically symmetric or not. Controls whether
        the entire transverse field is saved or only a 1D slice.
    """
    keys = ['Nx',
            'Ny',
            'X',
            'Y',
            'lam',
            'path',
            'name',
            'threads',
            'cyl']
    
    # Initialization functions
    #--------------------------------------------------------------------------
    
    def __init__(self, params):
        super().__init__(params)
        self.k = 2*np.pi / self.params['lam']
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
        """ Create the array to store the electric field values in. 
        
        Parameters
        ----------
        e : array-like, optional
            The array of field values to initialize the field to.
        """
        if e is None:
            self.e = np.zeros((self.Nx, self.Ny), dtype='complex128')
        else:
            self.e = e
        self.saveInd = 0
        self.z = []
        self.save_field(self.e, 0.0)
        
    # Getters and setters
    #--------------------------------------------------------------------------
        
    def set_field(self, e):
        """ Set the value of the electric field. """
        self.e = np.array(e, dtype='complex128')
        self.save_field(self.e, self.z[-1])
        
    # File managment
    #--------------------------------------------------------------------------
        
    def save_initial(self):
        """ Save the initial params object and the grid. """
        super().save_initial()
        np.save(self.filePre + '_x.npy', self.x)
        np.save(self.filePre + '_y.npy', self.y)
    
    def save_field(self, e, z):
        """ Save the current electric field to file and adavnce z. """
        if self.cyl:
            e = e[:, int(self.Ny/2)]
        np.save(self.filePre + '_field_' + str(self.saveInd) + '.npy', e)
        self.saveInd += 1
        self.z.append(z)
        
    def save_z(self):
        """ Save the z array. """
        np.save(self.filePre + '_z.npy', self.z)
        
    def load_field(self, ind):
        """ load the electric field at the specified index. 
        
        Parameters
        ----------
        ind : int
            The save index to load the field at.
        
        Returns
        -------
        e : array-like
            The electric field at the specified index.
        z : double
            The z coordinate of the field.
        """
        e = np.load(self.filePre + '_field_' + str(ind) + '.npy')
        z = self.z[ind]
        return e, z
        
    # Physics functions
    #--------------------------------------------------------------------------
        
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
        self.e = np.array(self.e, dtype='complex128')
    
    # Visualization functions
    #--------------------------------------------------------------------------
    
    def plot_current_intensity(self):
        """ Plots the current intensity of the beam. """
        im = self.plot_intensity(self.e, self.z[-1])
        plt.show(im)
        
    def plot_intensity_at(self, ind):
        """ Plots the intensity at a particular z distance.
        
        Parameters
        ----------
        ind : int
            The save index to plot the field at, see the _z file to find z.
        """
        e, z = self.load_field(ind)
        if not self.cyl:
            im = self.plot_intensity(e, z)
            plt.show(im)
        # TODO add in plots for cyl beams
    
    def plot_intensity(self, e, z):
        """ Create an intensity plot. """
        X = self.X
        Y = self.Y
        
        I = self.intensity_from_field(e)
        I = self.prep_data(I)
        im = plt.imshow(I, aspect='auto', extent=[-X/2, X/2, -Y/2, Y/2])
        cb = plt.colorbar()
        cb.set_label(r'Intensity')
        plt.set_cmap('viridis')
        plt.xlabel(r'x')
        plt.ylabel(r'y')
        plt.title('Transverse intensity at z='+str(z))
        return im
        


class GaussianLaser(Laser):
    """ A laser beam class that creates a Gaussian electric field. 
    
    Parameters
    ----------
    E0 : double
        The peak value of the electric field at the Gaussian waist. 
    waist : double
        The spot size of the Gaussian waist.
    z : double
        The position relative to the waist to start the beam at. +z is after
        the waist, -z is before the waist.
    """
    
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
                 * np.exp(1j*(k*z + k*r2/(2*Rz) - psi))
        super().initialize_field(e)


class SuperGaussianLaser(Laser):
    """ A laser beam class that creates a super-Gaussian electric field. 
    
    Parameters
    ----------
    E0 : double
        The peak value of the electric fieldon the flattop. 
    waist : double
        The spot size of the flattop region.
    order : int
        The order of the super Gaussian.
    """
    
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
