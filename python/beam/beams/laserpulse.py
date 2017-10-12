#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:15:01 2017

@author: robert
"""

import pyfftw
import numpy as np
from beam.beams import beam
from beam.calc import laser
import matplotlib.pyplot as plt


class Pulse(beam.Beam):
    """ A laser pulse class that stores the field for each transverse slice.
    
    This class stores a three dimensional grid, two transverse spatial
    directions and a temporal direction. The temporal coordinate is measured so
    that t=0 is at the center of the grid, Nt/2. The temporal component is
    stored in complex notation, it can be used as either an envelope or
    explicit field. For this class, z tracks the location of t=0.
    
    Parameters
    ----------
    Nx : int
        Number of grid points in the x direction, a power of 2 is recommended.
    Ny : int
        Number of grid points in the y direction, a power of 2 is recommended.
    Nt : int
        Number of grid points in the x direction, a power of 2 is recommended.
    X : double
        Width of the grid in the x direction, the grid goes from [-X/2, X/2).
    Y : double
        Width of the grid in the y direction, the grid goes from [-Y/2, Y/2).
    T : double
        Length of the grid in the temporal direction.
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
            'Nt',
            'X',
            'Y',
            'T',
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
        """ Create an x-y rectangular grid and temporal grid. """
        T = self.T
        X = self.X
        Y = self.Y
        self.t = np.linspace(-T/2, T/2, self.Nt, False, dtype='double')
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
        # Temporal fft
        tfft = pyfftw.empty_aligned(self.Nt, dtype='complex128')
        self.tfft = pyfftw.builders.fft(tfft, overwrite_input=True,
                                         avoid_copy=True, threads=threads)
        self.tifft = pyfftw.builders.ifft(tfft, overwrite_input=True, 
                                           avoid_copy=True, threads=threads)
    
    def initialize_field(self, e=None):
        """ Create the array to store the electric field values in.
        
        Parameters
        ----------
        e : array-like, optional
            The array of field values to initialize the field to.
        """
        if e is None:
            self.e = np.zeros((self.Nt, self.Nx, self.Ny,), dtype='complex128')
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
            e = e[:, :, int(self.Ny/2)]
        np.save(self.filePre + '_field_' + str(self.saveInd) + '.npy', e)
        self.saveInd += 1
        self.z.append(z)
        
    def save_z(self):
        """ Save the z array. """
        np.save(self.filePre + '_z.npy', self.z)
        
    def load_field(self, ind):
        """ Load the electric field at the specified index. 
        
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
        # TODO implement this function - maybe include dispersion?
    
    # Visualization functions
    #--------------------------------------------------------------------------
    
    def plot_current_tran_intensity(self):
        """ Plots the current intensity at the center of the pulse. """
        im = self.plot_tran_intensity(self.e[int(self.Nt/2), :, :], self.z[-1])
        plt.show(im)
        
    def plot_tran_intensity_at(self, ind):
        """ Plots the intensity at a particular z distance.
        
        Parameters
        ----------
        ind : int
            The save index to plot the field at, see the _z file to find z.
        """
        e, z = self.load_field(ind)
        if not self.cyl:
            data = e[int(self.Nt/2), :, :]
        else:
            x = self.x
            y = self.y
            data = self.reconstruct_from_cyl(x, e[int(self.Nt/2), :], x, y)
        im = self.plot_tran_intensity(data, z)
        plt.show(im)
    
    def plot_tran_intensity(self, e, z):
        """ Create a transverse intensity plot. """
        X = self.X
        Y = self.Y
        
        I = self.intensity_from_field(e)
        I = self.prep_data(I)
        im = plt.imshow(I, aspect='auto', extent=[-X/2, X/2, -Y/2, Y/2])
        cb = plt.colorbar()
        cb.set_label(r'Intensity ($\mathrm{10^{14}W/cm^2}$)')
        plt.set_cmap('viridis')
        plt.xlabel(r'x')
        plt.ylabel(r'y')
        plt.title('Transverse intensity at z='+str(z))
        return im
    
    def plot_current_long_intensity(self):
        """ Plots the current intensity at the in the x-t plane. """
        im = self.plot_long_intensity(self.e[:, :, int(self.Ny/2)], self.z[-1])
        plt.show(im)
        
    def plot_long_intensity_at(self, ind):
        """ Plots the intensity in x-t at a particular z distance.
        
        Parameters
        ----------
        ind : int
            The save index to plot the field at, see the _z file to find z.
        """
        e, z = self.load_field(ind)
        if not self.cyl:
            e = e[:, :, int(self.Ny/2)]
        im = self.plot_long_intensity(e, z)
        plt.show(im)
    
    def plot_long_intensity(self, e, z):
        """ Create an longitudinal intensity plot. """
        T = self.T
        X = self.X
        
        I = self.intensity_from_field(e)
        I = self.prep_data(I)
        im = plt.imshow(I, aspect='auto', extent=[-T/2, T/2, -X/2, X/2])
        cb = plt.colorbar()
        cb.set_label(r'Intensity ($\mathrm{10^{14}W/cm^2}$)')
        plt.set_cmap('viridis')
        plt.xlabel(r't')
        plt.ylabel(r'x')
        plt.title('Longitudinal intensity at z='+str(z))
        return im


class GaussianPulse(Pulse):
    """ A laser pulse class that creates a Gaussian electric field. 
    
    The pulse is Gaussian in both space and time.
    
    Parameters
    ----------
    E0 : double
        The peak value of the electric field at the Gaussian waist in GV/m. 
    waist : double
        The spot size of the Gaussian waist.
    z : double
        The position relative to the waist to start the beam at. +z is after
        the waist, -z is before the waist.
    tau : double
        The RMS duration of the pulse.
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['E0',
                 'waist',
                 'z0',
                 'tau'])
        super().__init__(params)
    
    def initialize_field(self):
        """ Create the array to store the electric field values in. 
        
        Fills the field array with the field of a Gaussian pulse.
        """
        k = self.k
        w0 = self.waist
        z0 = self.z0
        E0 = self.E0
        t2 = np.reshape(self.t, (self.Nt, 1, 1))**2
        x2 = np.reshape(self.x, (1, self.Nx, 1))**2
        y2 = np.reshape(self.y, (1, 1, self.Ny))**2
        # Calculate all the parameters for the Gaussian beam
        r2 = x2 + y2
        zr = np.pi*w0**2 / self.lam
        if z0 != 0:
            wz = w0 * np.sqrt(1+(z0/zr)**2)
            psi = np.arctan(z0/zr)
            Rz = z0 * (1 + (zr/z0)**2)
            # Create the Gaussian field
            e = E0 * w0 / wz * np.exp(-r2/wz**2-t2*np.pi/(2*self.tau**2)) \
                 * np.exp(1j*(k*z0 + k*r2/(2*Rz) - psi))
        else:
            e = E0 * np.exp(-r2/w0**2-t2*np.pi/(2*self.tau**2))
        super().initialize_field(e)
        