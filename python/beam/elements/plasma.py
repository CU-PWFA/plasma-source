#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 15:08:43 2017

@author: robert
"""

import numpy as np
from beam.elements import element
import matplotlib.pyplot as plt


class Plasma(element.Element):
    """ A class that represents a gas/plasma that a laser can ionize. 
    
    This class stores the gas and plasma density for a region where a laser
    pulse can pass through. It contains information about the gas density and
    the state of ionization.
    
    Parameters
    ----------
    Nx : int
        Number of grid points in the x direction, must match beam.
    Ny : int
        Number of grid points in the y direction, must match beam.
    Nz : int
        Number of grid points in the z direction.
    X : int
        Width of the grid in the x direction, the grid goes from [-X/2, X/2).
        This must match the beam parameters.
    Y : int
        Width of the grid in the y direction, the grid goes from [-Y/2, Y/2).
        This must match the beam parameters.
    Z : int
        Length of the grid in the z direction, the grid goes from[0, Z].
    n0 : double
        Nominal gas number density in 10^17 cm^-3.
    atom : dictionary
        EI : double
            Ionization energy in eV.
        Z : double
            Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
        l : double
            Orbital quantum number of the electron being ionized.
        m : double
            Magnetic quantum number of the electron being ionized.
        alpha : double
            Atomic polarizability of the gas in A^3.
    path : string
        The path for the calculation. This class will create a folder inside
        the path to store all output data in.
    name : string
        The name of the beam, used for naming files and folders.
    cyl : bool
        Whether the plasma is cylindrically symmetric or not. Controls whether
        the entire transverse density is saved or only a 1D slice.
    """
    keys = ['Nx',
            'Ny',
            'Nz',
            'X',
            'Y',
            'Z',
            'n0',
            'atom',
            'path',
            'name',
            'cyl']
    
    # Initialization functions
    #--------------------------------------------------------------------------
    
    def __init__(self, params):
        super().__init__(params)
        self.create_grid()
        self.save_initial()
        
    def create_grid(self):
        """ Create an x-y rectangular grid. """
        X = self.X
        Y = self.Y
        Z = self.Z
        self.x = np.linspace(-X/2, X/2, self.Nx, False, dtype='double')
        self.y = np.linspace(-Y/2, Y/2, self.Ny, False, dtype='double')
        self.z = np.linspace(0.0, Z, self.Nz, dtype='double')
    
    #File managment
    #--------------------------------------------------------------------------
    
    def save_initial(self):
        """ Save the initial params object and the grid. """
        super().save_initial()
        np.save(self.filePre + '_x.npy', self.x)
        np.save(self.filePre + '_y.npy', self.y)
        np.save(self.filePre + '_z.npy', self.z)
    
    def save_plasma_density(self, ne, ind):
        """ Save the plasma density to file at the given z ind. """
        if self.cyl:
            ne = ne[:, int(self.Ny/2)]
        np.save(self.filePre + '_plasmaDensity_' + str(ind) + '.npy', ne)
        
    def load_plasma_density(self, ind):
        """ Load the plasma density at the specified index. 
        
        Parameters
        ----------
        ind : int
            The save index to load the field at.
        
        Returns
        -------
        ne : array-like
            The plasma density at the specified index.
        z : double
            The z coordinate of the density.
        """
        ne = np.load(self.filePre + '_plasmaDensity_' + str(ind) + '.npy')
        z = self.z[ind]
        return ne, z
        
    # Visualization functions
    #--------------------------------------------------------------------------
    
    def plot_long_density_center(self):
        """ Plots the plasma density in an x-z plane at y=0. """
        Nz = self.Nz
        en = np.zeros((Nz, self.Nx))
        if not self.cyl:
            for i in range(1, Nz):
                en[i, :], z = self.load_plasma_density(i)[:, int(self.Ny/2)]
        else:
            for i in range(1, Nz):
                en[i, :], z = self.load_plasma_density(i)
        im = self.plot_long_density(en)
        plt.show(im)
    
    def plot_long_density(self, ne):
        """ Create a longitudinal plasma density plot. """
        Z = self.Z
        X = self.X
        
        ne = self.prep_data(ne)
        im = plt.imshow(ne, aspect='auto', extent=[0, Z, -X/2, X/2])
        cb = plt.colorbar()
        cb.set_label(r'Plasma density ($\mathrm{10^{17}cm^{-3}}$)')
        plt.set_cmap('plasma')
        plt.xlabel(r'z')
        plt.ylabel(r'x')
        plt.title('Longitudinal plasma density')
        return im


class UniformPlasma(Plasma):
    """ A uniform density gas that can be ionized into a plasma. """
    
    def __init__(self, params):
        super().__init__(params)
    
    def load_density(self, i):
        """ Returns a 2D array of the number density, always constant. """
        return np.full((self.Nx, self.Ny), self.n0, dtype='double')
