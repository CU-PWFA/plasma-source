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
        Gas number density for a uniform plasma in 10^17 cm^-3.
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
    
    def save_density(self, n, ind):
        """ Save the plasma density to file at the given z ind. """
        if self.cyl:
            n = n[:, int(self.Ny/2)]
        np.save(self.filePre + '_plasmaDensity_' + str(ind) + '.npy', n)
    
    # Visualization functions
    #--------------------------------------------------------------------------
    
    
