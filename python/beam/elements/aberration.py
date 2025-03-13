#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:28:09 2022

@author: valentinalee
"""


import numpy as np
from beam.elements import element
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import interpolate
from scipy.ndimage import gaussian_filter

class Phase(element.Element):
    """ A class that represents some type of phase mask on a laser beam.
    
    This class is meant to serve as a base class for more complex phase masks
    such as lens, diffraction gratings, etc. Note, phase is in radians.
    
    Parameters
    ----------
    Nx : int
        Number of grid points in the x direction, must match beam.
    Ny : int
        Number of grid points in the y direction, must match beam.
    X : int
        Width of the grid in the x direction, the grid goes from [-X/2, X/2).
        This must match the beam parameters.
    Y : int
        Width of the grid in the y direction, the grid goes from [-Y/2, Y/2).
        This must match the beam parameters.
    path : string
        The path for the calculation. This class will create a folder inside
        the path to store all output data in.
    name : string
        The name of the beam, used for naming files and folders.
    """
    keys = ['Nx',
            'Ny',
            'X',
            'Y',
            'path',
            'name',
            'lam']
    
    # Initialization functions
    #--------------------------------------------------------------------------
    
    def __init__(self, params):
        super().__init__(params)
        self.k = 2*np.pi/self.lam
        if self.load is False:
            self.create_grid()
            self.initialize_phase()
            self.save_initial()        
    
    def create_grid(self):
        """ Create an x-y rectangular grid. """
        X = self.X
        Y = self.Y
        self.x = np.linspace(-X/2, X/2, self.Nx, False, dtype='double')
        self.y = np.linspace(-Y/2, Y/2, self.Ny, False, dtype='double')
        
    def initialize_phase(self, phi=None):
        """ Create the array to store the phase in. 
        
        Parameters
        ----------
        phi : array-like, optional
            The array of phase to initialize the mask.
        """
        if phi is None:
            self.phi = np.zeros((self.Nx, self.Ny), dtype='complex128')
        else:
            self.phi = phi
        self.save_phase()
        
    def load_element(self):
        """ Load the phase mask. """
        self.create_grid()
        self.load_phase()
    
    #File managment
    #--------------------------------------------------------------------------
    
    def save_initial(self):
        """ Save the initial params object and the grid. """
        super().save_initial()
        np.save(self.filePre + '_x.npy', self.x)
        np.save(self.filePre + '_y.npy', self.y)
    
    def save_phase(self):
        """ Save the phase mask to file. """
        np.save(self.filePre + '_phase.npy', self.phi)
        
    def load_phase(self):
        """ Load the phase of the mask. """
        self.phi = np.load(self.filePre + '_phase.npy')
        
        
class Intensity(element.Element):
    """ A class that represents some type of transmission mask on a laser beam.
    
    This class is meant to serve as a base class for more complex transmission masks
    such as such as apertures.
    
    Parameters
    ----------
    Nx : int
        Number of grid points in the x direction, must match beam.
    Ny : int
        Number of grid points in the y direction, must match beam.
    X : int
        Width of the grid in the x direction, the grid goes from [-X/2, X/2).
        This must match the beam parameters.
    Y : int
        Width of the grid in the y direction, the grid goes from [-Y/2, Y/2).
        This must match the beam parameters.
    path : string
        The path for the calculation. This class will create a folder inside
        the path to store all output data in.
    name : string
        The name of the beam, used for naming files and folders.
    """
    keys = ['Nx',
            'Ny',
            'X',
            'Y',
            'path',
            'name',
            'lam']
    
    # Initialization functions
    #--------------------------------------------------------------------------
    
    def __init__(self, params):
        super().__init__(params)
        self.create_grid()
        self.initialize_t()
        self.save_initial()
    
    def create_grid(self):
        """ Create an x-y rectangular grid. """
        X = self.X
        Y = self.Y
        self.x = np.linspace(-X/2, X/2, self.Nx, False, dtype='double')
        self.y = np.linspace(-Y/2, Y/2, self.Ny, False, dtype='double')
        
    def initialize_t(self, t=None):
        """ Create the array to store the mask in. 
        
        Parameters
        ----------
        t : array-like, optional
            The array of transmission fraction to initialize the mask.
        """
        if t is None:
            self.t = np.zeros((self.Nx, self.Ny), dtype='double')
        else:
            self.t = t
        self.save_t()
    
    #File managment
    #--------------------------------------------------------------------------
    
    def save_initial(self):
        """ Save the initial params object and the grid. """
        super().save_initial()
        np.save(self.filePre + '_x.npy', self.x)
        np.save(self.filePre + '_y.npy', self.y)
    
    def save_t(self):
        """ Save the transmission mask to file. """
        np.save(self.filePre + '_t.npy', self.t)


class Noise(Phase):
    """ A phase mask that simulates a noise phase. Each pixel has a random value. 
    
    Parameters
    ----------
    RMS : double
        RMA of the random phase in length unit. 
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['RMS'])
        super().__init__(params)
    
    def initialize_phase(self):
        phi = self.k* np.random.normal(scale= self.RMS, size= (self.Nx, self.Ny))
        super().initialize_phase(phi)

class Random(Phase):
    """ A phase mask that simulates a smooth random phase variation.
    
    Parameters
    ----------
    RMS : double
        RMA of the random phase in length unit.
    scale : double
        Characteristic length scale of a "chuck" of random phase in length unit. 
        (e.g. if the length scale is equal to the pixel size, that is equivalent to "Noise")
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['RMS', 
                 'scale'])
        super().__init__(params)
    
    def initialize_phase(self):
        xb= np.arange(-self.X/2, self.X/2, self.scale)
        yb= np.arange(-self.Y/2, self.Y/2, self.scale)
        bigPhase= np.random.normal(scale= self.RMS, size= (len(xb), len(yb)))
        Smooth= gaussian_filter(bigPhase, sigma=2)
        f= interpolate.interp2d(xb, yb, Smooth, kind='cubic')
        phi = self.k* f(self.x, self.y)
        super().initialize_phase(phi)

class Defocus(Phase):
    """ A phase mask that simulates a defocus phase.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r= np.sqrt(self.x[:, None]**2+self.y[None, :]**2)
        rho= r/self.R
        phi= self.k* self.PV/2 *(2*rho**2-1)
        super().initialize_phase(phi)
        
class Astigmatism0(Phase):
    """ A phase mask that simulates a 0 degree astigmatism.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r2= self.x[:, None]**2+self.y[None, :]**2
        theta= np.arctan2(self.y[None, :], self.x[:, None])
        phi= self.k* self.PV/2 *r2/(self.R**2)*np.cos(2* theta)
        super().initialize_phase(phi)

class Astigmatism45(Phase):
    """ A phase mask that simulates a 45 degree astigmatism.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r2= self.x[:, None]**2+self.y[None, :]**2
        theta= np.arctan2(self.y[None, :], self.x[:, None])
        phi= self.k* self.PV/2 *r2/(self.R**2)*np.sin(2* theta)
        super().initialize_phase(phi)

class Coma0(Phase):
    """ A phase mask that simulates a 0 degree coma.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r= np.sqrt(self.x[:, None]**2+self.y[None, :]**2)        
        theta= np.arctan2(self.y[None, :], self.x[:, None])
        rho= r/self.R
        phi= self.k* self.PV/2 *(3*rho**3-2*rho)*np.sin(theta)
        super().initialize_phase(phi)

class Coma90(Phase):
    """ A phase mask that simulates a 90 degree coma.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r= np.sqrt(self.x[:, None]**2+self.y[None, :]**2)        
        theta= np.arctan2(self.y[None, :], self.x[:, None])
        rho= r/self.R
        phi= self.k* self.PV/2 *(3*rho**3-2*rho)*np.cos(theta)
        super().initialize_phase(phi)

class Spherical3rd(Phase):
    """ A phase mask that simulates a thrid spherical phase.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r= np.sqrt(self.x[:, None]**2+self.y[None, :]**2)
        rho= r/self.R
        phi= self.k* self.PV/2 *(6*rho**4-6*rho**2+1)
        super().initialize_phase(phi)
        
class Trefoil0(Phase):
    """ A phase mask that simulates a 0 degree trefoil.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r= np.sqrt(self.x[:, None]**2+self.y[None, :]**2)
        theta= np.arctan2(self.y[None, :], self.x[:, None])
        rho= r/self.R
        phi= self.k* self.PV/2 *(rho**3)*np.cos(3*theta)
        super().initialize_phase(phi)
        
class Trefoil90(Phase):
    """ A phase mask that simulates a 90 degree trefoil.
    
    Parameters
    ----------
    PV : double
        PV of the Zernike polynomial in length unit. 
    R : double
        Radias of the beam.
    
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['PV', 
                 'R'])
        super().__init__(params)
    
    def initialize_phase(self):
        r= np.sqrt(self.x[:, None]**2+self.y[None, :]**2)
        theta= np.arctan2(self.y[None, :], self.x[:, None])
        rho= r/self.R
        phi= self.k* self.PV/2 *(rho**3)*np.sin(3*theta)
        super().initialize_phase(phi)