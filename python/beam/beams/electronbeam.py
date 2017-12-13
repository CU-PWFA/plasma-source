#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:49:43 2017

@author: robert
"""

import numpy as np
from beam.beams import beam
import matplotlib.pyplot as plt

class ElectronBeam(beam.Beam):
    """ An electron beam class that stores a collection of macroparticles.
    
    Parameters
    ----------
    N : int
        The number of macroparticles in the beam.
    """
    keys = ['N']
    
    # Initialization functions
    #--------------------------------------------------------------------------
    
    def __init__(self, params):
        super().__init__(params)
        # Create internal variables
        if self.load is False:
            self.initialize_particles()
            self.save_initial()
    
    def initialize_particles(self, ptcls=None):
        """ Create the arrays to store particle positions and momentum in. 
        
        Parameters
        ----------
        ptcls : array-like, optional
            The array of particle position momentum to initialize.
        """
        if ptcls is None:
            self.ptcls = np.zeros((self.N, 6), dtype='double')
        else:
            self.ptcls = np.array(ptcls, dtype='double')
        self.saveInd = 0
        self.z = []
        self.save_ptcls(self.ptcls, 0.0)
    
    def load_beam(self):
        """ Load the beam, specfically load the z-grid and saveInd. """
        self.z = np.load(self.filePre + '_z.npy')
        self.saveInd = len(self.z)
        self.ptcls = self.load_ptcls(self.saveInd - 1)[0]
    
    # Getters and setters
    #--------------------------------------------------------------------------
    
    def get_x(self, ptcls):
        return ptcls[:, 0]
    
    def get_xp(self, ptcls):
        return ptcls[:, 1]
    
    def get_y(self, ptcls):
        return ptcls[:, 2]
    
    def get_yp(self, ptcls):
        return ptcls[:, 3]
    
    def get_z(self, ptcls):
        return ptcls[:, 4]
    
    def get_gamma(self, ptcls):
        return ptcls[:, 5]
    
    #File managment
    #--------------------------------------------------------------------------
    
    def save_ptcls(self, ptcls, z):
        """ Save the current particle distribution to file and advance z. """
        np.save(self.filePre + '_ptcls_' + str(self.saveInd) + '.npy', ptcls)
        self.saveInd += 1
        self.z.append(z)
        self.save_z()
    
    def save_z(self):
        """ Save the z array. """
        np.save(self.filePre + '_z.npy', self.z)
    
    def load_field(self, ind):
        """ Load the particle distribution at the specified index. 
        
        Parameters
        ----------
        ind : int
            The save index to load the field at.
        
        Returns
        -------
        ptcls : array-like
            The particle corrdinates and momenta at the specified index.
        z : double
            The z coordinate of the field.
        """
        e = np.load(self.filePre + '_ptcls_' + str(ind) + '.npy')
        z = self.z[ind]
        return e, z
    
    # Physics functions
    #--------------------------------------------------------------------------
        
    def propagate(self, z, n):
        """ Propagate the field to an array of z distances. """
        #TODO implement this function
        
    # Visualization functions
    #--------------------------------------------------------------------------
    
    def plot_current_phase(self):
        """ Plots a scatterplot of the particles in the current beam. """
        self.plot_phase(self.ptcls, self.z[-1])
        plt.show()
    
    def plt_phase_at(self, ind):
        """ Plots the particles at a particular z distance.
        
        Parameters
        ----------
        ind : int
            The save index to plot the particles at, see the _z file to find z.
        """
        ptcls, z = self.load_field(ind)
        self.plot_phase(ptcls, z)
        plt.show()
    
    def plot_phase(self, ptcls, z):
        """ Create an x-y plot of the particles. """        
        fig = plt.figure(figsize=(10, 4))
        plt.subplot(121)
        plt.scatter(ptcls[:, 0]*1e6, ptcls[:, 1]*1e3, 1)
        plt.title('x phase space')
        plt.xlabel('x (um)')
        plt.ylabel("x' (mrad)")
        plt.subplot(122)
        plt.scatter(ptcls[:, 2]*1e6, ptcls[:, 3]*1e3, 1)
        plt.title('y phase space')
        plt.xlabel('y (um)')
        plt.ylabel("y' (mrad)")

        plt.tight_layout()
        return fig


class GaussianElectronBeam(ElectronBeam):
    """ A electron beam with a Gaussian transverse profile. 
    
    Parameters
    ----------
    gamma : double
        The relativistic factor of the beam.
    emittance : double
        The normalized emittance of the beam in m*rad.
    betax : double
        The beta function in the x-direction at z=0, in m.
    betay : double
        The beta function in the y direction at z=0, in m.
    alphax : double
        The alpha parameter of the beam in the x-direction at z=0.
    alphay : double
        The alpha parameter of the beam in the y-direction at z=0.
    sigmaz : double
        The RMS size of the beam in z in m.
    dE : double
        The energy spread of the beam as a fraction, 0.01 = +-1% energy spread.
    """
    
    def __init__(self, params):
        self.keys.extend(
                ['gamma',
                 'emittance',
                 'betax',
                 'betay',
                 'alphax',
                 'alphay',
                 'sigmaz',
                 'dE'])
        super().__init__(params)
    
    def initialize_particles(self):
        """ Initialize the particles in a 6D distribution. """
        N = self.N
        gamma = self.gamma
        emittance = self.emittance
        betax = self.betax
        betay = self.betay
        ptcls = np.zeros((N, 6), dtype='double')
        # Calculate arrays of random numbers
        x1r = np.random.uniform(0, 1, N)
        x2r = np.random.uniform(0, 1, N)
        y1r = np.random.uniform(0, 1, N)
        y2r = np.random.uniform(0, 1, N)
        # Choose the particles from a distribution
        Jx = -emittance * np.log(x1r) / gamma
        Jy = -emittance * np.log(y1r) / gamma
        phix = 2*np.pi*x2r
        phiy = 2*np.pi*y2r
        ux = np.sqrt(2*Jx)*np.cos(phix)
        vx = -np.sqrt(2*Jx)*np.sin(phix)
        uy = np.sqrt(2*Jy)*np.cos(phiy)
        vy = -np.sqrt(2*Jy)*np.sin(phiy)
        # Calculate the coordinates
        ptcls[:, 0] = ux*np.sqrt(betax)
        ptcls[:, 1] = (vx-self.alphax*ux) / np.sqrt(betax)
        ptcls[:, 2] = uy*np.sqrt(betay)
        ptcls[:, 3] = (vy-self.alphay*uy) / np.sqrt(betay)
        ptcls[:, 4] = np.random.normal(0.0, self.sigmaz, N)
        ptcls[:, 5] = gamma * (1 + self.dE*np.random.uniform(-1, 1, N))
        super().initialize_particles(ptcls) 
