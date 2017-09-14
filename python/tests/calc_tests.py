#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:59:22 2017

@author: robert
"""

import sys
# Add the python folder to the path
sys.path.insert(0, "../../python")
import unittest
import numpy as np
from numpy.fft import fftfreq
import pyfftw
import beam.calc.laser as laser


#------------------------------------------------------------------------------
class fourier_prop_test_cases(unittest.TestCase):
    """ Test cases for fourier_prop in beam/calc/laser.py """
    
    def setUp(self):
        # Create the grid
        self.Nx = 2**8
        X = 100
        self.x = np.linspace(-X/2, X/2, self.Nx, False)
        self.Ny = 2**8
        Y = 100
        self.y = np.linspace(-Y/2, Y/2, self.Ny, False)
        
        self.Nz = 2**8
        self.Z = 100
        self.z = np.linspace(0, self.Z, self.Nz)
        # Setup fftw plans
        e = pyfftw.empty_aligned((self.Nx, self.Ny), dtype='complex128')
        self.pfft = pyfftw.builders.fft2(e, overwrite_input=True,
                                         avoid_copy=True, threads=4)
        self.pifft = pyfftw.builders.ifft2(e, overwrite_input=True, 
                                           avoid_copy=True, threads=4)
    
    def test_gaussian_diffraction(self):
        """ Test if a gaussian beam diffracts properly, runs an array of z """
        # Create the Gaussian beam
        w0 = 5
        zR = np.pi * w0**2
        wZ = w0 * np.sqrt(1 + (self.Z/zR)**2)
        r2 = np.reshape(self.x-10, (self.Nx, 1))**2 +\
             np.reshape(self.y, (1, self.Ny))**2
        E = np.array(np.exp(- r2 / w0**2), dtype='complex128')
        # The paraxial answer
        EE = (w0 / wZ) * np.exp(-r2 / wZ**2)
        # Run the fourier propagator
        Etest = laser.fourier_prop(E, self.x, self.y, self.z, 1.0, 1.0,
                                   self.pfft, self.pifft, 'test/')
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertAlmostEqual(abs(Etest[i, j]), EE[i, j], delta=1e-3)
                
    def test_zeros_input(self):
        """ If the function is passed all zeros it should return all zeros """
        E = np.zeros((self.Nx, self.Ny), dtype='complex128')
        # Run the fourier propagator
        Etest = laser.fourier_prop(E, self.x, self.y, 1.0, 1.0, 1.0,
                                   self.pfft, self.pifft, 'test/')
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertEqual(Etest[i, j], 0.0)

#------------------------------------------------------------------------------
class fourier_step_test_cases(unittest.TestCase):
    """ Test cases for fourier_step in beam/calc/laser.py """
    
    def setUp(self):
        # Create the grid
        self.Nx = 2**8
        X = 100
        self.x = np.linspace(-X/2, X/2, self.Nx, False)
        self.Ny = 2**8
        Y = 100
        self.y = np.linspace(-Y/2, Y/2, self.Ny, False)
        
        # Setup fftw plans
        e = pyfftw.empty_aligned((self.Nx, self.Ny), dtype='complex128')
        self.pfft = pyfftw.builders.fft2(e, overwrite_input=True,
                                         avoid_copy=True, threads=4)
        self.pifft = pyfftw.builders.ifft2(e, overwrite_input=True, 
                                           avoid_copy=True, threads=4)
        # Calculate the RS spatial frequencies
        dx = X / (self.Nx-1)
        dy = Y / (self.Ny-1)
        fx2 = fftfreq(self.Nx, dx)**2
        fy2 = fftfreq(self.Ny, dy)**2
        pre = 1j*2*np.pi
        self.ikz = pre*np.sqrt(1.0 - np.reshape(fx2, (self.Nx, 1)) -\
                             np.reshape(fy2, (1, self.Ny)))
        
    def test_gaussian_diffraction(self):
        """ Test if a gaussian beam diffracts properly"""
        # Create the Gaussian beam
        Z = 100
        w0 = 5
        zR = np.pi * w0**2
        wZ = w0 * np.sqrt(1 + (Z/zR)**2)
        r2 = np.reshape(self.x-10, (self.Nx, 1))**2 +\
             np.reshape(self.y, (1, self.Ny))**2
        E = np.array(np.exp(- r2 / w0**2), dtype='complex128')
        # The paraxial answer
        EE = (w0 / wZ) * np.exp(-r2 / wZ**2)
        # Run the fourier propagator
        Etest = laser.fourier_step(E, self.ikz, Z, self.pfft, self.pifft)
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertAlmostEqual(abs(Etest[i, j]), EE[i, j], delta=1e-3)
    
    def test_zeros_input(self):
        """ If the function is passed all zeros it should return all zeros """
        Z = 100
        E = np.zeros((self.Nx, self.Ny), dtype='complex128')
        # Run the fourier propagator
        Etest = laser.fourier_step(E, self.ikz, Z, self.pfft, self.pifft)
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertEqual(Etest[i, j], 0.0)

#------------------------------------------------------------------------------
class ikz_RS_test_cases(unittest.TestCase):
    """ Test cases for ikz_RS in beam/calc/laser.py """
    
    def setUp(self):
        # Create the grid
        self.Nx = 2**8
        X = 100
        self.Ny = 2**8
        Y = 100
        
        # Calculate the RS spatial frequencies
        dx = X / (self.Nx-1)
        dy = Y / (self.Ny-1)
        self.fx = fftfreq(self.Nx, dx)
        self.fy = fftfreq(self.Ny, dy)
        fx2 = self.fx**2
        fy2 = self.fy**2
        pre = 1j*2*np.pi
        self.ikz = pre*np.sqrt(1.0 - np.reshape(fx2, (self.Nx, 1)) -\
                             np.reshape(fy2, (1, self.Ny)))
    
    def test_ikz_values(self):
        ikz = laser.ikz_RS(self.fx, self.fy, 1.0, 1.0)
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertEqual(ikz[i, j], self.ikz[i, j])

#------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()
