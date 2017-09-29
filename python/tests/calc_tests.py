#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:59:22 2017

@author: robert
"""

import sys
# Add the python folder to the path
sys.path.insert(0, "../python")
import unittest
import numpy as np
from numpy.fft import fftfreq
import pyfftw
from beam.calc import laser

#------------------------------------------------------------------------------
# laser.pyx
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
        # Create a meaningless save function
        def save(e, z): return None
        self.save = save
    
    def test_gaussian_diffraction(self):
        """ Test if a gaussian beam diffracts properly, runs an array of z """
        # Create the Gaussian beam
        w0 = 10
        zR = np.pi * w0**2
        wZ = w0 * np.sqrt(1 + (self.Z/zR)**2)
        r2 = np.reshape(self.x-10, (self.Nx, 1))**2 +\
             np.reshape(self.y, (1, self.Ny))**2
        E = np.array(np.exp(- r2 / w0**2), dtype='complex128')
        # The paraxial answer
        EE = (w0 / wZ) * np.exp(-r2 / wZ**2)
        # Run the fourier propagator
        Etest = laser.fourier_prop(E, self.x, self.y, self.z, 1.0, 1.0,
                                   self.pfft, self.pifft, self.save)
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertAlmostEqual(abs(Etest[i, j]), EE[i, j], delta=1e-3)
                
    def test_zeros_input(self):
        """ If the function is passed all zeros it should return all zeros """
        E = np.zeros((self.Nx, self.Ny), dtype='complex128')
        z = np.array(1.0, ndmin=1)
        # Run the fourier propagator
        Etest = laser.fourier_prop(E, self.x, self.y, z, 1.0, 1.0,
                                   self.pfft, self.pifft, self.save)
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertEqual(Etest[i, j], 0.0)
                
    def test_save_function_call(self):
        """ Test that the save function is called correctly """
        E = np.zeros((self.Nx, self.Ny), dtype='complex128')
        z = np.array(1.5, ndmin=1)
        called = False
        def save(e, z): 
            nonlocal called
            self.assertEqual(z, 1.5)
            called = True
        laser.fourier_prop(E, self.x, self.y, z, 1.0, 1.0, self.pfft,
                           self.pifft, save)
        self.assertTrue(called)

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
        fz2 = 1.0 - np.reshape(fx2, (self.Nx, 1))-np.reshape(fy2, (1, self.Ny))
        fz2 = np.array(fz2, dtype='complex128')
        self.ikz = pre*np.sqrt(fz2)
        
    def test_gaussian_diffraction(self):
        """ Test if a gaussian beam diffracts properly """
        # Create the Gaussian beam
        Z = 100
        w0 = 10
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
    
    def test_ikz_values(self):
        """ Check an entire array of ikz values """
        fx2 = self.fx**2
        fy2 = self.fy**2
        pre = 1j*2*np.pi
        fz2 = 1.0 - np.reshape(fx2, (self.Nx, 1))-np.reshape(fy2, (1, self.Ny))
        fz2 = np.array(fz2, dtype='complex128')
        ikz = pre*np.sqrt(fz2)
        ikzTest = laser.ikz_RS(self.fx, self.fy, 1.0, 1.0)
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.assertEqual(ikzTest[i, j], ikz[i, j])
    
    def test_complex_frequency(self):
        """ Test to ensure it returns a complex frequency """
        fx = np.array([1.0])
        fy = np.array([1.0])
        ikz = laser.ikz_RS(fx, fy, 1.0, 1.0)
        self.assertEqual(ikz[0, 0], -2*np.pi)        

#------------------------------------------------------------------------------
class fresnel_axis_test_cases(unittest.TestCase):
    """ Test cases for fresnel_axis in beam/calc/laser.py """
    
    def setUp(self):
        # Create the grid, both r and z start at 0
        self.Nr = 2**8
        R = 100
        self.r = np.linspace(0, R, self.Nr)
        self.Nz = 2**8
        z0 = 500
        Z = 100
        self.z = np.linspace(z0, z0+Z, self.Nz)
    
    def test_gaussian_propagation(self):
        """ Check that the power of a Gaussian decreases correctly """
        w0 = 10
        zR = np.pi * w0**2
        wZ = w0 * np.sqrt(1 + (self.z/zR)**2)
        E = np.array(np.exp(- self.r**2 / w0**2), dtype='complex128')
        EE = (w0 / wZ)
        Etest = laser.fresnel_axis(E, self.r, self.z, 1.0, 1.0)
        for i in range(self.Nz):
            self.assertAlmostEqual(abs(Etest[i]), EE[i], delta=1e-3)
    
    def test_zero_z(self):
        """ Check that an exception is raised for z=0 """
        z = np.linspace(0, 100, self.Nz)
        w0 = 10
        E = np.array(np.exp(- self.r**2 / w0**2), dtype='complex128')
        self.assertRaises(ValueError, laser.fresnel_axis, E, self.r, z, 1.0,
                          1.0)

#------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()
