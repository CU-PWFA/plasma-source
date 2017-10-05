#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:09:17 2017

@author: robert
"""

import sys
# Add the python folder to the path
sys.path.insert(0, "../python")
import unittest
import numpy as np
from beam.beams import beam
from beam.beams import laserbeam
from beam.beams import laserpulse

#------------------------------------------------------------------------------
# beam.py
#------------------------------------------------------------------------------
class beam_class_test_cases(unittest.TestCase):
    """ Test cases for the beam class in beam/beams/beam.py """
    
    def setUp(self):
        params = {'name' : 'testBeam',
                  'path' : 'tests/'}
        self.beam = beam.Beam(params)
    
    def test_missing_check_params(self):
        """ Test if check_params catches a missing parameter """
        testBeam = self.beam
        testBeam.keys = ['test']
        params = {}
        self.assertRaises(ValueError, testBeam.check_params, params)
    
    def test_additional_check_params(self):
        """ Ensure check_params allows additional keys in params """
        params = {'name' : 'testBeam',
                  'path' : 'tests/',
                  'test' : 15}
        testBeam = beam.Beam(params)
        self.assertEqual(testBeam.test, params['test'])
    
    def test_params_to_attrs(self):
        """ Test if the parameter attributes are added as attributes """
        testBeam = self.beam
        params = {'name' : 'testBeam',
                  'path' : 'tests/',
                  'test1' : 12,
                  'test2' : 'string'}
        testBeam.params_to_attrs(params)
        self.assertEqual(testBeam.test1, params['test1'])
        self.assertEqual(testBeam.test2, params['test2'])
                
#------------------------------------------------------------------------------
# laserbeam.py
#------------------------------------------------------------------------------
class laser_class_test_cases(unittest.TestCase):
    """ Test cases for the laser class in beam/beams/laserbeam.py """
    
    def setUp(self):
        self.params = {'Nx' : 2**8,
                  'Ny' : 2**7,
                  'X' : 200,
                  'Y' : 100,
                  'lam' : 1.0,
                  'path' : 'tests/',
                  'name' : 'testBeam',
                  'threads' : 4,
                  'cyl' : False
                  }
        self.beam = laserbeam.Laser(self.params)
        
    def test_create_grid(self):
        """ Test if the grid was created correctly """
        testBeam = self.beam
        self.assertEqual(testBeam.x[0], -self.params['X']/2)
        self.assertEqual(testBeam.y[0], -self.params['Y']/2)
        self.assertEqual(testBeam.x[int(self.params['Nx']/2)], 0.0)
        self.assertEqual(testBeam.y[int(self.params['Ny']/2)], 0.0)
        
    def test_initialize_field_zero(self):
        """ Test that the field is initialized to 0 """
        for i in range(self.params['Nx']):
            for j in range(self.params['Ny']):
                self.assertEqual(self.beam.e[i, j], 0.0)
    
    def test_initialize_field_nonzero(self):
        """ Test initializing the field to some values """
        testBeam = self.beam
        e = np.ones((self.params['Nx'], self.params['Ny']), dtype='complex128')
        testBeam.initialize_field(e)
        for i in range(self.params['Nx']):
            for j in range(self.params['Ny']):
                self.assertEqual(self.beam.e[i, j], 1.0)
                
    def test_initialize_field_z(self):
        """ Test that z is set to 0 """
        self.beam.z = [5.0]
        self.beam.initialize_field()
        self.assertEqual(self.beam.z, [0.0])
        
    def test_initialize_field_saveInd(self):
        """ Test that the save index is reset- goes to 1 because of save """
        self.beam.saveInd = 5
        self.beam.initialize_field()
        self.assertEqual(self.beam.saveInd, 1)
    
    def test_set_field(self):
        """ Test if the field is set correctly and the z array updated """
        e = np.ones((self.params['Nx'], self.params['Ny']), dtype='complex128')
        self.beam.z = [0.0, 5.0]
        self.beam.set_field(e)
        for i in range(self.params['Nx']):
            for j in range(self.params['Ny']):
                self.assertEqual(self.beam.e[i, j], 1.0)
        self.assertEqual(self.beam.z, [0.0, 5.0, 5.0])
        
    def test_save_field_saveInd(self):
        """ Test that the saveInd is incremented """
        testBeam = self.beam
        testBeam.saveInd = 5
        testBeam.save_field(testBeam.e, 0.0)
        self.assertEqual(testBeam.saveInd, 6)
    
    def test_save_field_a(self):
        """ Test that the z array has a value appended to it """
        testBeam = self.beam
        testBeam.z = [0.0, 5.0]
        testBeam.save_field(testBeam.e, 8.0)
        self.assertEqual(testBeam.z, [0.0, 5.0, 8.0])

    def test_propagate_zeros(self):
        """ Test if propagating a zero field returns zeros """
        self.beam.propagate(100, 1.0)
        for i in range(self.params['Nx']):
            for j in range(self.params['Ny']):
                self.assertEqual(self.beam.e[i, j], 0.0)
                
#------------------------------------------------------------------------------
# laserpulse.py
#------------------------------------------------------------------------------
class pulse_class_test_cases(unittest.TestCase):
    """ Test cases for the laser class in beam/beams/laserbeam.py """
    
    def setUp(self):
        self.params = {'Nx' : 2**8,
                  'Ny' : 2**7,
                  'Nt' : 2**6,
                  'X' : 200,
                  'Y' : 100,
                  'T' : 100,
                  'lam' : 1.0,
                  'path' : 'tests/',
                  'name' : 'testBeam',
                  'threads' : 4,
                  'cyl' : False
                  }
        self.pulse = laserpulse.Pulse(self.params)  
        
    def test_create_grid(self):
        """ Test if the grid was created correctly """
        testPulse = self.pulse
        self.assertEqual(testPulse.t[0], -self.params['T']/2)
        self.assertEqual(testPulse.x[0], -self.params['X']/2)
        self.assertEqual(testPulse.y[0], -self.params['Y']/2)
        self.assertEqual(testPulse.t[int(self.params['Nt']/2)], 0.0)
        self.assertEqual(testPulse.x[int(self.params['Nx']/2)], 0.0)
        self.assertEqual(testPulse.y[int(self.params['Ny']/2)], 0.0)

#------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()