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
from beam.beams import beam
from beam.beams import laserbeam

#------------------------------------------------------------------------------
# beam.py
#------------------------------------------------------------------------------
class beam_class_test_cases(unittest.TestCase):
    """ Test cases for the beam class in beam/beams/beam.py """
    
    def setUp(self):
        params = {}
        self.beam = beam.Beam(params)
    
    def test_missing_check_params(self):
        """ Test if check_params catches a missing parameter """
        testBeam = self.beam
        testBeam.keys = ['test']
        params = {}
        self.assertRaises(ValueError, testBeam.check_params, params)
    
    def test_additional_check_params(self):
        """ Ensure check_params allows additional keys in params """
        params = {'test' : 15}
        testBeam = beam.Beam(params)
        self.assertEqual(testBeam.test, params['test'])
    
    def test_params_to_attrs(self):
        """ Test if the parameter attributes are added as attributes """
        testBeam = self.beam
        params = {'test1' : 12,
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
                  'threads' : 4
                  }
        self.beam = laserbeam.Laser(self.params)
        
    def test_create_grid(self):
        """ Test if the grid was created correctly """
        testBeam = self.beam
        self.assertEqual(testBeam.x[0], -self.params['X']/2)
        self.assertEqual(testBeam.y[0], -self.params['Y']/2)
        self.assertEqual(testBeam.x[int(self.params['Nx']/2)], 0.0)
        self.assertEqual(testBeam.y[int(self.params['Ny']/2)], 0.0)

    def test_propagate_zeros(self):
        """ Test if propagating a zero field returns zeros """
        self.beam.propagate(100, 1.0)
        for i in range(self.params['Nx']):
            for j in range(self.params['Ny']):
                self.assertEqual(self.beam.e[i, j], 0.0)

#------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()