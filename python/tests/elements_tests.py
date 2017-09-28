#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:10:47 2017

@author: robert
"""

import sys
# Add the python folder to the path
sys.path.insert(0, "../python")
import unittest
import numpy as np
from beam.elements import element

#------------------------------------------------------------------------------
# element.py
#------------------------------------------------------------------------------
class element_class_test_cases(unittest.TestCase):
    """ Test cases for the beam class in beam/elements/element.py """
    
    def setUp(self):
        params = {}
        self.elem = element.Element(params)
    
    def test_missing_check_params(self):
        """ Test if check_params catches a missing parameter """
        testElem = self.elem
        testElem.keys = ['test']
        params = {}
        self.assertRaises(ValueError, testElem.check_params, params)
    
    def test_additional_check_params(self):
        """ Ensure check_params allows additional keys in params """
        params = {'test' : 15}
        testElem = element.Element(params)
        self.assertEqual(testElem.test, params['test'])
    
    def test_params_to_attrs(self):
        """ Test if the parameter attributes are added as attributes """
        testElem = self.elem
        params = {'test1' : 12,
                  'test2' : 'string'}
        testElem.params_to_attrs(params)
        self.assertEqual(testElem.test1, params['test1'])
        self.assertEqual(testElem.test2, params['test2'])

#------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()