#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:08:32 2017

@author: robert
"""

class Element:
    """ The base class for elements.
    
    Implements base methods to check class construction.
    """
    keys = []
    
    def __init__(self, params):
        self.params = params
        self.check_params(params)
        self.params_to_attrs(params)
    
    def check_params(self, params):
        """ Check to ensure all required keys are in the params dictionary. """
        for key in self.keys:
            if key not in params:
                raise ValueError('The params object has no key: %s' % key)
    
    def params_to_attrs(self, params):
        """ Add all params as attributes of the class. """
        for key in params:
            setattr(self, key, params[key])
    