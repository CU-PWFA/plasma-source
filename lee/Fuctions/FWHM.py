#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 15:05:15 2021

@author: valentinalee
"""

import numpy as np

def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res

def FWHM(X,Y):
    '''
    Args:
        X: 1D array, x variable
        Y: 1D array, y= f(x)
    
    Returns:
        value: the FWHM of signal f(x)
    
    '''

    half_max = max(Y) / 2.
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    return X[right_idx] - X[left_idx] 