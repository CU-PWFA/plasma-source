#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 14:56:58 2017

@author: litos
"""
 
## common libraries
from numpy import *

## define constants
c = 3e8 # m/s, speed of light

## define function
def calcwaist(T):
    [beta,alpha,gamma] = T
    beta0 = beta/(1+(alpha**2))
    z0 = -alpha*beta0
    
    return [beta0,z0]