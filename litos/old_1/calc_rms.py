#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:44:57 2017

@author: litos
"""

import numpy as np
from calc_mean import calc_mean

def calc_rms(x,frac=1.0):
        
    mu  = calc_mean(x,frac)
    dx = x-mu
    dx.sort()
    ndx = round(frac*len(dx))
    dx = dx[:ndx]
    rms = np.sqrt(np.sum(dx**2)/(len(dx)-1))
     
    return rms