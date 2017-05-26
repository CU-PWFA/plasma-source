#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 21:08:32 2017

@author: litos
"""

import numpy as np

def calc_mean(x,frac=1.0):
        
    med   = np.median(x)
    dmed  = x-med
    idmed = sorted(range(len(dmed)),key=lambda idmed: dmed[idmed])
    nx = round(frac*len(dmed))
    x = x[:nx]
    mu = np.mean(x)
    dmu = x-mu
    rms = np.sqrt(np.sum(dmu**2)/(len(dmu)-1))
     
    return rms