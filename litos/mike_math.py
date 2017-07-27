#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:49:33 2017

@author: mike
"""

import numpy as np

def calc_mean(x,frac=1.0):
    if frac<1.0:
        med   = np.median(x)
        dmed  = x-med
        idmed = sorted(range(len(dmed)),key=lambda idmed: dmed[idmed])
        nx = round(frac*len(idmed))
        idmed = idmed[:nx]
        mu = np.mean(x[idmed])
    else:
        mu = np.mean(x)
    return mu

def calc_rms(x,frac=1.0):
    mu  = calc_mean(x,frac)
    dx = x-mu
    if frac<1.0:
        dx.sort()
        ndx = round(frac*len(dx))
        dx = dx[:ndx]
    rms = np.sqrt(np.sum(dx**2)/(len(dx)-1))
    return rms