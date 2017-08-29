#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:47:08 2017

@author: chris
"""

import numpy as np

def CalcFWHM(axis, vals):
    return CalcTopWidth(axis, vals, 0.50)
    
def CalcTopWidth(axis, vals, thresh = 0.95):
    axis = np.array(axis)
    vals = np.array(vals)
    
    top = max(vals) * thresh
    left = 0;  right = 0
    for i in range(len(vals)):
        if vals[i] >= top:
            right = i
            if left == 0:
                left = i
    return axis[right] - axis[left]

def CalcUpRamp(axis, vals, start = 0.10, stop = 0.95):
    axis = np.array(axis)
    vals = np.array(vals)
    
    bot = max(vals) * start
    top = max(vals) * stop
    lower = 0;  upper = 0
    for i in range(len(vals)):
        if ((vals[i] >= bot) & (lower == 0)):
                lower = i
        if ((vals[i] >= top) & (upper == 0)):
                upper = i
    return axis[upper] - axis[lower]

def CalcDownRamp(axis, vals, start = 0.10, stop = 0.95):
    return -1 * CalcUpRamp(axis[::-1],vals[::-1])