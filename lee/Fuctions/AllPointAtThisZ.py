#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 19:28:29 2021

@author: valentinalee
"""
#%%
import numpy as np

def AllPointAtThisZ (z, Zpoint, error= 0.1):
    '''
    Args:
        z: 2D array, z value, z(x, y)
        Zpoint: float, given point Zpoint, this fn will return points that are close to Zpoint
        error= 0.01: acceptance error, the function will return all the points where abs(z-Zpoint)<error
    
    Returns:
        x_acpt: 1D array, x value of the accepted points
        y_acpt: 1D array, y value of the accepted points
        z_acpt: 1D array, z value of the accepted points
    
    '''
    x_acpt= []
    y_acpt= []
    z_acpt= []
    for xidx in range (z.shape[1]):
        for yidx in range (z.shape[0]):
            if abs(z[yidx, xidx]-Zpoint)/Zpoint< error:
                x_acpt.append(xidx)
                y_acpt.append(yidx)
                z_acpt.append(z[yidx, xidx])
    return np.array(x_acpt), np.array(y_acpt), np.array(z_acpt)

