#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 07:57:58 2021

@author: valentinalee
"""

#%%
import numpy as np

def polar_to_cartesian(radius_grid, angle_grid, data):
    '''
    Args:
        radius_grid: 1D array, radius grid of the polar data (Usually goes from 0 to 1)
        angle_grid: 1D array, angle grid of the polar data (Usually goes from 0 to 2pi)
        data: 2D array, x axis is the r value and y axis is the angle value
    
    Returns:
        NewArray: 2D array, in cartiesian coordinate 
    '''
    r_grid, a_grid = np.meshgrid(radius_grid, angle_grid)

    new= np.zeros_like(data)+0
    x= np.linspace(-1, 1, new.shape[1])
    y= np.linspace(-1, 1, new.shape[0])
    for i in range(new.shape[0]):
        for j in range(new.shape[1]):
            x0, y0= x[j], y[i]
            r, a= np.sqrt(x0**2 + y0**2), np.arctan2(y0, x0)
            data_i= np.argmin(np.abs(a_grid[:, 0] - a))
            data_j= np.argmin(np.abs(r_grid[0, :] - r))
            val= data[data_i, data_j]

            if r<=1:
                new[i, j]= val
            print(i, j)

    return new
