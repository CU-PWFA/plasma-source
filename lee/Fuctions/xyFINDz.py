#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 21:29:44 2021

@author: valentinalee
"""

#%%
import numpy as np


def xyFINDz (Xaxis, Yaxis, CorrZ, GivenX, GivenY):
    '''
    Args:
        Xaxis: 1D array, the x axis of the 2D array
        Yaxis: 1D array, the y axis of the 2D array
        CorrZ: 2D array, z values
        GivenX: float, the given x where you want to find its z(x, y)
        GivenY: float, the given y where you want to find its z(x, y)
    
    Returns:
        ZValue: float, z(GivenX, GivenY)
    
    '''
    idxX= (np.abs(Xaxis-GivenX)).argmin()
    idxY= (np.abs(Yaxis-GivenY)).argmin()
    return(CorrZ[idxY, idxX])

    
