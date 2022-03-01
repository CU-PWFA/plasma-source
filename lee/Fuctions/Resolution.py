#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:19:18 2021

@author: valentinalee
"""

import numpy as np
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2

def Resolution(Array):
    '''
    Args:
        Array: 2D array that you test its resolution
    
    Returns:
        value: resolution with arbitrary unit (Use to compair two or more imagies)
    
    '''
    FT= np.real(fft2(Array))
    Reso= 0
    for yidx in range (FT.shape[0]):
        for xidx in range (FT.shape[1]):
            Reso= Reso+ FT[yidx, xidx]*np.sqrt(xGrid[xidx]**2+yGrid[yidx]**2)
    return Reso