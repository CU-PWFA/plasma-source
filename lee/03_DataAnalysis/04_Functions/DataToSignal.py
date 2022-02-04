#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:53:10 2022

@author: valentinalee
"""
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy import optimize
from scipy.ndimage.interpolation import shift

def DataToSignal(avgBG, DataFile):
    '''
    avgBG: 2D array of the avg background
    Data: tiff file data
    --------
    Returns:
    Signal: 2D array of the signal (filtered)
    '''
    def shiftBG(s1, s2):
        ShiftedBG= shift(avgBG, (s1, s2))
        return ShiftedBG

    data= np.array(Image.open(DataFile))
    errorFn= lambda p1: np.ravel(shiftBG(*p1)- data)
    p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
    data= shift(data, (-p[0], -p[1]))
    signal= data-avgBG
    signalFiltered= gaussian_filter(signal, sigma=6)

    return signalFiltered
