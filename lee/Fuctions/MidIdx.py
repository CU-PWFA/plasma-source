#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:18:19 2021

@author: valentinalee
"""

#%%

def MidIdx(array, axis=0):
    '''
    Args:
        array: n-D array that you want to find the middle index of the given axis
        axis: the axis you want to find the middle index
    
    Returns:
        Middle index of a given axis of a n-d array
        If the axis has an even len, return Mid-Right Index
    
    '''
    print(len(array.shape))
    assert (len(array.shape)> axis), 'Axis out of range!'
    if (array.shape[axis])%2==0:
        print('This axis has an even len. Return Mid-Right Index')
        return int((array.shape[axis])/2)
    else:
        return int(((array.shape[axis])-1)/2)

