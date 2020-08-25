# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 16:50:02 2020

@author: Robert
"""

import numpy as np

def load_flow_data(file, ncol):
    """ Load fluid data from a flow sim output text file.
    
    Parameters
    ----------
    file : string
        Path/filename of the solidworks output file.
    ncol : int
        Number of columns in the exported file.
        
    Returns
    -------
    data : array of floats
        Data from the file with non fluid cells set to nan.
    """
    converters = {}
    for i in range(3, ncol):
        converters[i] = convert_element
    data = np.loadtxt(file, delimiter='\t', skiprows=1, converters=converters)
    print('Data size:', np.shape(data))
    return data

def convert_element(s):
    """ Converter for np.loadtxt that sets empty elements to nan. """
    if s==b'':
        return 0.0
    else:
        return s