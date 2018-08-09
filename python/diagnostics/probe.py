#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 09:28:25 2018

@author: robert
"""

import sys
sys.path.insert(0, "../../../photoDAQ")

import numpy as np
import file

def get_shots(dataSets):
    """ Load the number of shots in each data set.
    
     Parameters
    ----------
    dataSets : array-like
        The array of data set numbers.
        
    Returns
    -------
    shots : array-like
        The number of shots in each data set.
    """
    meta = {}
    M = len(dataSets)
    shots = np.zeros(M, dtype='int')
    for i in range(M):
        meta[i] = file.load_meta(dataSets[i])
        shots[i] = meta[i]['Shots']
    return shots


def load_data(dataSets, sign, shots):
    """ Load the scope traces from the data set files. 
    
    Parameters
    ----------
    dataSets : array-like
        The array of data set numbers.
    sign : array-like
        Array of the voltage sign for each data set, 1 positive, 0 negative.
    shots : array-like
        The number of shots in each data set.
        
    Returns
    -------
    data : array-like
        Data array where each trace is a row.
    volt : array-like
        Power supply voltage array, each voltage corresponds to a trace.
    time : array-like
        The time array for a single trace.
    """
    N = np.sum(shots) # Number of traces taken
    M = len(dataSets)
    data = np.zeros((N, 2500), dtype='double')
    volt = np.zeros(N, dtype='double')
    # All the traces should have the same time scale
    time = np.array(file.load_TRACE('TDS2024C', dataSets[0], 1)['t'],
                    dtype='double')
    ind = 0
    for i in range(M):
        for j in range(shots[i]):
            data[ind, :] = file.load_TRACE('TDS2024C', dataSets[i], j+1)['y']
            volt[ind] = file.load_SET('KA3005P', dataSets[i], j+1)['voltage']
            if sign[i] == 0:
                volt[ind] *= -1
            ind += 1
    return data, volt, time


def shift_data(data, trigger, sign):
    """ Adjust the data to account for different trigger levels. 
    
    Shifts the data in the array so the signal starts at the first element.
    
    Parameters
    ----------
    data : array-like
        Data array where each row is a trace.
    trigger : array-like
        The trigger level to use for finding the start.
    sign : int
        Trigger on 0 rising, 1 falling, or 2 either.
        
    Returns
    -------
    dataShifted : array-like
        The array of shifted data.
    """
    shape = np.shape(data)
    dataShifted = np.zeros(shape, dtype='double')
    for i in range(shape[0]):
        for j in range(2500):
            if data[i, j] > trigger and sign == 0:
                ind = j
                break
            if data[i, j] < -trigger and sign == 1:
                ind = j
                break
            if abs(data[i, j]) > trigger and sign == 2:
                ind = j
                break
        dataShifted[i, :] = np.roll(data[i, :], -ind)
    return dataShifted


def shift_dataSets(shots, data, triggers, signs):
    """ Uses a different trigger for different datasets.
    
    Parameters
    ----------
    shots : array-like
        The number of shots in each data set.
    data : array-like
        Data array where each row is a trace.
    trigger : array-like
        The an array of different trigger values
    signs : array-like
        Trigger on 0 rising, 1 falling, or 2 either.
        
    Returns
    -------
    dataShifted : array-like
        The array of shifted data.
    """
    shape = np.shape(data)
    dataShifted = np.zeros(shape, dtype='double')
    start = 0
    for i in range(len(shots)):
        end = start + shots[i]
        temp = data[start:end]
        dataShifted[start:end] = shift_data(temp, triggers[i], signs[i])
        start += shots[i]
    return dataShifted


def shift_time(dataSets):
    """ Adjust the time to match the shifted data.
    
    Parameters
    ----------
    dataSets : array-like
        The array of data set numbers.
        
    Returns
    -------
    timeShifted : array-like
        The array of shifted time.
    """
    meta = file.load_TRACE('TDS2024C', dataSets[0], 1)['meta']
    timeShifted = np.arange(0.0, meta["Sampling interval"]*meta["Data points"],
                            meta["Sampling interval"])
    return timeShifted

