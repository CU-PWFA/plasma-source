#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:43:09 2019
Module for interpolating particle velocity and position in a plasma for betatron
radiation calculations
@author: keenan
"""

# Standard imports
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import sys
# WARGSim imports
sys.path.insert(0, "/home/keenan/WARGSim")
from beams import electronbeam
# Set path for dumping
path = "/media/keenan/Data_Storage/WARGSim/Dump"
# pretty plots
plt.style.use('huntstone')


def velocity_interp(v, tau, tau_int):
    """ 
    Function to perform a linear 1D interpolation of a particle's velocity
        
    Parameters
    ----------
    v : array_like
        Array of particle's velocity
    tau : array_like
          Proper time array corresponding to v
    tau_int : array_like
        Array of query proper times for interpolation
        
    Returns
    -------
    v1n : array_like
        Array of interpolation coefficients for every time step in tau
    v_int : array_like
        Array of interpolated velocities for each time step in tau
    """
    
    # Get indices for points in tau surrounding points in tau_int
    
    ind1 = np.squeeze(np.array([np.argwhere(tau >=i)[0] - 1 for i in tau_int]))
    ind1[np.argwhere(ind1 < 0)] = 0
    ind2 = ind1 + 1
    v0   = v[ind1]
    v1   = (v[ind2] - v[ind1]) / (tau[ind2] - tau[ind1])
    v1n  = np.append(np.unique(v1), v1[-1])
    
    v_int = v0 + v1 * (tau_int - tau[ind1]);
    
    return v1n, v_int

def position_interp(x, tau, tau_int):
    """
    Function to perform a quadratic 1D interpolation of a particle's position
    
    Parameters
    ----------
    x : array_like
        Array of particle's position
    tau : array_like
        Array of proper times corresponding to x
    tau_int : array_like
        Array of query proper times for the interpolation
        
    Returns
    -------
    x1n : array_like
        Array of 1st interpolation coefficient for every time step in tau
    x2n : array_like
        Array of 2nd interpolation coefficient for every time step in tau
    x_int : array_like
        Array of interpolated positions corresponding to tau_int
    """

    x1n   = np.zeros(len(tau_int))
    x2n   = np.zeros(len(tau_int))
    x_int = np.zeros(len(tau_int))
    for i in range(len(tau_int)):
        tau_in = tau_int[i]
        inds   = np.argwhere(tau>= tau_in)[0]
        if inds > 0:
            ind1 = inds[0] - 1
            ind2 = inds[0]
            ind3 = inds[0] + 1
        else:
            ind1 = inds[0]
            ind2 = inds[0] + 1
            ind3 = inds[0] + 2
            
        x1n[i] = (x[ind2] - x[ind1]) / (tau[ind2] - tau[ind1])
        if ind3 >= len(x):
            x_additional = x[-1] + x1n[i]
            tau_additional = tau[-1] + (abs(tau[1] - tau[2]))
            x2n[i] = ((x_additional - x[ind2]) / (tau_additional - tau[ind2])) - x1n[i]
        else:
            x1n[i] = (x[ind2] - x[ind1]) / (tau[ind2] - tau[ind1])
            x2n[i] = ((x[ind3] - x[ind2]) / (tau[ind3] - tau[ind2])) - x1n[i]
            
        x_int[i] = x[ind1] + x1n[i] * (tau_in - tau[ind1]) + x2n[i] * \
                   (tau_in - tau[ind1])**2
                   
    x1n = np.append(np.unique(x1n), x1n[-1])
    x2n = np.append(np.unique(x2n), x2n[-1])
    return x1n, x2n, x_int
    
    
    
    

