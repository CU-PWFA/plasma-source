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
    
    # Calculate differentials
    dv = np.diff(v) / np.diff(tau)
    dv = np.append(dv, dv[-1])

    # Preallocate for loop
    v_int    = np.zeros(len(tau_int))
    v1n      = np.zeros(len(v))
    v1n[0]   = dv[0]
    v_int[0] = v[0] 
    # Initialize indexing variables
    ind0 = 0
    k    = 0
    # Loop and interpolate 
    for i in range(1, len(tau_int)):
        tau_in = tau_int[i]
        ind    = np.argwhere(tau <= tau_in)[-1]
        if ind != ind0:
            v1n[k]  = dv[ind]
            ind0 = ind
            k    += 1 
        tn = tau[ind]
        v_int[i] = v[ind] + dv[ind] * (tau_in - tn)
    v1n[-1] = v1n[-2]
    return v1n, v_int
def position_interp(x, tau, tau_int, n_pad = 2):
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
    n_pad : int
        The number of padded entries to tau

    Returns
    -------
    x1n : array_like
        Array of 1st interpolation coefficient for every time step in tau
    x2n : array_like
        Array of 2nd interpolation coefficient for every time step in tau
    x_int : array_like
        Array of interpolated positions corresponding to tau_int
    """
    # Calculate differentials
    dx = np.diff(x) / np.diff(tau)
    dx = np.append(dx, dx[-1])
    dx2 = np.diff(dx) / np.diff(tau)
    dx2 = np.append(dx2, dx2[-1])
    
    # Preallocate for loop
    x_int    = np.zeros(len(tau_int))
    x1n      = np.zeros(len(x))
    x2n      = np.zeros(len(x))
    x1n[0]   = dx[0]
    x2n[0]   = dx2[0]
    x_int[0] = x[0]
    # Initialize indexing variables
    ind0 = 0
    k    = 0

    # Loop and interpolate
    for i in range(1, len(tau_int)):
        tau_in = tau_int[i]
        ind    = np.argwhere(tau <= tau_in)[-1]
        if ind != ind0:
            x1n[k] = dx[ind]
            x2n[k] = 0.5 * dx2[ind]
            ind0   = ind
            k     += 1
        tn = tau[ind]
        x_int[i] = x[ind] + dx[ind] * (tau_in - tn) \
                          + 0.5 * dx2[ind] * (tau_in - tn)**2

    x1n[-1] = x1n[-2]
    x2n[-1] = x2n[-2]
    return x1n, x2n, x_int
