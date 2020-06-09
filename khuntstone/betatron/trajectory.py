# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:14:47 2020
Module to interpolate particle velocity through a PWFA. Used in betatron 
radiation calculations. 
@author: keenan
"""

# Python libraries
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import sys

def x_interp(x, tau, tau_int):
    """
    Function to perform quadratic interpolation of a particle's position
    
    Parameters:
    -----------
    x       : array_like
              Array of particles 4-position (ct, x, y, z)
    tau     : array_like
              Array of proper time corresponding to x
    tau_int : array_like
              Arrray of proper times to interpolate over
              
    Returns:
    --------
    x1n : array_like
          Array of linear coefficients (ct1n, x1n, y1n, z1n)
    x2n : array_like
          Array of quadratic coefficients (ct2n, x2n, y2n, z2n)
    x_int : array_like
            Array of the interpolated trajectory (ct_int, x_int, y_int, z_int)
    """
    
    def single_interp(x, tau, tau_int):
        # Interpolates in one-dimension
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
        return x1n, x2n, x_int
        
    # Interpolate each component of the 4-position
    ct1n, ct2n, ct_int = single_interp(x[0], tau, tau_int)
    xp1n, xp2n, xp_int = single_interp(x[1], tau, tau_int)
    y1n, y2n, y_int    = single_interp(x[2], tau, tau_int)
    z1n, z2n, z_int    = single_interp(x[3], tau, tau_int)
    
    
    x1n   = (ct1n, xp1n, y1n, z1n)
    x2n   = (ct2n, xp2n, y2n, z2n)
    x_int = (ct_int, xp_int, y_int, z_int)
        
    return x1n, x2n, x_int

def v_interp(v, tau, tau_int):
    """
    Function to perform a linear interpolation of particle velocity
    
    Parameters:
    -----------
    v       : array_like
              Array of particle 4-velocity gamma * (vx, vy, vz)
    tau     : array_like
              Proper time array corresponding to v
    tau_int : array_like
              Array of proper times for interpolation
    
    Returns:
    --------
    v1n : array_like
          Array of interpolation coefficients
    v_int : array of interpolated velocities. 
    """
    
    def single_interp(v, tau, tau_int):
        # Interps a single velocity
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
        
    # Interpolate for each velocity
    vx1n, vx_int = single_interp(v[0], tau, tau_int)
    vy1n, vy_int = single_interp(v[1], tau, tau_int)
    vz1n, vz_int = single_interp(v[2], tau, tau_int)
    
    v1n = (vx1n, vy1n, vz1n)
    v_int = (vx_int, vy_int, vz_int)

    
    return v1n, v_int



