#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:33:36 2019

Module to compute the spectral intensity emitted by charged particles following
an arbitray trajectory.

@author: keenan
"""
import numpy as np
import v_x_interp as interp

def get_coeffs(x_arr, v_arr, tau_arr, tau_int):
    """
    Function to get the interpolation coefficients of a particle's four-position
    and four-velocities
    
    Parameters
    ----------
    x_arr : array_like
        A 4 x n array of the particle's 4-position where n is the number of time
        steps
    v_arr : array_like
        A 4 x n array of the particle's 4-velocity
    tau_arr : array_like
        A 1 x n array of the proper time corresponding to x_arr and v_arr
    tau_int : array_like
        An array of query proper times to interpolate for
        
    Returns
    -------
    x1n : array_like
        A 4 x n array of the 1st order interpolation coefficient of the 
        4-position
    x2n : array_like
        A 4 x n array of the 2nd order interpolation coefficient of the 
        4-position
    v1n : array_like
        A 4 x n array of the 1st order interpolation coefficient of the 
        4-velocity
    
    x_int : array_like
        4 x len(tau_int) array of the interpolated 4 position vector
    
    v_int : array_line
        4 x len(tau_int) array of the interpolated 4 velocity vector
    """
    
    # Interpolate the position
    
    x10n, x20n, x0_int = interp.position_interp(x_arr[0], tau_arr, tau_int)
    x11n, x21n, x1_int = interp.position_interp(x_arr[1], tau_arr, tau_int)
    x12n, x22n, x2_int = interp.position_interp(x_arr[2], tau_arr, tau_int)
    x13n, x23n, x3_int = interp.position_interp(x_arr[3], tau_arr, tau_int)
    
    x1n   = (x10n, x11n, x12n, x13n)
    x2n   = (x20n, x21n, x22n, x23n)
    x_int = (x0_int, x1_int, x2_int, x3_int)
    
    # Interpolate velocity
    v10n, v0_int = interp.velocity_interp(v_arr[0], tau_arr, tau_int)
    v11n, v1_int = interp.velocity_interp(v_arr[1], tau_arr, tau_int)
    v12n, v2_int = interp.velocity_interp(v_arr[2], tau_arr, tau_int)
    v13n, v3_int = interp.velocity_interp(v_arr[3], tau_arr, tau_int)
    
    v1n   = (v10n, v11n, v12n, v13n)
    v_int = (v0_int, v1_int, v2_int, v3_int)
    
    return x1n, x2n, v1n, x_int, v_int
    
def get_chis(s_hat, x1jn, z2jn):
    """
    Computes reduced variables chi_1jn and chi_2jn from observation vector
    and interpolation coeffs
    
    Parameters
    ----------
    s_hat : array_like
        3 entry array of the normalized radiation observation vector (i.e. 
        [0, 0, 1] for on-axis radiation)
    x1jn : array_like
        Array of 1st order interpolation coeffs for the jth particle at 
        all n time steps, size is 3 x n (n steps in x, y, and z)
    x2jn : array_like 
        Array of 2nd order interpolation coeffs for the jth particle at all n
        time steps, size is the same as x1jn
    
    Returns
    -------
    chi_1jn : array_like
        Reduced variable chi_1jn (the dot product of the 4-wave vector and the
        4-position vector 1st interpolation coefficient)
    chi_2jn : array_like
        Reduced variable chi_2jn (the dot product of the 4-wave vector and the 
        4-position vecotr 2nd interpolation coefficient)
    """
    
    

