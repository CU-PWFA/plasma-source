#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:33:36 2019

Module to compute the spectral intensity emitted by charged particles following
an arbitray trajectory.

@author: keenan
"""
import numpy as np
import v_x_interp as vx

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
        A 3 x n array of the particle's velocity at each time step
    tau_arr : array_like
        A 1 x n array of the proper time corresponding to x_arr and v_arr
    tau_int : array_like
        An array of query proper times to interpolate for
        
    Returns
    -------
    x0n : array_like
        A 4 x n array of the zeroth order interpolation coefficient of the 
        4-position
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
    x00n, x01n, x02n, x0_int = vx.position_interp(x_arr[0], tau_arr, tau_int)
    x10n, x11n, x12n, x1_int = vx.position_interp(x_arr[1], tau_arr, tau_int)
    x20n, x21n, x22n, x2_int = vx.position_interp(x_arr[2], tau_arr, tau_int)
    x30n, x31n, x32n, x3_int = vx.position_interp(x_arr[3], tau_arr, tau_int)
    
    x0n = (x00n, x10n, x20n, x30n)
    x1n = (x01n, x11n, x21n, x31n)
    x2n = (x02n, x12n, x22n, x32n)
    x_int = (x0_int, x1_int, x2_int, x3_int)
    
    # Interpolate velocity
    v11n, v1_int = vx.velocity_interp(v_arr[0], tau_arr, tau_int)
    v12n, v2_int = vx.velocity_interp(v_arr[1], tau_arr, tau_int)
    v13n, v3_int = vx.velocity_interp(v_arr[2], tau_arr, tau_int)
    
    v1n   = (v11n, v12n, v13n)
    v_int = (v1_int, v2_int, v3_int)
    
    return x0n, x1n, x2n, v1n, x_int, v_int
    
def get_vars(kappa, x1n, x2n, dtau, v_arr):
    """
    Computes reduced variables in Thomas equations 14-19 from 4-wave vector, 
    interpolation coefficients, and initial time step
    
    Parameters
    ----------
    kappa : array_like
        4 entry array of the wave vector omega * (1/c, s_hat/c) where s_hat
        is the observation unit vector
    x1n : array_like
        Array of 1st order interpolation coeffs for the jth particle at 
        all n time steps, size is 4 x n (n steps in ct, x, y, and z)
    x2n : array_like 
        Array of 2nd order interpolation coeffs for the jth particle at all n
        time steps, size is the same as x1jn
    dtau : float
        The time step of the original particle trajectory [s]
    v_arr : array_like
        A 3 x n array of the particle's velocity at each time step
    
    Returns
    -------
    chi_1jn : array_like
        Reduced variable chi_1jn (the dot product of the 4-wave vector and the
        4-position vector 1st interpolation coefficient) (s^-1)
    chi_2jn : array_like
        Reduced variable chi_2jn (the dot product of the 4-wave vector and the 
        4-position vecotr 2nd interpolation coefficient) (s^-2)
    phi_p : array_like
        3 x n vector reduced variable (m/s^2)
    phi_m  : array_like
        3 x n vector reduced variable (m/s^2)
    """
    
    chi_1n = kappa[0] * x1n[0] - kappa[1] * x1n[1] - kappa[2] * x1n[2] \
             - kappa[3] * x1n[3]
    chi_2n = kappa[0] * x2n[0] - kappa[1] * x2n[1] - kappa[2] * x2n[2] \
             - kappa[3] * x2n[3]
    
    phi_p = (dtau**2 / 4) * chi_2n + (dtau / 2) * chi_1n
    phi_m = (dtau**2 / 4) * chi_2n + (dtau / 2) * chi_1n
    
    
    
    
    return chi_1n, chi_2n

