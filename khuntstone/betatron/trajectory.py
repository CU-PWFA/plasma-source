# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:14:47 2020
Module to interpolate particle velocity through a PWFA. Used in betatron 
radiation calculations. 
@author: keenan
"""

# Python libraries
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as ius

def get_ders(tau, var):
    s = ius(tau, var)
    xp = s.derivative(n = 1)
    xpp = s.derivative(n = 2)
    return xp(tau), xpp(tau)

def quad_interp(var, var1n, var2n, tau, tau_int):
    var_int = np.zeros(len(tau_int))
    for i in range(len(tau_int)):
        ind = np.argmin(abs(tau_int[i] - tau))
        taun = tau[ind]
        var_int[i] = var[ind]\
                     + var1n[ind] * (tau_int[i] - taun) \
                     + var2n[ind] * (tau_int[i] - taun)**2
    return var_int

def lin_interp(var, var1n, tau, tau_int):
    var_int = np.zeros(len(tau_int))
    for i in range(len(tau_int)):
        ind  = np.argmin(abs(tau_int[i] - tau))
        taun = tau[ind]
        var_int[i] = var[ind] + var1n[ind] * (tau_int[i] - taun)
    return var_int

def traj_interp(x, v, tau, tau_int):
    """
    Function to perform quadratic interpolation of a particle's position
    
    Parameters:
    -----------
    x       : array_like
              Array of particles 4-position (ct, x, y, z)
    v       : array_like
              Array of particles 4 velocity (c*gamma, vx, vy, vz)
    tau     : array_like
              Array of proper time corresponding to x
    tau_int : array_like
              Arrray of proper times to interpolate over
              
    Returns:
    --------
    x41n  : array_like
            Array of linear coefficients (ct1n, x1n, y1n, z1n)
    x42n  : array_like
            Array of quadratic coefficients (ct2n, x2n, y2n, z2n)
    v1n   : Array of linear velocity coefficients (vx1n, vy1n, vz1n)
    x_int : array_like
            Array of the interpolated trajectory (ct_int, x_int, y_int, z_int)
    v_int : array_like
            Array of the interpolated velocity (vx_int, vy_int, vz_int)
    """
    # Position Interpolation
    # ct
    ct1n, ct2n = get_ders(tau, x[0])
    ct2n   = 0.5 * ct2n
    ct_int = quad_interp(x[0], ct1n, ct2n, tau, tau_int)
    # x
    x1n, x2n = get_ders(tau, x[1])
    x2n      = 0.5 * x2n
    x_int    = quad_interp(x[1], x1n, x2n, tau, tau_int)
    # y 
    y1n, y2n = get_ders(tau, x[2])
    y2n      = 0.5 * y2n
    y_int    = quad_interp(x[2], y1n, y2n, tau, tau_int)
    # z
    z1n, z2n = get_ders(tau, x[3])
    z2n      = 0.5 * z2n
    z_int    = quad_interp(x[3], z1n, z2n, tau, tau_int)
    
    x41n   = (ct1n, x1n, y1n, z1n)
    x42n   = (ct2n, x2n, y2n, z2n)
    x_int = (ct_int, x_int, y_int, z_int)
    
    # Velocity interpolation
    # vx
    vx1n = get_ders(tau, v[1])[0]
    vx_int = lin_interp(v[1], vx1n, tau, tau_int)
    # vy
    vy1n = get_ders(tau, v[2])[0]
    vy_int = lin_interp(v[2], vy1n, tau, tau_int)
    # vz
    vz1n = get_ders(tau, v[3])[0]
    vz_int = lin_interp(v[3], vz1n, tau, tau_int)
    
    v_int = (vx_int, vy_int, vz_int)
    v1n   = (vx1n, vy1n, vz1n)
    
    return x41n, x42n, v1n, x_int, v_int
        



