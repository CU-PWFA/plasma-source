#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:31:51 2017

@author: robert
"""

import numpy as np
from scipy.interpolate import interp1d
from numpy.fft import fft, ifft
from scipy.special import jn


# !!! Work in progress! Don't use for anything important.
def fht(r, f):
    """ Fast Hankel transform, zero order, of x.

    Parameters
    ----------
    r : array-like
        Values of the radius f is specified at. Must start at 0.
    f : array-like
        Function of r to be transformed.

    Returns
    -------
    k : array-like
        Spatial frequencies of the transform.
    h : array-like
        Zero order Hankel transform.
    """
    N = np.size(f)
    # Calculate the sampled regions in both spaces
    R = r[N-1]
    r0 = r[1]-r[0]
    k0 = 2*np.pi/R
    alpha  = np.log(R/r0)/N
    # Create the nonlinear sampling grids
    f = interp1d(r, f)
    e = np.exp(alpha*np.arange(0, N))
    r = r0*e
    k = k0*e
    h = fft(fft(r*f(r)) * ifft(alpha*r0*k0*e*jn(0, r0*k0*e))) / k  
    return k, h


def ifht(k, h):
    """ Inverse fast Hankel transform, zero order, of k.
    
    Parameters
    ----------
    r : array-like
        Values of the radius f is specified at. Must start at 0.
    f : array-like
        Function of r to be transformed.

    Returns
    -------
    k : array-like
        Spatial frequencies of the transform.
    H : array-like
        Zero order Hankel transform.
    """
    N = np.size(h)
    # Calculate the sampled regions in both spaces
    K = k[N-1]
    k0 = k[1]
    r0 = 2*np.pi/K
    alpha = np.log(K/k0)/N
    # Create the nonlinear sampling grids
    h = interp1d(k, h)
    e = np.exp(alpha*np.arange(0, N))
    r = r0*e
    k = k0*e
    f = fft(fft(k*h(k)) * ifft(alpha*r0*k0*e*jn(0, r0*k0*e))) / r  
    return r, f
