#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:22:51 2017

@author: robert
"""

import numpy as np
from scipy import integrate
from scipy.special import jn


def ht(f, x, k):
    """ Forward Hankel transform through direct integration.

    Calculates the forward Hankel transform by directly evaluating the
    integrals over the Bessel function kernels.

    Parameters
    ----------
    f : array-like
        Array of function values in normal space.
    x : array-like
        Array of points the function is defined at.
    k : array-like
        Array of frequency space points to find the transform at.

    Returns
    -------
    h : array-like
        The Hankel transform at each point in k.
    """
    N = np.size(f)
    M = np.size(k)
    x = np.reshape(x, (1, N))
    f = np.reshape(f, (1, N))
    k = np.reshape(k, (M, 1))
    h = integrate.simps(f*jn(0, x*k)*x, x)
    h = np.reshape(h, (M))
    return h


def iht(h, k, x):
    """ Inverse Hankel transform through direct integration.

    Calculates the inverse Hankel transform by directly evaluating the
    integrals over the Bessel function kernels.

    Parameters
    ----------
    h : array-like
        Array of function values in frequency space.
    k : array-like
        Array of frequency values the function is defined at.
    x : array-like
        Array of normal space points to find the transform at.

    Returns
    -------
    f : array-like
        The function at each point in x.
    """
    N = np.size(h)
    M = np.size(x)
    k = np.reshape(k, (1, N))
    h = np.reshape(h, (1, N))
    x = np.reshape(x, (M, 1))
    f = integrate.simps(h*jn(0, k*x)*k, k)
    f = np.reshape(f, (M))
    return f
