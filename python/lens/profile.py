#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 08:49:19 2017

@author: robert
"""

import numpy as np
from ionization import ionization


def intensity_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, ion, N, zf=0):
    """ Creates an intensity profile for a uniform plasma with Gaussian ends.

    Returns the on-axis intensity profile necessary to create a fully ionized
    on-axis plasma with Gaussian entrance and exit density ramps. The
    calculation assumes a 

    Parameters
    ----------
    z0 : double
        The distance at which the uniform fully ionized plasma starts.
    dz : double
        The length of the fully ionized plasma.
    sigmaIn : double
        The length of the input ramp, sigma for the Gaussian.
    sigmaOut : double
        The length of the output ramp, sigma for the Gaussian.
    ion : dictionary
        atom : dictionary
            See the description in ionization.ionization for details.
        tau : double
            Temporal length of the pulse in fs.
        type : string
            The type of pulse, i.e. gaussian, flattop, etc.
    N : int
        Number of grid points in z.

    Returns
    -------
    z : array-like
        Array of distances along the optical axis the intensity is returned at.
    I : array-like
        Intensity profile given at each point in z.
    """
    # Create the z grid
    d1 = 0.0
    d2 = z0
    d3 = z0 + dz
    if zf is not 0:
        d4 = zf
    else:
        d4 = d3+10*sigmaOut
    Z = d4 - d1
    z = np.linspace(d1, d1+Z, N)
    # Create the density profile
    frac = np.zeros(N)
    # 0.999 prevents going outside of the interpolating functions range
    peak = 0.999
    sel = z <= d2
    frac[sel] = peak*np.exp(-(z[sel]-d2)**2/(2*sigmaIn**2))
    sel = np.array(z > d2) * np.array(z < d3)
    frac[sel] = peak
    sel = z >= d3
    frac[sel] = peak*np.exp(-(z[sel]-d3)**2/(2*sigmaOut**2))
    # Get the required intensity
    I = ionization.intensity_from_density(ion, frac)
    return z, I
    