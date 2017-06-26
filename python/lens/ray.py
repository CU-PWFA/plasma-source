#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:44:21 2017

@author: rariniello
"""

import numpy as np
from scipy import integrate


def uniform_phase(I, z, R, r0=0):
    """ Generates a phase mask to create a on axis intensity profile.

    Generates a phase mask to produce a passed longitudinal on axis intensity
    profile assuming a uniform (flat top) mask illumination. The calculation
    uses pure ray tracing to create the phase mask and does not account for
    diffraction effects.

    Parameters
    ----------
    I : array-like
        Array of desired intensity profile along the optical axis.
    z : array-like
        Array of the distances in z from the phase mask.
    R : double
        Radius of the phase mask / input beam.
    r0 : double, optional
        Starting radius of the phase mask, places a hole of radius r0 in the
        middle of the phase mask.

    Returns
    -------
    I0 : double
        The required input intensity to achieve the desired output.
    r : array-like
        Radius vector where the phase is specified at. Each radius refracts
        rays to the corresponding element in z.
    phi : array-like
        Phase of the mask at each r. Note that this phase must be multiplied by
        k to get the actual phase delay that goes in the argument of the
        complex exponential.
    """
    # Set up the arrays to store everything
    Nz = np.size(z)
    r = np.zeros(Nz)
    phi = np.zeros(Nz)
    r[0] = r0
    # Calculate the required input intensity
    I0 = integrate.trapz(I, z) / (np.pi * R**2)
    # Calculate r and phi incrementally
    sinnew = r0/np.sqrt(r0**2 + z[0]**2)
    for i in range(1, Nz):
        Iavg = (I[i] + I[i-1]) / 2
        dz = z[i] - z[i-1]
        r[i] = np.sqrt(Iavg*dz/(np.pi*I0) + r[i-1]**2)
        dr = r[i] - r[i-1]
        sinold = sinnew
        sinnew = r[i] / np.sqrt(r[i]**2 + z[i]**2)
        phi[i] = phi[i-1] - (sinnew + sinold)*dr/2
    # Return everything
    return I0, r, phi
