# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 21:14:28 2017

@author: robert
"""

import numpy as np
from scipy import integrate
from ionization import ionization


def phase_function(I0, r, I, z, lam):
    """ Generates a phase mask to create an on axis intensity profile.

    Parameters
    ----------
    I0 : array-like
        Intensity at the lens, each element corresponds to an element in r, in
        10^14 W/cm^2. (Not sure the units matter, as long as I0 and I are the
        same I think it will work)
    r : array-like
        Array of radius values the input beam is specified at. Must be evenly
        spaced for the Fourier transform and start at 0.
    I : array-like
        Array of desired intensity profile along the optical axis specified
        at each point in z, in 10^14 W/cm^2.
    z : array-like
        Array of the distances in z from the phase mask. Must be evenly spaced
        for the Fourier transform.
    lam : double
        Wavelength of the incident light.
    """
    N = np.size(r)
    E0 = ionization.field_from_intensity(I0)
    Ez = ionization.field_from_intensity(I)
    # Find the functions that are a Fourier-real space pair
    chi = r**2
    eta = np.append(np.zeros(N), E0)
    eps = 1j*lam*z*np.exp(-1j*2*np.pi*z/lam)/np.pi * Ez
    # Frequency of the Fourier transform
    f = np.pi**2 / (2*z)
    