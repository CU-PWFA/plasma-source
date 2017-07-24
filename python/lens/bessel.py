#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:14:24 2017

@author: robert
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq, fftshift
from ionization import ionization
from scipy.interpolate import interp1d
from ht import intht


def spectrum_from_axis(E, z):
    """ Returns the spatial spectrum in kz from the electric field on-axis.

    Returns the spatial frequencies, kz, and the spatial spectrum in terms of
    kz, of the electric field along the optical axis.

    Parameters
    ----------
    E : array_like
        Array of complex electric field values along the optical axis.
    z : array-like
        Array of z coordinates along the optical axis. Must be evenly spcaed
        for the FFT.

    Returns
    -------
    kz : array-like
        Array of spatial frequencies in z.
    S : array-like
        Spatial spectrum in terms of kz of the electromagnetic field.
    """
    N = np.size(E)
    dz = z[1] - z[0]
    S = fftshift(fft(E)) / N
    kz = 2*np.pi * fftshift(fftfreq(np.size(z), dz))
    return kz, S


def kr_from_kz(kz, lam, S=None):
    """ Returns the spatial frequencies in r from the frequencies in z.
    
    Parameters
    ----------
    kz : array-like
        Spatial frequencies in z.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    S : array-like, optional
        The spatial spectrum expressed in terms of kz.

    Returns
    -------
    kr : array-like
        The spatial spectrum in terms of kr. Note may be a different length
        than kz, decaying (imaginary) frequencies are removed.
    S : array-like
        Spatial spectrum in terms of kr with imaginary frequency terms removed.
    """
    k = 2*np.pi / lam
    sel = kz <= k
    kr = np.sqrt(k**2 - kz[sel]**2)
    if S is not None:
        S = S[sel]
        return kr, S
    else:
        return kr


def uniform_bessel(params, E, z, n=0):
    """ Calculate the required electric field to create the passed intensity.

    Calculates the electric field necessary to create the passed intensity
    distribution along the optical axis. Unifrom Bessel means that it will
    assign a linear phase to the on axis electric field so that it exists in
    approximatly a single Bessel mode. This means the intensity distribution
    will have a constant width. 

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            N : int
                Number of integration steps in Bessel integrals. Try 1000.
            M : int
                Number of r values to calculate, each requires an integral.
            R : double
                Radius to the first zero of the Bessel function.
            lam : double
                Wavelength of the electromagnetic wave in vacuum.
            rmax : double
                Maximum radius to return the electric field at.
    E : array-like
        Desired on axis field in GV/m.
    z : array-like
        Array of z coordinates along the optical axis. Must be evenly spcaed
        for the FFT.
    n : int, optional
        Order of the Bessel function, defaults to order zero.

    Returns
    -------
    r : array-like
        Array of radius coordinates the electric field is given at.
    E : array-like
        Electric field as a function of r on the boundary.
    """
    #TODO figure out if the normalization is just not enough steps or an error
    lam = params['lam']
    k = 2*np.pi/lam
    kz, S = spectrum_from_axis(E, z)
    # Add the linear phase based of the width of the Bessel beam
    kr0 = 2.4048 / params['R']
    kz0 = np.sqrt(k**2 - kr0**2)
    kz = kz + kz0
    # Calculate the true spatial spectrum
    S = S / kz
    # Calculate the radial spatial frequencies
    kr, S = kr_from_kz(kz, lam, S)
    # Inverse Hankel transform
    krn = np.linspace(0, np.amax(kr), params['N'])
    Sn = interp1d(kr, S, fill_value=(S[-1], 0.0), bounds_error=False)
    Sn = Sn(krn)
    r = np.linspace(0, params['rmax'], params['M'])
    E = intht.ihtn(Sn, krn, r, n)
    return r, E
