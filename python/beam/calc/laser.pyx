#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:15:19 2017

@author: robert
"""

#import numpy as np
cimport numpy as np
from cython.parallel import prange

# Load necessary C functions
cdef extern from "complex.h" nogil:
    double complex cexp(double complex)
    double complex csqrt(double complex)


cdef double complex I = 1j


def fourier_prop(double complex[:, :] E, double[:] x, double[:] y, double[:] z,
                 double lam, double n, fft, ifft, path):
    """ Propagates an electromagnetic wave from a 2D boundary to an array of z.

    Uses the Rayleigh-Sommerfeld transfer function to propagate an
    electromagnetic wave from a 2D boundary. The calculation assumes a
    homogeneous index of refraction in the region of propagation.

    Parameters
    ----------
    E : double complex[:, :]
        Array of E field values at points (x, y) along the boundary.
    x : double[:]
        Array of transverse locations in x on the boundary. Elements must be
        evenly spaced for fft.
    y : double[:]
        Array of transverse locations in y on the boundary. Elements must be
        evenly spaced for fft.
    z : double[:]
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    n : double, optional
        Index of refraction of the medium the wave is propagating through.
    fft : function
        The fft scheme, an object of the pyfftw.FFTW class.
    ifft : function
        The ifft scheme, an object of the pyfftw.FFTW class.

    Returns
    -------
    e : double complex[:, :]
        The electric field at position (Z, x, y).
    """


cdef double complex[:, :] fourier_step(double complex[:, :] E,
                   double complex[:, :] ikz, double dz, fft, ifft):
    """ Propagates a field across a single step of length dz.
    
    A lightweight version of fourier_prop meant to be integrated into other
    algorithms such as the split step method.
    
    Parameters
    ----------
    E : double complex[:, :]
        Array of initial E field values at points (x, y).
    ikz : double complex[:, :]
        Array of spatial frequencies in z, see ikz_RS.
    dz : double
        Size of the step to propagate.
    fft : function
        The fft scheme, an object of the pyfftw.FFTW class.
    ifft : function
        The ifft scheme, an object of the pyfftw.FFTW class.

    Returns
    -------
    e : double complex[:, :]
        The electric field after propagating a distance dz.
    """


cdef double complex[:, :] ikz_RS(double[:] fx, double[:] fy, double lam,
                   double n):
    """ Calculates i*kz for the Rayleigh-Sommerfeld transfer function.
    
    Parameters
    ----------
    fx : double[:, :]
        The spatial frequencies in the x direction.
    fy : double[:, :]
        The spatial frequencies in the y direction.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    n : double
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    fz : double complex[:, :]
        The spatial frequency in the z direction.
    """


cdef double complex[:] fresnel_axis(double complex[:], double[:] r,
                   double[:] z, double lam, double n):
    """ Returns the electric field along the optical axis in the Fresnel limit.

    Uses the Fresnel diffraction integral to calculate the electric field along
    the optical axis resulting from a cylindrically symmetric felectric field.
    Not this is only valid in the paraxial limit. 

    Parameters
    ----------
    E : double complex[:]
        Array of E field values at points r on the boundary, E(r).
    r : double[:]
        Array of radial locations in on the boundary. Doesn't need to be evenly
        spaced.
    z : double[:]
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    n : double
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : double complex[:]
        The electric field at position (z).
    """
