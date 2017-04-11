#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 21:28:31 2017

@author: rariniello
"""

import numpy as np
from numpy.lib.scimath import sqrt
from numpy.fft import fft, ifft, fft2, ifft2, fftfreq


def fourier_prop(E, x, z, lam, n=1):
    """ Propogates an electromagnetic wave from a 1D boundary.

    Uses the Rayleigh-Sommerfeld transfer function to propagate an
    electromagnetic wave from a 1D boundary. The calculation assumes a
    homogeneous index of refraction in the region of propagation.

    Parameters
    ----------
    E : array_like
        Array of E field values at points x along the boundary.
    x : array_like
        Array of transverse locations in x on the boundary. Elements must be
        evenly spaced for fft.
    z : array_like
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    n : double, optional
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : array_like
        The electric field at position (z, x).
    """
    Nx = np.size(x)
    Nz = np.size(z)
    X = x[Nx-1]-x[0]
    dx = X / (Nx-1)
    # Fourier transform the electric field on the boundary
    eb = fft(E)
    # Calculate the transfer function at every point in Fourier space
    fx = fftfreq(Nx, dx)
    fz = sqrt((n/lam)**2 - fx**2)
    z = np.reshape(z, (Nz, 1))
    e = np.exp(1j*2*np.pi*z*fz) * eb
    # Inverse Fourier transform to get the real-space field
    e = ifft(e, axis=1)
    return e


def fourier_prop2(E, x, y, z, lam, n=1):
    """ Propagates an electromagnetic wave from a 2D boundary.

    Uses the Rayleigh-Sommerfeld transfer function to propagate an
    electromagnetic wave from a 2D boundary. The calculation assumes a
    homogeneous index of refraction in the region of propagation.

    Parameters
    ----------
    E : array_like
        Array of E field values at points (x, y) along the boundary.
    x : array_like
        Array of transverse locations in x on the boundary. Elements must be
        evenly spaced for fft.
    y : array_like
        Array of transverse locations in y on the boundary. Elements must be
        evenly spaced for fft.
    z : array_like
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    n : double, optional
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : array_like
        The electric field at position (z, x, y).
    """
    # Need these to calculate frequencies later on
    Nx = np.size(x)
    Ny = np.size(y)
    Nz = np.size(z)
    X = x[Nx-1] - x[0]
    Y = y[Ny-1] - y[0]
    dx = X / (Nx-1)
    dy = Y / (Ny-1)
    # Fourier transform the electric field on the boundary
    eb = fft2(E)
    # Calculate the transfer function at every point in Fourier space
    fx2 = np.reshape(fftfreq(Nx, dx)**2, (Nx, 1))
    fy2 = fftfreq(Ny, dy)**2
    fz = sqrt((n/lam)**2 - fx2 - fy2)
    z = np.reshape(z, (Nz, 1, 1))
    e = np.exp(1j*2*np.pi*z*fz) * eb
    # Inverse fourier transform to get the real-space field
    e = ifft2(e)
    return e


def beam_prop(E, nih, x, z, lam, nh):
    """ Propogates an electromagnetic wave from a 1D boundary.

    Propogates a beam from a 1D boundary through a region with an inhomogenous
    index of refraction. This algorithm assumes that the variation in index is
    small enough that their is not significant refraction between grid points
    and that diffraction can be calculated solely with nh.
    The function uses the Rayleigh-Sommerfeld transfer function to propagate an
    electromagnetic wave between grid points and then refracts off of the index
    of refraction change at the grid point.

    The total index of refraction is given by n = nh + nih, where nh is the
    homogenous index of refraction and nih is the inhomogenous variation from
    nh.

    Parameters
    ----------
    E : array_like
        Array of E field values at points x along the boundary.
    nih : array_like
        Variation in index of refraction from nh passed as (x, z).
    x : array_like
        Array of transverse locations in x on the boundary. Elements must be
        evenly spaced for fft.
    z : array_like
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    nh : double
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : array_like
        The electric field at position (z, x).
    """
    Nx = np.size(x)
    Nz = np.size(z)
    e = np.zeros((Nz, Nx), dtype=np.complex)
    e[0, :] = E
    T = 1
    for i in range(1, Nz):
        dz = z[i] - z[i-1]
        # Fourier propogate to the next z point
        e[i] = fourier_prop(e[i-1], x, [dz], lam, nh)
        e[i] = e[i] * T
        # Phase transmission function for refraction
        T = np.exp(1j * np.pi * (nih[:, i-1]+nih[:, i]) * dz)
    return e
