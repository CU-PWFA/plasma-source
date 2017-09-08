#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, nonecheck=False, overflowcheck=False
#cython: linetrace=True, binding=True
"""
Created on Fri Sep  1 14:07:46 2017

@author: robert
"""

import os
import numpy as np
cimport numpy as np
from numpy.fft import fftfreq
from cython.parallel import prange
from libc.stdio cimport FILE, fwrite, fopen, fclose
import pyfftw


# Load necessary C functions
cdef extern from "complex.h" nogil:
    double complex cexp(double complex)
    double complex csqrt(double complex)


cdef double complex I = 1j


cdef void write_data(double complex *data, char *filename, int N) nogil:
    """ Creates a new binary file and writes the data to it.

    Parameters
    ----------
    data : 
        Pointer to the data to be written.
    filename : string
        The name of the file to be written.
    N : int
        The number of elements in the array to be written.
    """
    cdef FILE *file = fopen(filename, 'a')
    fwrite(data, sizeof(double complex), N, file)
    fclose(file)


def fourier_prop2(double complex[:, :] E, double[:] x, double[:] y, 
                  double[:] z, double lam, path, double n=1):
    """ Propagates an electromagnetic wave from a 2D boundary.

    Uses the Rayleigh-Sommerfeld transfer function to propagate an
    electromagnetic wave from a 2D boundary. The calculation assumes a
    homogeneous index of refraction in the region of propagation. This
    implementation will save a file for each z-slice in the folder params.

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
        The electric field at position (Z, x, y).
    """
    # Need these to calculate frequencies later on
    cdef int i, j, k
    cdef int Nx = np.size(x)
    cdef int Ny = np.size(y)
    cdef int Nz = np.size(z)
    cdef int N = Nx * Ny
    cdef double X = x[Nx-1] - x[0]
    cdef double Y = y[Ny-1] - y[0]
    cdef double dx = X / (Nx-1)
    cdef double dy = Y / (Ny-1)
    # Initialize fftw for the FFT
    e = pyfftw.empty_aligned((Nx, Ny), dtype='complex128')
    cdef double complex[:, :] e_view = e
    pfft = pyfftw.builders.fft2(e, overwrite_input=True, avoid_copy=True, threads=4)
    pifft = pyfftw.builders.ifft2(e, overwrite_input=True, avoid_copy=True, threads=4)
    # Fourier transform the electric field on the boundary
    cdef double complex[:, :] eb = np.zeros((Nx, Ny), dtype='complex128')
    eb = pfft(E)
    # Pre-calculate the spatial frequencies
    cdef double complex[:, :] fz = np.zeros((Nx, Ny), dtype='complex128')
    cdef double[:] fx2 = fftfreq(Nx, dx)**2
    cdef double[:] fy2 = fftfreq(Ny, dy)**2
    cdef double complex pre = I*2*np.pi
    cdef double fmax2 = (n/lam)**2
    with nogil:
        for i in prange(Nx):
            for j in range(Ny):
                fz[i, j] = pre*csqrt(fmax2 - fx2[i] - fy2[j])
    for i in range(Nz):
        with nogil:
            for j in prange(Nx):
                for k in range(Ny):
                    e_view[j, k] = eb[j, k] * cexp(fz[j, k]*z[i])
        e = pifft()
        np.save('electricField'+str(i)+'.npy', e)
    return e
