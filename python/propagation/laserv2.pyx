#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, nonecheck=False, overflowcheck=False
#cython: linetrace=False
"""
Created on Fri Sep  1 14:07:46 2017

@author: robert
"""

import os
import numpy as np
from numpy.fft import fftfreq
from cython.parallel import prange
from libc.stdio cimport FILE, fwrite, fopen, fclose


# Load necessary C functions
cdef extern from "complex.h" nogil:
    double complex cexp(double complex)
    double complex csqrt(double complex)


cdef extern from "fftw3.h":
    # Load constants from fftw
    cdef int FFTW_FORWARD
    cdef int FFTW_BACKWARD
    cdef unsigned FFTW_MEASURE
    cdef unsigned FFTW_ESTIMATE
    # Memory allocation
    void *fftw_malloc(size_t)
    void fftw_free(void *)
    # Define cython structure that point to the fftw_plan structure
    ctypedef struct _fftw_plan:
        pass
    ctypedef _fftw_plan *fftw_plan
    # Load fftw functions: plan, execute, destroy
    fftw_plan fftw_plan_dft_2d(int, int, double complex*, double complex*, int,
                               unsigned)
    void fftw_execute(fftw_plan) nogil
    void fftw_destroy_plan(fftw_plan)


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
        The electric field at position (z, x, y).
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
    filename_py = bytes(path + 'electricField.bin', 'ascii')
    if os.path.isfile(filename_py):
        os.remove(filename_py) # Otherwise we'll append to it
    cdef char *filename = filename_py
    # Initialize fftw for the FFT
    cdef double complex[:, :] e = <double complex[:Nx, :Ny]> fftw_malloc(sizeof(double complex) * N)
    cdef fftw_plan pfft = fftw_plan_dft_2d(Nx, Ny, &(e[0, 0]), &(e[0, 0]), FFTW_FORWARD,
                                           FFTW_ESTIMATE)
    cdef fftw_plan pifft = fftw_plan_dft_2d(Nx, Ny, &(e[0, 0]), &(e[0, 0]), FFTW_BACKWARD,
                                           FFTW_ESTIMATE)
    # Fourier transform the electric field on the boundary
    with nogil:
        for i in prange(Nx):
            for j in range(Ny):
                e[i, j] = E[i, j]
    fftw_execute(pfft)
    cdef double complex[:, :] eb = np.zeros((Nx, Ny), dtype='complex128')
    with nogil:
        for i in prange(Nx):
            for j in range(Ny):
                eb[i, j] = e[i, j]
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
            for j in prange(Nx):
                for k in range(Ny):
                    e[j, k] = eb[j, k] * cexp(fz[j, k]*z[i])
            fftw_execute(pifft)
            write_data(&(e[0, 0]), filename, N)
    fftw_destroy_plan(pfft)
    fftw_destroy_plan(pifft)
    Eret = np.zeros((Nx, Ny), dtype='complex128')
    ebret = np.zeros((Nx, Ny), dtype='complex128')
    eRet = np.zeros((Nx, Ny), dtype='complex128')
    for i in range(Nx):
        for j in range(Ny):
            Eret[i, j] = E[i, j]
            ebret[i, j] = eb[i, j]
            eRet[i, j] = e[i, j]
    fftw_free(&(e[0, 0]))
    return Eret, ebret, eRet
