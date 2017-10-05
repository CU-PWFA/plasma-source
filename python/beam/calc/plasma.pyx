#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, nonecheck=False
#cython: overflowcheck=False, cdivision=True
#cython: linetrace=True, binding=True
"""
Created on Mon Sep 18 16:59:45 2017

@author: robert
"""

import numpy as np
cimport numpy as np
from numpy.fft import fftfreq
from scipy import integrate
from cython.parallel import prange
from beam.calc import laser
from beam.calc.ionization cimport adk_rate_linear

# Load necessary C functions
cdef extern from "complex.h" nogil:
    double complex cexp(double complex)
    double complex csqrt(double complex)
    double cabs(double complex)

cdef extern from "math.h" nogil:
    double exp(double)
    double sqrt(double)


def plasma_refraction(double complex[:, :, :] E, double[:] x, double[:] y,
                      double[:] z, double[:] t, double lam, double n0, 
                      double z0, fft, ifft, saveE, saven, atom):
    """ Propagate a laser pulse through a plasma accounting for refraction.

    Propogates a laser pulse through a region of partially ionized gas. This
    function accounts for refraction from the plasma. It determines the plasma
    density by calculating the ionization that has resulted from each temporal
    piece of the pulse. The results are stored in a file, only the central x-z
    plane is recorded.
    """
    cdef int i, j, k, l
    # TODO abstract this into its own function
    cdef int Nx = len(x)
    cdef int Ny = len(y)
    cdef int Nz = len(z)
    cdef int Nt = len(t)
    cdef double dx = x[1] - x[0]
    cdef double dy = y[1] - y[0]
    cdef double dt = t[1] - t[0]
    # Get the ionization information from atom
    cdef double EI = atom['EI']
    cdef int Z = atom['Z']
    cdef int ll = atom['l']
    cdef int m = atom['m']
    # Plasma density and index of refraction arrays
    cdef double[:, :] n = np.zeros((Nx, Ny), dtype='double')
    cdef double[:, :] nih = np.zeros((Nx, Ny), dtype='double')
    cdef double ngas = atom['alpha'] * 5.0e-8
    cdef double nplasma = plasma_index(1.0, lam) - 1.0
    cdef double nh = 1.0 + ngas*n0
    # Pre-calculate the spatial frequencies
    cdef double[:] fx = fftfreq(Nx, dx)
    cdef double[:] fy = fftfreq(Ny, dy)
    cdef double complex[:, :] ikz = laser.ikz_RS(fx, fy, lam, nh)
    cdef double complex arg
    cdef double rate
    for i in range(1, Nz):
        dz = z[i] - z[i-1]
        arg = 1j*2*np.pi*dz / lam
        for j in range(Nt):
            # Propagate the beam through
            Etemp = laser.fourier_step(E[j, :, :], ikz, dz, fft, ifft)
            for k in range(Nx):
                for l in range(Ny):
                    E[j, k, l] = Etemp[k, l]
            with nogil:
                for k in prange(Nx):
                    for l in range(Ny):
                        E[j, k, l] *= cexp(arg*nih[k, l])
                        # Ionize the gas
                        rate = adk_rate_linear(EI, cabs(E[j, k, l]), Z, ll, m)
                        n[k, l] = n0 - (n0 - n[k, l])*exp(-rate * dt)
                        nih[k, l] = n[k, l]*(nplasma - ngas) + n0*ngas
        saveE(E, z[i]+z0)
        saven(n, i)
        # Reset the plasma density for the next slice
        with nogil:
            for k in prange(Nx):
                for l in range(Ny):
                    n[k, l] = 0.0
                    nih[k, l] = 0.0                           
    return E

# TODO, add a plasma refraction function for plasmas with variable density

# TODO, move this to a helper function file
cpdef double plasma_index(double n, double lam):
    """ Calculates the index of refraction of a plasma.

    Parameters
    ----------
    n : double
        Density of the plasma in 10^17 cm^-3.
    lam : double
        Wavelength of the incident light in um.
    """
    return 1.0 - n * lam*lam * 4.47869e-5
