#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, nonecheck=False
#cython: overflowcheck=False, cdivision=True
#cython: linetrace=False, binding=False
"""
Created on Mon Dec 11 10:58:22 2017

@author: robert
"""

import numpy as np
cimport numpy as np
from numpy.fft import fftfreq
from scipy import integrate
from cython.parallel import prange

cdef extern from "math.h" nogil:
    double sin(double)
    double cos(double)
    double sqrt(double)

def electron_propagation_plasma(double[:, :] ptcls, double[:] z, double z0, 
                                double[:] ne, int dumpPeriod, saveP, dgammadz):
    """ Propagate an electron beam through an ion column.
    
    Propagates a collection of macro particles through a full blowout plasma
    wake. The calculation essentially assumes the electrons are propagating
    through a pure ion column and allows for energy increase/decrease.
    """
    cdef int i, j
    cdef int N = np.shape(ptcls)[0]
    cdef int Nz = len(z)
    cdef double dgamma
    cdef double kp, kb, dz
    cdef double coskb, sinkb, angle
    cdef double R11, R12, R21, R22
    # Calculate parameters for each z-slice
    for i in range(Nz-1):
        kp = 5.95074e4 * sqrt(ne[i])
        # Pre-calculate the energy gain per slice
        dz = z[i+1] - z[i]
        dgamma = dgammadz(ne[i]) * dz
        with nogil:
            for j in prange(N):
        #if True:
        #    for j in range(N):
                ptcls[j, 5] += 0.5*dgamma
                kb = kp/sqrt(2*ptcls[j, 5])
                coskb = cos(kb*dz)
                sinkb = sin(kb*dz)
                angle = 1 - dgamma / ptcls[j, 5]
                # Calculate the components of the transfer matrix
                if ne[i] == 0.0:
                    R11 = 1.0
                    R12 = dz
                    R21 = 0.0
                    R22 = 1.0
                else:
                    R11 = coskb
                    R12 = sinkb / kb
                    R21 = -angle * kb * sinkb
                    R22 = angle * coskb
                ptcls[j, 0] = R11 * ptcls[j, 0] + R12 * ptcls[j, 1]
                ptcls[j, 1] = R21 * ptcls[j, 0] + R22 * ptcls[j, 1]
                ptcls[j, 2] = R11 * ptcls[j, 2] + R12 * ptcls[j, 3]
                ptcls[j, 3] = R21 * ptcls[j, 2] + R22 * ptcls[j, 3]
                ptcls[j, 5] += 0.5*dgamma
        if (i % dumpPeriod) == 0:
            saveP(ptcls, z[i]+z0)
    return np.array(ptcls)
