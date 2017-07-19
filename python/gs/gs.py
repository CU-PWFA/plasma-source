# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 15:49:40 2017

@author: robert
"""

import numpy as np
from numpy.fft import fft, ifft, fft2, ifft2
from numpy.random import rand


def gs(source, target, n=100, err=None, zeroPhase=None):
    """ The 1D Gerchberg-Saxton algorithm between a source and a target.

    Given two amplitude distributions related by a Fourier transform, this will
    find the phase required at the source such that the Fourier transform of
    the source is the approximatly the target amplitude with arbitrary phase.
    Note that the algorithm is guaranteed to converge but not necessarily to
    the global minimum.

    Parameters
    ----------
    source : array-like
        The source amplitude distribution. If complex, the algorithm will use
        the starting phase.
    target : array-like
        The target amplitude distribution. Only the amplitude will be used.
    n : double, optional
        Number of iterations.
    err : double, optional
        Specify a level of convergence, the algorithm will converge once it has
        reached this level or if n is reached.
    zeroPhase : bool, optional
        Pass as True to start the phase on the source at 0.

    Returns
    -------
    phi : array-like
        The phase to be applied to the source.
    """
    size = np.size(source)
    if zeroPhase is not None:
        theta = np.zeros(size)
    elif np.iscomplexobj(source):
        theta = np.angle(source)  # Initial angle
    else:
        theta = 2*np.pi*rand(size)
    source = abs(source) * np.exp(1j*theta)
    target = abs(target)
    # Run for the specified number of iterations
    if err is None:
        for i in range(0, n):
            targetPhi = np.angle(fft(source))
            target = abs(target) * np.exp(1j*targetPhi)
            phi = np.angle(ifft(target))
            source = abs(source) * np.exp(1j*phi)
    # Run until the phase converges to a specified level
    else:
        i = 0
        delta = 2*np.pi
        phiOld = np.zeros(size)
        while delta > err and i < n:
            targetPhi = np.angle(fft(source))
            target = abs(target) * np.exp(1j*targetPhi)
            phi = np.angle(ifft(target))
            source = abs(source) * np.exp(1j*phi)
            delta = np.amax(abs(phiOld-phi))
            phiOld = phi
            i += 1
            if i == n:
                print("Convergence not reached in n=", n, "iterations.")
    return np.unwrap(phi)


def gs2(source, target, n=100, err=None, zeroPhase=None):
    """ The 2D Gerchberg-Saxton algorithm between a source and a target.

    Given two amplitude distributions related by a Fourier transform, this will
    find the phase required at the source such that the Fourier transform of
    the source is the approximatly the target amplitude with arbitrary phase.
    Note that the algorithm is guaranteed to converge but not necessarily to
    the global minimum.

    Parameters
    ----------
    source : array-like
        The source amplitude distribution. If complex, the algorithm will use
        the starting phase.
    target : array-like
        The target amplitude distribution. Only the amplitude will be used.
    n : double, optional
        Number of iterations.
    err : double, optional
        Specify a level of convergence, the algorithm will converge once it has
        reached this level or if n is reached.
    zeroPhase : bool, optional
        Pass as True to start the phase on the source at 0.

    Returns
    -------
    phi : array-like
        The phase to be applied to the source.
    """
    shape = np.shape(source)
    if zeroPhase is not None:
        theta = np.zeros(shape)
    elif np.iscomplexobj(source):
        theta = np.angle(source)  # Initial angle
    else:
        theta = 2*np.pi*rand(shape[0], shape[1])
    source = abs(source) * np.exp(1j*theta)
    target = abs(target)
    # Run for the specified number of iterations
    if err is None:
        for i in range(0, n):
            targetPhi = np.angle(fft2(source))
            target = abs(target) * np.exp(1j*targetPhi)
            phi = np.angle(ifft2(target))
            source = abs(source) * np.exp(1j*phi)
    # Run until the phase converges to a specified level
    else:
        i = 0
        delta = 2*np.pi
        phiOld = np.zeros(shape)
        while delta > err and i < n:
            targetPhi = np.angle(fft2(source))
            target = abs(target) * np.exp(1j*targetPhi)
            phi = np.angle(ifft2(target))
            source = abs(source) * np.exp(1j*phi)
            delta = np.amax(abs(phiOld-phi))
            phiOld = phi
            i += 1
            if i == n:
                print("Convergence not reached in n=", n, "iterations.")
    return np.unwrap(phi)
