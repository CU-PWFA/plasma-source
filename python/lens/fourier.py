# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 21:14:28 2017

@author: robert
"""

import numpy as np
from scipy import integrate
from ionization import ionization
from scipy.interpolate import interp1d
from numpy.fft import fft, ifft, fftfreq, ifftshift


#XXX Don't use this, it doesn't converge
def phase_function(I0, r, I, z, lam, M):
    """ Generates a phase mask to create an on axis intensity profile.

    Generates a phase mask to produce a passed longitudinal on axis intensity
    profile from an input intensity profile. The calculation use the Fresnel
    diffraction integral to solve the problem. Only the functional form of the
    the input intensity I0 matters, it will be rescaled to satisfy Parseval's
    theorem.

    Parameters
    ----------
    I0 : array-like
        Intensity at the lens, each element corresponds to an element in r. The
        magnitude is irrelevant, it will be rescaled to satisfy Parseval's
        theorem.
    r : array-like
        Array of radius values the input beam is specified at. Must start at 0.
        This is used to create an interpolating function with I0. The length of
        r is used as the size of the fft.
    I : array-like
        Array of desired intensity profile along the optical axis specified
        at each point in z, in 10^14 W/cm^2.
    z : array-like
        Array of the distances in z from the phase mask. Doesn't have to be
        evenly spaced, but must start at 0. This is used to create an
        interpolating function with I.
    lam : double
        Wavelength of the incident light.
    M : integer
        Number of iterations of Gerchberg-Saxton.

    Returns
    -------
    r : array-like
        Radius vector where the phase is specified at. 
    phi : array-like
        Phase of the mask at each r. 
    """
    N = np.size(r)
    Nz = np.size(z)
    Z0 = z[1]
    Z = z[Nz-1]
    R = r[N-1]
    E0 = ionization.field_from_intensity(I0)
    Ez = ionization.field_from_intensity(I)
    # Create interpolating functions
    E0 = interp1d(r, E0)
    Ez = interp1d(z, Ez)
    # Effective coordinate and frequency of the Fourier transform
    F = z_to_f(Z0, lam)
    f0 = z_to_f(Z, lam)
    f = np.linspace(-F, F, N, False)
    f = ifftshift(f)
    chi = np.linspace(0, R**2, N)
    # Find the functions that are a Fourier-real space pair
    eta = E0(np.sqrt(chi))
    eps = np.zeros(N, dtype=np.complex)
    sel = f > f0
    eps[sel] = 1j*2*np.pi * Ez(f_to_z(abs(f[sel]), lam)) / f[sel]
    # Scale the frequency space function to satisfy Parseval's Theorem
    scale = integrate.trapz(abs(eps)**2, f) / integrate.trapz(abs(eta)**2, chi)
    eta = eta * scale
    # Gerchberg-Saxton algorithm
    # TODO add in error criterion instead of iteration number
    A = fft(eps);
    for i in range(0, M):
        B = abs(eta) * np.exp(1j*np.angle(A))
        C = ifft(B)
        D = abs(eps) * np.exp(1j*np.angle(C))
        A = fft(D)
    phi = np.unwrap(np.angle(A))
    # Return in r coordinates
    return np.sqrt(chi), phi


def f_to_z(f, lam):
    return 2*np.pi**2/(lam*f)


def z_to_f(z, lam):
    return 2*np.pi**2/(lam*z)
