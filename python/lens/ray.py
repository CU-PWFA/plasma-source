#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:44:21 2017

@author: rariniello
"""

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d


def uniform_phase(I, z, R, r0=0):
    """ Generates a phase mask to create a on axis intensity profile.

    Generates a phase mask to produce a passed longitudinal on axis intensity
    profile assuming a uniform (flat top) mask illumination. The calculation
    uses pure ray tracing to create the phase mask and does not account for
    diffraction effects.

    Parameters
    ----------
    I : array-like
        Array of desired intensity profile along the optical axis.
    z : array-like
        Array of the distances in z from the phase mask.
    R : double
        Radius of the phase mask / input beam.
    r0 : double, optional
        Starting radius of the phase mask, places a hole of radius r0 in the
        middle of the phase mask.

    Returns
    -------
    I0 : double
        The required input intensity to achieve the desired output.
    r : array-like
        Radius vector where the phase is specified at. Each radius refracts
        rays to the corresponding element in z.
    phi : array-like
        Phase of the mask at each r. Note that this phase must be multiplied by
        k to get the actual phase delay that goes in the argument of the
        complex exponential.
    """
    # Set up the arrays to store everything
    Nz = np.size(z)
    r = np.zeros(Nz)
    phi = np.zeros(Nz)
    r[0] = r0
    # Calculate the required input intensity
    I0 = integrate.trapz(I, z) / (np.pi * R**2)
    # Calculate r and phi incrementally
    sinnew = r0/np.sqrt(r0**2 + z[0]**2)
    for i in range(1, Nz):
        Iavg = (I[i] + I[i-1]) / 2
        dz = z[i] - z[i-1]
        r[i] = np.sqrt(Iavg*dz/(np.pi*I0) + r[i-1]**2)
        dr = r[i] - r[i-1]
        sinold = sinnew
        sinnew = r[i] / np.sqrt(r[i]**2 + z[i]**2)
        phi[i] = phi[i-1] - (sinnew + sinold)*dr/2
    # Return everything
    return I0, r, phi


def arbitrary_phase(I0, rin, I, z, r0=0, m=1.0):
    """ Generates a phase mask to create a on axis intensity profile.

    Generates a phase mask to produce a passed longitudinal on axis intensity
    profile with arbitrary mask illumination. The calculation uses pure ray
    tracing to create the phase mask and does not account for diffraction
    effects.

    Parameters
    ----------
    I0 : array-like
        Intensity at the lens, each element corresponds to an element in rin.
    rin : array-like
        Array of radius values the input beam is specified at.
    I : array-like
        Array of desired intensity profile along the optical axis.
    z : array-like
        Array of the distances in z from the phase mask.
    r0 : double, optional
        Starting radius of the phase mask, places a hole of radius r0 in the
        middle of the phase mask.
    m : double, optional
        Correction factor since 2D->1D energy mapping isn't exact, decrease if
        the output doesn't cover the entire lens. If a outside of interpolation
        error is thrown, increase it.

    Returns
    -------
    I0 : double
        The required input intensity to achieve the desired output.
    r : array-like
        Radius vector where the phase is specified at. Each radius refracts
        rays to the corresponding element in z.
    phi : array-like
        Phase of the mask at each r. Note that this phase must be multiplied by
        k to get the actual phase delay that goes in the argument of the
        complex exponential.
    """
    # Set up the arrays to store everything
    Nz = np.size(z)
    r = np.zeros(Nz)
    phi = np.zeros(Nz)
    r[0] = r0
    # Calculate the required input intensity
    Iamp = integrate.trapz(I, z) / integrate.trapz(2*np.pi*I0*rin, rin)
    I0 = I0*Iamp*m
    I0 = interp1d(rin, I0)
    # Calculate r and phi incrementally
    sinnew = r0/np.sqrt(r0**2 + z[0]**2)
    for i in range(1, Nz):
        Iavg = (I[i] + I[i-1]) / 2
        dz = z[i] - z[i-1]
        r[i] = np.sqrt(Iavg*dz/(2*np.pi*I0(r[i-1])) + r[i-1]**2)
        dr = r[i] - r[i-1]
        sinold = sinnew
        sinnew = r[i] / np.sqrt(r[i]**2 + z[i]**2)
        phi[i] = phi[i-1] - (sinnew + sinold)*dr/2
    # Return everything
    return Iamp, r, phi


def lens_design(I0, rin, I, rp, L):
    """ Generates a phase only lens to produce a desired radial intensity.

    Uses ray tracing to produce a phase only lens that will produce a desired
    radial intensity pattern. Note the intensity patterns must contain the same
    amount of power.

    Parameters
    ----------
    I0 : array-like
        Intensity at the lens, each element corresponds to an element in rin.
    rin : array-like
        Array of radius values the input beam is specified at.
    I : array-like
        Array of desired intensity.
    rp : array-like
        Array of radiuses the target intensity is specified at.
    L : double
        Distance between the phase mask and the target.

    Returns
    -------
    r : array-like
        Radius vector where the phase is specified at. Each radius refracts
        rays to the corresponding element in rp.
    phi : array-like
        Phase of the mask at each r. Note that this phase must be multiplied by
        k to get the actual phase delay that goes in the argument of the
        complex exponential.
    """
    N = np.size(rp)
    r = np.zeros(N)
    phi = np.zeros(N)
    I0 = interp1d(rin, I0, bounds_error=False, fill_value=0.0)
    # Calculate r and phi incrementally
    sinnew = rp[0]/np.sqrt(rp[0]**2 + L**2)
    for i in range(1, N):
        dr = rp[i] - rp[i-1]
        # XXX this will break for donut beams, the annulus will be too big
        if I0(r[i-1]) == 0.0:
            r[i] = r[i-1]
            continue
        r[i] = np.sqrt((I[i]*rp[i] + I[i-1]*rp[i-1])*dr/I0(r[i-1]) + r[i-1]**2)
        dr = r[i] - r[i-1]
        sinold = sinnew
        sinnew = (r[i] - rp[i]) / np.sqrt((r[i]-rp[i])**2 + L**2)
        phi[i] = phi[i-1] - (sinnew + sinold)*dr/2
    return r, phi
