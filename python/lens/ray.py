#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:44:21 2017

@author: rariniello
"""

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
from ionization import ionization


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
        the output doesn't cover the entire lens. 

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
    I0 += np.amax(I0)*0.0001 # we need some value everywhere
    I += np.amax(I)*0.0001
    I0 = interp1d(rin, I0, bounds_error=False, fill_value=0.0)
    # Calculate r and phi incrementally
    if r0 == 0:
        sinnew = 0.0
    else:
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


def super_gaussian_phase(params, I, z, ind=None):
    """ Calculates the phase required to create an on-axis intensity profile.
    
    Calculates the phase required to focus a super Gaussian to a specified on-
    axis intensity profile.
    
    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            lam : double
                Wavelength of the electromagnetic wave in vacuum.
            rmax : array-like
                Maximum radius to return the electric field at. Pass an array
                of length L, one for each Bessel mode. Must be a factor of root
                2 larger than the grid in prop.
            rc : array-like
                Array of center locations for the super Gaussians.
            r0 : array-like
                Beginning radius of the higher order Beams.
            w : array-like
                Array of widths for the super Gaussians.
            nGauss : array-like
                Order of the super Gaussians.
            m : array-like
                Multiplier for the arbitrary phase function, test before use.
    z : array-like
        Array of on axis z values the intensity is specified at.
    I : array-like
        Desired Intensity profile along the optical axis.
    ind : int, optional
        The index of the parameter arrays to use. Default is 0.
        
    Returns
    -------
    E : array-like
        The electric field with phase after the lens.
    r : array-like
        The radiuses the electric field is returned at.
    """
    if ind is None:
        ind = 0
    rmax = params['rmax'][ind]
    rc = params['rc'][ind]
    r0 = params['r0'][ind]
    w = params['w'][ind]
    m = params['m'][ind]
    n = params['nGauss'][ind]
    k = 2*np.pi/params['lam']
    # Create the super-Gaussian
    rin = np.linspace(r0, rmax, 10000)
    I0 = np.exp(-2*((rin-rc)/w)**n)
    Iamp, r, phi = arbitrary_phase(I0, rin, I, z, r0=r0, m=m)
    E = Iamp * np.exp(-2*((r-rc)/w)**n)
    E = ionization.field_from_intensity(E) * np.exp(1j*k*phi)
    return E, r


def lens_design(Iin, rin, Iout, rout, L):
    """ Generates a phase only lens to produce a desired radial intensity.

    Uses ray tracing to produce a phase only lens that will produce a desired
    radial intensity pattern. Note the intensity patterns must contain the same
    amount of power.

    Parameters
    ----------
    Iin : array-like
        Intensity at the lens, each element corresponds to an element in rin.
    rin : array-like
        Array of radius values the input beam is specified at.
    Iout : array-like
        Array of desired intensity.
    rout : array-like
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
    N = np.size(rout)
    r = np.zeros(N)
    phi = np.zeros(N)
    Iin = interp1d(rin, Iin, bounds_error=False, fill_value='extrapolate')
    # Calculate r and phi incrementally
    sinnew = rout[0]/np.sqrt(rout[0]**2 + L**2)
    for i in range(1, N):
        dr = rout[i] - rout[i-1]
        # XXX this will break for donut beams, the annulus will be too big
        if Iin(r[i-1]) == 0.0:
            r[i] = r[i-1]
            continue
        r[i] = np.sqrt((Iout[i]*rout[i] + Iout[i-1]*rout[i-1])*dr/Iin(r[i-1]) + r[i-1]**2)
        dr = r[i] - r[i-1]
        sinold = sinnew
        sinnew = (r[i] - rout[i]) / np.sqrt((r[i]-rout[i])**2 + L**2)
        phi[i] = phi[i-1] - (sinnew + sinold)*dr/2
    return r, phi
