#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 08:49:19 2017

@author: robert
"""

import numpy as np
from ionization import ionization
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.special import erf


def plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, N, zf=0):
    """ Creates a plasma density with Gaussian ends.

    Returns a plasma profile with a flattop and Gaussian ramps on either side.

    Parameters
    ----------
    z0 : double
        The distance at which the uniform fully ionized plasma starts.
    dz : double
        The length of the fully ionized plasma.
    sigmaIn : double
        The length of the input ramp, sigma for the Gaussian.
    sigmaOut : double
        The length of the output ramp, sigma for the Gaussian.
    N : int
        Number of grid points in z.
    zf : double, optional
        The end of the grid in z.

    Returns
    -------
    z : array-like
        Array of distances along the optical axis the intensity is returned at.
    I : array-like
        Intensity profile given at each point in z.
    """
    # Create the z grid
    d1 = 0.0
    d2 = z0
    d3 = z0 + dz
    if zf is not 0:
        d4 = zf
    else:
        d4 = d3+10*sigmaOut
    Z = d4 - d1
    z = np.linspace(d1, d1+Z, N)
    # Create the density profile
    frac = np.zeros(N)
    # 0.999 prevents going outside of the interpolating functions range
    peak = 0.999
    sel = z <= d2
    frac[sel] = peak*np.exp(-(z[sel]-d2)**2/(2*sigmaIn**2))
    sel = np.array(z > d2) * np.array(z < d3)
    frac[sel] = peak
    sel = z >= d3
    frac[sel] = peak*np.exp(-(z[sel]-d3)**2/(2*sigmaOut**2))
    return z, frac


def plasma_linear_ramps(z0, dz, sigmaIn, sigmaOut, N, zf=0):
    """ Creates a plasma density with Gaussian ends.

    Returns a plasma profile with a flattop and Gaussian ramps on either side.

    Parameters
    ----------
    z0 : double
        The distance at which the uniform fully ionized plasma starts.
    dz : double
        The length of the fully ionized plasma.
    sigmaIn : double
        The length of the input ramp, halfwidth for linear.
    sigmaOut : double
        The length of the output ramp, halfwidth for linear.
    N : int
        Number of grid points in z.
    zf : double, optional
        The end of the grid in z.

    Returns
    -------
    z : array-like
        Array of distances along the optical axis the intensity is returned at.
    I : array-like
        Intensity profile given at each point in z.
    """
    # Create the z grid
    d1 = 0.0
    d2 = z0
    d3 = z0 + dz
    if zf is not 0:
        d4 = zf
    else:
        d4 = d3+2*sigmaOut
    Z = d4 - d1
    z = np.linspace(d1, d1+Z, N)
    # Create the density profile
    frac = np.zeros(N)
    # 0.999 prevents going outside of the interpolating functions range
    peak = 0.999
    sel = np.array(z <= d2) * np.array(z >= d2-2*sigmaIn)
    frac[sel] = peak*(z[sel]-d2)/(2*sigmaIn) + peak
    sel = np.array(z > d2) * np.array(z < d3)
    frac[sel] = peak
    sel = np.array(z >= d3) * np.array(z <= d3+2*sigmaOut)
    frac[sel] = peak*(-(z[sel]-d3)/(2*sigmaOut)) + peak
    return z, frac


def intensity_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, ion, N, zf=0):
    """ Creates an intensity profile for a uniform plasma with Gaussian ends.

    Returns the on-axis intensity profile necessary to create a fully ionized
    on-axis plasma with Gaussian entrance and exit density ramps. The
    calculation assumes a 

    Parameters
    ----------
    z0 : double
        The distance at which the uniform fully ionized plasma starts.
    dz : double
        The length of the fully ionized plasma.
    sigmaIn : double
        The length of the input ramp, sigma for the Gaussian.
    sigmaOut : double
        The length of the output ramp, sigma for the Gaussian.
    ion : dictionary
        atom : dictionary
            See the description in ionization.ionization for details.
        tau : double
            Temporal length of the pulse in fs.
        type : string
            The type of pulse, i.e. gaussian, flatend, etc.
    N : int
        Number of grid points in z.
    zf : double, optional
        The end of the grid in z.

    Returns
    -------
    z : array-like
        Array of distances along the optical axis the intensity is returned at.
    I : array-like
        Intensity profile given at each point in z.
    """
    z, frac = plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, N, zf)
    # Get the required intensity
    I = ionization.intensity_from_density(ion, frac)
    return z, I


def smoothed_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, ion, N, order, alpha, 
                            zf=0):
    """ This function smooths out the hard corners of the Gaussian ramps.
    
    First, this function corrects for the DC baseline due to the slow asymptote
    of the inverse ADK model. Second it domes the flatend intensity region to
    the sharp corners at the transition to 100% ionization fraction.
    
    Parameters
    ----------
    z0 : double
        The distance at which the uniform fully ionized plasma starts.
    dz : double
        The length of the fully ionized plasma.
    sigmaIn : double
        The length of the input ramp, sigma for the Gaussian.
    sigmaOut : double
        The length of the output ramp, sigma for the Gaussian.
    ion : dictionary
        atom : dictionary
            See the description in ionization.ionization for details.
        tau : double
            Temporal length of the pulse in fs.
        type : string
            The type of pulse, i.e. gaussian, flatend, etc.
    N : int
        Number of grid points in z.
    order : int
        The order of the super Gaussian flattop.
    alpha : double
        The perctentage higher the super-Gaussian is than the flattop.

    Returns
    -------
    z : array-like
        Array of distances along the optical axis the intensity is returned at.
    I : array-like
        Intensity profile given at each point in z.
    """
    # First create the unsmoothed ramp
    z, I = intensity_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, ion, N, zf)
    # Remove the dc by fitting a parabola below 5% of peak intensity
    Imax = np.amax(I)
    smallI = np.zeros(N)
    smallI[I > Imax*0.05] = 1.0
    deltaZ = (zf-1)/N
    grad = np.gradient(I, deltaZ)
    # Find where the non-small region begins and ends
    smallInd = np.nonzero(smallI)
    beg = np.amin(smallInd)
    end = np.amax(smallInd)
    # Find the parabola for the beginning
    begI1 = I[beg]
    begzp = grad[beg]
    begA = begzp**2 / (4*begI1)
    begB = z[beg] - 2*begI1/begzp
    begPar = begA * (z[:beg]-begB)**2
    sel = z[:beg] < begB
    begPar[sel] = 0.0
    # Find the parabola for the end
    endI1 = I[end]
    endzp = grad[end]
    endA = endzp**2/(4*endI1)
    endB = z[end] - 2*endI1/endzp
    endPar = endA*(z[end:]-endB)**2
    sel = z[end:] > endB
    endPar[sel] = 0.0
    # Create a smoothing curve for the center.
    sel = I > 0.999*Imax
    zc = 0.5*(z[sel][-1] + z[sel][0])
    w = (z[sel][0] - zc)**(2*order)/0.693
    I[sel] += alpha*Imax*(np.exp(-(z[sel] - zc)**(2*order)/w)-0.5)
    # Update the intensity with the new curves
    I[:beg] = begPar
    I[end:] = endPar
    return z, I


def cutoff_gaussian_intensity(Nz, Z, z0, length, l_ent, s_ent, l_ext, s_ext, ion,
                              order, alpha, path, xlim=None, z_start=0.0):
    """ Calculate the required intensity for a plasma with Gaussian ramps.
    
    Creates a plasma density with Gaussian ramps that end at a hard cutoff.
    Calculates the on-axis intensity required to ionize such a distribution
    and fits parabolas to the hard cutoffs to smooth out the intensity. Adds a
    super Gaussian intensity to the flattop part of the profile.
    
    
    Parameters
    ----------
    Nz : int
        The number of points in the z grid.
    Z : double
        The extent of the z grid, needs to be large for FFT resolution.
    z0 : double
        The distance at which the uniform fully ionized plasma starts.
    length : double
        The length of the fully ionized plasma.
    l_ent : double
        The distance at which the entrance ramp is set to 0.
    s_ent : double
        The length of the entrance ramp, sigma for the Gaussian.
    l_ext : double
        The distance at which the exit ramp is set to 0.
    s_ext : double
        The length of the exit ramp, sigma for the Gaussian.
    ion : dictionary
        atom : dictionary
            See the description in ionization.ionization for details.
        tau : double
            Temporal length of the pulse in fs.
        type : string
            The type of pulse, i.e. gaussian, flatend, etc.
    order : int
        The order of the super Gaussian flattop.
    alpha : double
        The perctentage higher the super-Gaussian is than the flattop.
    path : string
        Path to save the plasma parameters in.
    xlim : optional, tuple or array
        Bounds of z-axis plots.
    z_start : optional, double
        Start of the z grid, defaults to 0.

    Returns
    -------
    z : array-like
        Array of distances along the optical axis the intensity is returned at.
    I : array-like
        Intensity profile given at each point in z.
    """
    # Calculate the plasma density profile
    #--------------------------------------------------------------------------
    z, dz = np.linspace(z_start, Z, Nz, retstep=True)
    frac_l = np.zeros(Nz, dtype='double')
    
    # Uniform accelerating plasma
    sel_u = np.logical_and(z > z0, z < z0+length)
    frac_l[sel_u] = 1.0
    
    # Entrance ramp
    sel_ent = np.logical_and(z >= z0-l_ent, z <= z0)
    ramp_ent = np.exp(-(z-z0)**2/(2*s_ent**2))
    frac_l[sel_ent] = ramp_ent[sel_ent]
    
    #Exit ramp
    sel_ext = np.logical_and(z >= z0+length, z <= z0+length+l_ext)
    ramp_ext = np.exp(-(z-z0-length)**2/(2*s_ext**2))
    frac_l[sel_ext] = ramp_ext[sel_ext]
    
    
    plt.figure(figsize=(8, 2), dpi=150)
    plt.plot(z/1e6, frac_l)
    plt.xlabel(r'z (m)')
    plt.ylabel(r'$n_e/n_{e0}$')
    if xlim is not None:
        plt.xlim(xlim)
    plt.show()
    
    np.save(path+'plasma.npy', [z, frac_l])
    
    # Calculate the required intensity profile
    #--------------------------------------------------------------------------
    # Calculate the desired on axis intensity profile from the ionization fraction
    I_l = ionization.intensity_from_density(ion, 0.999*frac_l)
    
    plt.figure(figsize=(8, 2), dpi=150)
    plt.plot(z/1e6, I_l)
    plt.xlabel(r'z (m)')
    plt.ylabel(r'I ($\mathrm{10^{14}W/cm^2}$)')
    if xlim is not None:
        plt.xlim(xlim)
    plt.show()
    
    # Smooth out the intensity profile
    #--------------------------------------------------------------------------
    # The intensity profile needs to be smoothed to remove steps and hard corners

    # Smooth the center section and increase the intensity
    Imax = np.amax(I_l)
    z_u = z[sel_u]
    zc = 0.5*(z_u[-1] + z_u[0])
    w = (z_u[0] - zc)**(2*order)/0.6
    I_la = np.copy(I_l)
    I_new = alpha*Imax*(np.exp(-(z_u - zc)**(2*order)/w)-0.5)
    I_la[sel_u] += I_new
    
    #Smooth the beginning of the entrance ramp with a parabola
    I_ent = I_l[sel_ent]
    gradI = (I_ent[1] - I_ent[0]) / dz
    A = gradI**2 / (4*I_ent[0])
    B = z[sel_ent][0] - 2*I_ent[0]/gradI
    sel = np.logical_and(z >= B, z < z0-l_ent)
    I_la[sel] = A*(z[sel]-B)**2
    I_la[z < B] = 0.0
    plasma_start = B
    
    #Smooth the end of the exit ramp with a parabola
    I_ext = I_l[sel_ext]
    gradI = (I_ext[-1] - I_ext[-2]) / dz
    A = gradI**2 / (4*I_ext[-1])
    B = z[sel_ext][-1] - 2*I_ext[-1]/gradI
    sel = np.logical_and(z <= B, z > z0+length+l_ext)
    I_la[sel] = A*(z[sel]-B)**2
    I_la[z > B] = 0.0
    plasma_end = B
    
    print('Plasma starts at %0.2fm and ends at %0.2fm' % (plasma_start/1e6, plasma_end/1e6))
    plasma_size = [plasma_start, plasma_end]
    np.save(path+'plasma_size.npy', plasma_size)
    
    plt.figure(figsize=(8, 2), dpi=150)
    plt.plot(z/1e6, I_l, '--m')
    plt.plot(z/1e6, I_la)
    plt.legend(['Without smoothing', 'With smoothing'])
    plt.xlabel(r'z (m)')
    plt.ylabel(r'I ($\mathrm{10^{14}W/cm^2}$)')
    if xlim is not None:
        plt.xlim(xlim)
    plt.show()
    
    return z, I_la


def lithium_oven_profile(z, center, ne0):
    """ Return an interpolating function for the lithium oven gas density.
    
    Parameters
    ----------
    z : array of doubles
        The z array to claculate the density at.
    center : double
        The location for the center of the plasma density profile.
    ne0 : double
        Bulk plasma density.
    
    Returns
    -------
    start : double
        The start of the gas density.
    n : density
        The gas density in the lithium oven at the passed z positions
    n_func : function
        Interpolating function for the gas density as a function of z.
    """
    Nz = len(z)
    n = np.zeros(Nz, dtype='double')
    # Error function ramps?
    # Uniform accelerating plasma
    length = 30.8e4
    z0 = -length/2+center # start of the uniform plasma
    sel_u = np.logical_and(z > z0, z < z0+length)
    n[sel_u] = 1.0
    
    # Entrance ramp
    l_ent = 40e4-length/2 # length of the entrance ramp
    s_ent = 4.2e4
    sel_ent = np.logical_and(z >= z0-l_ent, z <= z0)
    ramp_ent = 0.5*(1+erf((z-z0+13e4)/(np.sqrt(2)*s_ent)))
    n[sel_ent] = ramp_ent[sel_ent]
    
    #Exit ramp
    l_ext = 40e4-length/2 # length of the entrance ramp
    s_ext = 4.2e4
    sel_ext = np.logical_and(z >= z0+length, z <= z0+length+l_ext)
    ramp_ext = 0.5*(1+erf(-(z-z0-13e4-length)/(np.sqrt(2)*s_ext)))
    n[sel_ext] = ramp_ext[sel_ext]
    # Find the start for the refraction code
    start = z0-l_ent
    n *= ne0
    n_func = interp1d(z, n)
    return start, n, n_func
    
    
