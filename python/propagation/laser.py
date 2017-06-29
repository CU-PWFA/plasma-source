#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 21:28:31 2017

@author: rariniello
"""

from propagation import C
import numpy as np
from numpy.lib.scimath import sqrt
from numpy.fft import fft, ifft, fft2, ifft2, fftfreq
from scipy import integrate


def fourier_prop(E, x, z, lam, n=1):
    """ Propogates an electromagnetic wave from a 1D boundary.

    Uses the Rayleigh-Sommerfeld transfer function to propagate an
    electromagnetic wave from a 1D boundary. The calculation assumes a
    homogeneous index of refraction in the region of propagation. This function
    is not equivalent to taking a slice through the full 2D boundary
    calculation. For example, it will not give the correct result for the
    amplitude of a Gaussian beam propagating in z.

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
    # Need these to calculate frequencies later on
    Nx = np.size(x)
    Nz = np.size(z)
    X = x[Nx-1]-x[0]
    dx = X / (Nx-1)
    # Fourier transform the electric field on the boundary
    eb = fft(E)
    # Calculate the transfer function at every point in Fourier space
    fx = fftfreq(Nx, dx)
    fz = sqrt((n/lam)**2 - fx**2)
    # Reshape to create the Nz X Nx output array
    z = np.reshape(z, (Nz, 1))
    e = np.exp(1j*2*np.pi*z*fz) * eb
    # Inverse Fourier transform to get the real-space field
    e = ifft(e, axis=1)
    return e


def scale_fourier(e):
    """ Scales the output of fourier_prop to approximate the real 2D field.

    Parameters
    ----------
    e : array_like
        Array of E field values from fourier_prop.

    Returns
    -------
    e : array_like
        Scaled electric field values.
    """
    s = np.amax(abs(e), axis=1)
    return e * np.reshape(s, (np.size(s), 1))


def fresnel_prop(E, x, z, lam, n=1):
    """ Propogates an electromagnetic wave from a 1D boundary.

    Uses the Fresnel transfer function to propagate an electromagnetic wave
    from a 1D boundary. The calculation assumes a homogeneous index of
    refraction in the region of propagation. This function is equivalent to
    taking a slice through the full 2D boundary calculation. However, it is
    only valid in the paraxial approximation. The beam must be cnetered in the
    x array in order for this function to work correctly.

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
    # Need these to calculate frequencies later on
    Nx = np.size(x)
    Nz = np.size(z)
    X = x[Nx-1]-x[0]
    dx = X / (Nx-1)
    # Fourier transform the electric field on the boundary
    eb = fft(E)
    # Reshape to create the Nz X Nx output array
    z = np.reshape(z, (Nz, 1))
    # Calculate the transfer function at every point in Fourier space
    fx = fftfreq(Nx, dx)
    H = np.exp(-1j*np.pi*n*z*fx**2/lam)
    e1 = H * eb
    # Inverse Fourier transform to get the real-space field
    e1 = ifft(e1, axis=1)
    # Multiply by the overal phase shift and handle y component
    e = np.exp(1j*2*np.pi*n*z/lam) * e1 * np.reshape(e1[:, int(Nx/2)], (Nz, 1))
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
    # Reshape to create Nz X Nx X Ny output array
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
        T = np.exp(1j * np.pi * (nih[:, i-1]+nih[:, i]) * dz / lam)
    return e


def beam_prop2(E, nih, x, y, z, lam, nh):
    """ Propogates an electromagnetic wave from a 2D boundary.

    Propogates a beam from a 2D boundary through a region with an inhomogenous
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
        Array of E field values at points (x, y) along the boundary.
    nih : array_like
        Variation in index of refraction from nh passed as (x, y, z).
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
    nh : double
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : array_like
        The electric field at position (z, x, y).
    """
    Nx = np.size(x)
    Ny = np.size(y)
    Nz = np.size(z)
    e = np.zeros((Nz, Nx, Ny), dtype=np.complex)
    # Fourier propagate to the first z grid point
    e[0, :, :] = fourier_prop2(E, x, y, z[0], lam, nh)
    T = 1
    for i in range(1, Nz):
        dz = z[i] - z[i-1]
        # Phase transmission function for refraction
        T = np.exp(1j * np.pi * (nih[:, :, i-1]+nih[:, :, i]) * dz / lam)
        # Fourier propogate to the next z point
        e[i] = fourier_prop2(e[i-1], x, y, [dz], lam, nh)
        e[i] = e[i] * T
    return e


def pulse_prop(E, Et, x, z, t, c=C, n=1):
    """ Propogates an electromagnetic wave pulse from a 1D boundary.

    Uses the Rayleigh-Sommerfeld transfer function to propagate each temporal
    frequency componenet of an electromagnetic wave pulse from a 1D boundary.
    The calculation assumes a homogeneous index of refraction in the region of
    propagation. Additionally it assumes the temporal pulse shape is the same
    throughout the spatial extent of the pulse.

    Parameters
    ----------
    E : array_like
        Array of E field values at points x along the boundary.
    Et : array_like
        Array of E field values at points t in time.
    x : array_like
        Array of transverse locations in x on the boundary. Elements must be
        evenly spaced for fft.
    z : array_like
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    t : array_like
        Array of t times the field is specified at. Elements must be evenly
        spaced in time for fft.
    c : double, optional
        Pass c in the units of time and lam to ensure correct calculation.
    n : double, optional
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : array_like
        The electric field for a given frequency at all x and z positions.
        e(z, x, ft)
    """
    # Calculate time frequencies from t array
    Nx = np.size(x)
    Nz = np.size(z)
    Nt = np.size(t)
    T = t[Nt-1]-t[0]
    dt = T / (Nt-1)
    ft = fftfreq(Nt, dt)
    # 2/Nt removes the normalization built into the FFT
    a = 2*fft(Et)/Nt
    # Array to store field from each frequency
    e = np.zeros((Nz, Nx, Nt), dtype=np.complex)
    # Propogate each frequency
    for i in range(Nt):
        if ft[i] > 0:
            l = c/ft[i]
            ea = E*a[i]
            e[:, :, i] = fourier_prop(ea, x, z, l, n)
        else:
            e[:, :, i] = np.zeros((Nz, Nx), dtype=np.complex)
    return e


def pulse_prop2(E, Et, x, y, z, t, c=C, n=1):
    """ Propogates an electromagnetic wave pulse from a 2D boundary.

    Uses the Rayleigh-Sommerfeld transfer function to propagate each temporal
    frequency componenet of an electromagnetic wave pulse from a 2D boundary.
    The calculation assumes a homogeneous index of refraction in the region of
    propagation. Additionally it assumes the temporal pulse shape is the same
    throughout the spatial extent of the pulse.

    Currently returns the field on a plane at y[Ny/2]. We don't return the full
    field for memory reasons and you can only animate a 2D slice anyways.

    Parameters
    ----------
    E : array_like
        Array of E field values at points x along the boundary.
    Et : array_like
        Array of E field values at points t in time.
    x : array_like
        Array of transverse locations in x on the boundary. Elements must be
        evenly spaced for fft.
    y : array_like
        Array of transverse locations in y on the boundary. Elements must be
        evenly spaced for fft.
    z : array_like
        Array of z distances from the boundary to calculate the field at. Does
        not need to be evenly spaced.
    t : array_like
        Array of t times the field is specified at. Elements must be evenly
        spaced in time for fft.
    c : double, optional
        Pass c in the units of time and lam to ensure correct calculation.
    n : double, optional
        Index of refraction of the medium the wave is propagating through.

    Returns
    -------
    e : array_like
        The electric field for a given frequency at all x and z positions.
        e(z, x, ft)
    """
    # Calculate time frequencies from t array
    Nx = np.size(x)
    Ny = np.size(y)
    Nz = np.size(z)
    Nt = np.size(t)
    T = t[Nt-1]-t[0]
    dt = T / (Nt-1)
    ft = fftfreq(Nt, dt)
    # 2/Nt removes the normalization built into the FFT
    a = 2*fft(Et)/Nt
    # Array to store field from each frequency
    e = np.zeros((Nz, Nx, Nt), dtype=np.complex)
    # Propogate each frequency
    for i in range(Nt):
        print('Frequency component ', i+1, ' of ', Nt)
        if ft[i] > 0:
            l = c/ft[i]
            ea = E*a[i]
            e[:, :, i] = fourier_prop2(ea, x, y, z, l, n)[:, :, int(Ny/2)]
        else:
            e[:, :, i] = np.zeros((Nz, Nx), dtype=np.complex)
    return e


def pulse_time(e, z, t, time, c=C, n=1):
    """ Returns the electric field at a given time for a pulse.

    Returns the electric field at a given time at places in the x-z plane. Uses
    the output of pulse_prop and the same t array passed to pulse_prop. Returns
    the electric field at a given time. General usage would be:

    e = pulse_prop(E, Et, x, z, t)
    pulse_time(e, z, t, time)

    Parameters
    ----------
    e : array_like
        Array returned by pulse_prop, electric field written as e(ft, z, x).
    z : array_like
        z array passed to pulse_prop.
    t : array_like
        t array passed to pulse_prop, ft=fftfreq(Nt, dt). Where Nt is the
        number of elements in t and dt is the spacing of elements in t.
    time : array_like
        Time to return the electric field at.

    Returns
    -------
    et : array_like
        The electric field at time, e(time, z, x).
    """
    # Cast time as a numpy array, handles imputs passed as naked doubles
    time = np.array(time, ndmin=1)
    Ntime = np.size(time)
    time = np.reshape(time, (Ntime, 1))
    # Get time frequencies
    Nt = np.size(t)
    T = t[Nt-1]-t[0]
    dt = T / (Nt-1)
    ft = fftfreq(Nt, dt)
    Nz = np.size(z)
    Z = z[Nz-1]-z[0]
    # Calculate complex expontential in time
    exp = np.exp(-1j*2*np.pi*(time-t[0])*ft)
    exp = np.reshape(exp, (Ntime, 1, Nt, 1))
    # matmul uses last two indexes of each argument returns a Ntime X Nz X Nx
    et = np.matmul(e, exp)
    et = np.squeeze(et, axis=3)
    # Spatial size of time window T, must be integer multiple of c*T
    # The +2 ensures the pulse fully leaves the screen before the next pulse
    dz = c*T/(2*n)
    size = (int(Z/(c*T)) + 2)*c*T

    def loop_variable(var, size):
        """places a variable, var, in the range 0-size regardless of sign."""
        return ((int(var/size)+1)*size + var) % size

    # Set anything not in the time frame to 0
    for i in range(Ntime):
        # Center, left, and right sides of time frame
        zc = c*time[i]/n
        zl = loop_variable(zc - dz, size) + z[0]
        zr = loop_variable(zc + dz, size) + z[0]
        # Check to make sure the left side hasn't looped around
        if (zr > zl):
            sel = np.array(z < zl) + np.array(z > zr)
        elif (zr < zl):
            sel = np.array(z < zl) * np.array(z > zr)
        et[i, sel, :] = 0
    return et


def fresnel_axis(E, r, z, lam, n=1):
    """ Returns the electric field along the optical axis in the Fresnel limit.

    Uses the Fresnel diffraction integral to calculate the electric field along
    the optical axis resulting from a cylindrically symmetric felectric field.
    Not this is only valid in the paraxial limit. 

    Parameters
    ----------
    E : array_like
        Array of E field values at points r on the boundary, E(r).
    r : array_like
        Array of radial locations in on the boundary. Doesn't need to be evenly
        spaced.
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
        The electric field at position (z).
    """
    Nz = np.size(z)
    k = 2*np.pi*n/lam
    # Reshape to create the Nz X Nr array for integration
    z = np.reshape(z, (Nz, 1))
    # Integral prefactor
    pre = k*np.exp(1j*k*z) / (1j*z)
    # Calculate the argument of the integral
    arg = E * np.exp(1j*k*r**2/(2*z)) * r  
    e = pre * integrate.simps(arg, r, axis=1)
    return e
