#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:39:10 2017

@author: rariniello
"""

import numpy as np


def keldysh(EI, I, wavelength):
    """ Calculates the Keldysh parameter.

        Calculates the Keldysh parameter from the ionization energy and laser
        intensity assuming a linearly polarized electric field. Tunneling
        ionization if gamma << 1, multi-photon ionization if gamma > 1.

        Parameters
        ----------
        EI
            Ionization energy of the electron in eV
        I
            Intensity in 10^14 W/cm^2
        wavelength
            Laser wavelength in um

        Returns
        -------
        gamma
            Keldysh parameter
    """
    gamma = np.sqrt(EI / (18.6*I*np.power(wavelength, 2)))
    return gamma


def field_from_intensity(I, n=1.0):
    """ Calculates the electric field from the intensity.

        Calculates the peak electric field from the intensity assuming a
        monochromatic propogating wave.

        Parameters
        ----------
        I
            Intensity in 10^14 W/cm^2.
        n
            Index of refraction of the medium the wave propogates through.

        Returns
        -------
        E
            Electric field in GV/m.
    """
    E = 27.4492 * np.sqrt(I/n)
    return E


def gaussian_envelope(I, t, tau, chirp=0):
    """ Returns the envelope of the electric field of a Gaussian pulse.

    Parameters
    ----------
    I
        Peak intensity of the pulse.
    t
        Array of time points to return the electric field at.
    tau
        RMS length of the pulse.
    chirp
        Frequency chirp of the pulse.

    Returns
    -------
    E
        Array of electric field vlues at times given by array t.
    """
    a = np.pi / (2*tau**2)
    Gamma = a + chirp*1j
    E0 = field_from_intensity(I)
    E = E0 * np.exp(-Gamma * t**2)
    return E


def gaussian_field(I, t, f, tau, chirp=0):
    """ Returns the electric field of a Gaussian pulse.

    Parameters
    ----------
    I
        Peak intensity of the pulse.
    t
        Array of time points to return the electric field at.
    f
        Carrier wave frequency.
    tau
        RMS length of the pulse.
    chirp
        frequency chirp of the pulse.

    Returns
    -------
    E
        Array of electric field vlues at times given by array t.
    """
    w0 = 2*np.pi * f
    E = gaussian_envelope(I, t, tau, chirp) * np.exp(-1j*w0 * t)
    return E
