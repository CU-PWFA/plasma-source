#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 14:25:13 2017

@author: rariniello
"""

import numpy as np
from scipy.special import factorial
from scipy.special import gamma
from scipy import integrate
# from ionization import ionization


def rate_static(EI, E, Z, l=0, m=0):
    """ Calculates the ionization rate of a gas using the ADK model.

    Calculates the tunneling ionization rate of a gas in a constant electric
    field using the ADK model.

    Parameters
    ----------
    EI : double
        Ionization energy of the electron in eV.
    E : array_like
        Electric field strength in GV/m.
    Z : int
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l : int, optional
        Orbital quantum number of the electron being ionized. Defaults to 0.
    m : int, optional
        Magnetic quantum number of the electron being ionized. Defaults to 0.

    Returns
    -------
    w : array_like
        Ionization rate in 1/fs.
    """
    # Cast E as an numpy array, handles imputs passed as naked doubles
    E = np.array(E, ndmin=1)
    w = np.zeros(np.size(E))
    n = 3.68859*Z / np.sqrt(EI)
    E0 = np.power(EI, 3/2)
    Cn2 = (np.power(4, n)) / (n*gamma(2*n))
    N = 1.51927 * (2*l+1) * factorial(l+abs(m)) \
        / (np.power(2, abs(m)) * factorial(abs(m)) * factorial(l-abs(m)))
    # Only calculate the ionization rate when z is nonzero
    Enz = np.array(E > 0)
    if np.size(Enz) > 0:
        w[Enz] = N * Cn2 * EI \
            * np.power(20.4927*E0/E[Enz], 2*n-abs(m)-1) \
            * np.exp(-6.83089*E0/E[Enz])
    return w


def rate_linear(EI, E, Z, l=0, m=0):
    """ Calculates the ionization rate of a gas using the ADK model.

    Calculates the average tunneling ionization rate of a gas in a linearly
    polarized electric field. Use this function in conjunction with the
    envelope of the pulse to find the ionization fraction.

    Parameters
    ----------
    EI : double
        Ionization energy of the electron in eV.
    E : array_like
        Electric field strength in GV/m.
    Z : int
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l : int, optional
        Orbital quantum number of the electron being ionized. Defaults to 0.
    m : int, optional
        Magnetic quantum number of the electron being ionized. Defaults to 0.

    Returns
    -------
    w : array_like
        Ionization rate in 1/fs
    """
    w = 0.305282 * np.sqrt(E/np.power(EI, 3/2)) * rate_static(EI, E, Z, l, m)
    return w


def ionization_frac(EI, E, t, Z, l=0, m=0, envelope=True):
    """ Calculates the ionization fraction from a time varying electric field.

    Calculates the ionization fraction for a time varying electric field. This
    function can take either the explicit electric field or the envelope of the
    electric field (ignoring the fast oscillation).

    Parameters
    ----------
    EI : double
        Ionization energy of the electron in eV.
    E : array_like
        Electric field strength in GV/m. Array of E at time t.
    t : array_like
        Array of time values, in fs, matching the electric field vector.
    Z : int
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l : int, optional
        Orbital quantum number of the electron being ionized. Defaults to 0.
    m : int, optional
        Magnetic quantum number of the electron being ionized. Defaults to 0.
    envelope : bool, optional
        Specifies whether the envelope of the electric field is passed in or
        the explicit electric field is passed. Controls whether the function
        uses the time averaged ADK model or the static ADK model.

    Returns
    -------
    frac : double
        Fraction of gas atoms ionized.
    """
    if envelope:
        rate = rate_linear(EI, E, Z, l, m)
    else:
        rate = rate_static(EI, E, Z, l, m)
    frac = 1 - np.exp(-integrate.simps(rate, t))
    return frac


def gaussian_frac(EI, E, tau, Z, l=0, m=0):
    """ Calculates the ionization fraction from a Gaussian laser pulse.

    Parameters
    ----------
    EI : double
        Ionization energy of the electron in eV.
    E : array_like
        Peak electric field strength in GV/m. This is the peak of the envelope,
        not necessairly the true peak field, it is the peak field assuming
        0 carrier envelope offset.
    tau : array_like
        RMS length of the pulse in fs.
    Z : int
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l : int, optional
        Orbital quantum number of the electron being ionized. Defaults to 0.
    m : int, optional
        Magnetic quantum number of the electron being ionized. Defaults to 0.

    Returns
    -------
    frac : array_like
        Fraction of gas atoms ionized.
    """
    rate = rate_static(EI, E, Z, l, m)
    frac = 1 - np.exp(-6.82909 * E * tau * rate / np.power(EI, 3/2))
    return frac
