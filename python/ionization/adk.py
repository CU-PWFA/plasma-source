#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 14:25:13 2017

@author: rariniello
"""

import numpy as np
from scipy.special import factorial
from scipy.special import gamma
from ionization import ionization


def adk_static(EI, E, Z, l, m):
    """ Calculates the ionization rate of a gas using the ADK model.

    Calculates the tunneling ionization rate of a gas in a constant electric
    field using the ADK model.

    Parameters
    ----------
    EI
        Ionization energy of the electron in eV.
    E
        Electric field strength in GV/m.
    Z
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l
        Orbital quantum number of the electron being ionized.
    m
        Magnetic quantum number of the electron being ionized.

    Returns
    -------
    w
        Ionization rate in 1/fs
    """
    n = 3.68859*Z / np.sqrt(EI)
    E0 = np.power(EI, 3/2)
    Cn2 = (np.power(4, n)) / (n*gamma(2*n))
    N = 1.51927 * (2*l+1) * factorial(l+abs(m)) \
        / (np.power(2, abs(m)) * factorial(abs(m)) * factorial(l-abs(m)))
    w = N * Cn2 * EI \
        * np.power(20.4927*E0/E, 2*n-abs(m)-1) \
        * np.exp(-6.83089*E0/E)
    return w


def adk_linear(EI, E, Z, l, m):
    """ Calculates the ionization rate of a gas using the ADK model.

    Calculates the average tunneling ionization rate of a gas in a linearly
    polarized electric field. Use this function in conjugtion with the envelope
    of the pulse to find the ionization fraction.

    Parameters
    ----------
    EI
        Ionization energy of the electron in eV.
    E
        Electric field strength in GV/m.
    Z
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l
        Orbital quantum number of the electron being ionized.
    m
        Magnetic quantum number of the electron being ionized.

    Returns
    -------
    w
        Ionization rate in 1/fs
    """
    E0 = np.power(EI, 3/2)
    w = 0.305282 * np.sqrt(E/E0) * adk_static(EI, E, Z, l, m)
    return w


def ionization_frac(EI, I, t, Z, l, m):
    """ Work in progress, ionization fraction from intensity and time

    t is in fs, I is 10^14 W/cm^2
    """
    E = ionization.field_from_intensity(I)
    frac = adk_linear(EI, E, Z, l, m) * t
    frac[frac >= 1.0] = 1.0
    return frac
