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
from ionization import ionization


def rate_static(EI, E, Z, l, m):
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
    # Cast E as an numpy array, handles imputs passed as doubles
    E = np.array(E)
    w = np.zeros(np.size(E))
    n = 3.68859*Z / np.sqrt(EI)
    E0 = np.power(EI, 3/2)
    Cn2 = (np.power(4, n)) / (n*gamma(2*n))
    N = 1.51927 * (2*l+1) * factorial(l+abs(m)) \
        / (np.power(2, abs(m)) * factorial(abs(m)) * factorial(l-abs(m)))
    # Only calculate the ionization rate when z is nonzero
    Enz = E[E > 0]
    if np.size(Enz) > 0:
        w[E > 0] = N * Cn2 * EI \
            * np.power(20.4927*E0/Enz, 2*n-abs(m)-1) \
            * np.exp(-6.83089*E0/Enz)
    return w


def rate_linear(EI, E, Z, l, m):
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
    w = 0.305282 * np.sqrt(E/E0) * rate_static(EI, E, Z, l, m)
    return w


def ionization_frac(EI, E, t, Z, l, m, envelope=True):
    """ Calculates the ionization fraction from a time varying electric field.

    Calculates the ionization fraction for a time varying electric field. This
    function can take either the explicit electric field or the envelope of the
    electric field (ignoring the fast oscillation).

    Parameters
    ----------
    EI
        Ionization energy of the electron in eV.
    E
        Electric field strength in GV/m. Array of E at time t.
    t
        Array of time values, in fs, matching the electric field vector.
    Z
        Atomic residue i.e. which electron is being ionizaed (1st, 2nd...).
    l
        Orbital quantum number of the electron being ionized.
    m
        Magnetic quantum number of the electron being ionized.
    envelope
        Specifies whether the envelope of the electric field is passed in or
        the explicit electric field is passed. Controls whether the function
        uses the time averaged ADK model or the static ADK model.
    """
    if envelope:
        rate = rate_linear(EI, E, Z, l, m)
    else:
        rate = rate_static(EI, E, Z, l, m)
    frac = 1-np.exp(-integrate.simps(rate, t))
    return frac
