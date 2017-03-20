#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:39:10 2017

@author: rariniello
"""

import numpy as np


def keldysh(EI, I, wavelength):
    """ Calculates the Keldysh parameter.
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
