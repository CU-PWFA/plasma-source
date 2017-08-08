#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:05:19 2017

@author: robert
"""

from vsim import C
import numpy as np


def get_gamma(data):
    """ Calculates the relativistic factor from the data object for the beam.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    gammaAvg : double
        The average gamma of all the particles in the beam.
    """
    ux = get_ux(data)
    uy = get_uy(data)
    weights = get_weights(data)
    gamma = np.sqrt(1 + (ux**2 + uy**2)/C**2)
    gammaAvg = np.average(gamma, weights=weights)
    return gammaAvg


def get_emittance(data):
    """ Calculates the emittance of the beam from the data object.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    e : double
        The emittance of the beam in mm*mrad.
    """
    y = get_y(data)
    ux = get_ux(data)
    uy = get_uy(data)
    weights = get_weights(data)
    yp = uy / ux
    dy = y - np.average(y, weights=weights)
    dyp = yp - np.average(yp, weights=weights)
    # Calculate the RMS sizes and the correlation
    sigmay2 = np.average(dy**2, weights=weights)
    sigmayp2 = np.average(dyp**2, weights=weights)
    sigmayyp = np.average(dy*dyp, weights=weights)
    # Calculate the emittance
    e = np.sqrt(sigmay2*sigmayp2 - sigmayyp**2)*1e6
    return e
    


def get_normemittance(data):
    """ Calculates the normalized emittance of the beam from the data object.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    en : double
        The normalized emittance of the beam.
    """
    e = get_emittance(data)
    gamma = get_gamma(data)
    return e * gamma


def get_x(data):
    """ Get the array of x positions from the data object.
    """
    return data[:, 0]


def get_y(data):
    """ Get the array of y positions from the data object.
    """
    return data[:, 1]

def get_ux(data):
    """ Get the array of x velocities from the data object.
    """
    return data[:, 2]

def get_uy(data):
    """ Get the array of y velocities from the data object.
    """
    return data[:, 3]

def get_weights(data):
    """ Get the array of particle weights from the data object.
    """
    return data[:, 6]
