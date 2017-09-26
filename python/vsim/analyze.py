#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:05:19 2017

@author: robert
"""

import scipy.constants as const
import numpy as np


def get_ptc_gamma(data):
    """ Calculates the relativistic factor for each particle in the beam.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    gamma : array-like
        The gamma of each particle in the beam.
    """
    ux = get_ux(data)
    uy = get_uy(data)
    u2 = ux**2 + uy**2
    gamma = 1/np.sqrt(1 - u2/(const.c**2 + u2))
    return gamma


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
    weights = get_weights(data)
    gamma = get_ptc_gamma(data)
    gammaAvg = np.average(gamma, weights=weights)
    return gammaAvg


def get_ptc_energy(data, mass):
    """ Calculates the energy of each particle from the data object for a beam.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.
    mass : double
        The mass of the particle in energy units. The units of energy determine
        the return units. me=0.511 MeV, returns energy in MeV.

    Returns
    -------
    energy : array-like
        The energy of each particle in the beam.
    """
    gamma = get_ptc_gamma(data)
    energy = gamma * mass
    return energy


def get_energy(data, mass):
    """ Calculates the energy of the beam from the data object for a beam.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.
    mass : double
        The mass of the particle in energy units. The units of energy determine
        the return units. me=0.511 MeV, returns energy in MeV.

    Returns
    -------
    energyAvg : double
        The average energy of all the particles in the beam.
    """
    gamma = get_gamma(data)
    energyAvg = gamma * mass
    return energyAvg


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


def get_density_temp(data, mass, attrs, mesh):
    """ Calculate the electron # density and temperature in each grid cell.
    
    Parameters
    ----------
    data : HDF5 dataset
        The data set for the particle species of interest, electrons generally.
        Use load.get_species_data to load the dataset object from a file.
    mass : double
        The mass of the particle in eV. 
    attrs : dictionary
        The attributes for the particle species, output of 
        load.get_species_attrs.
    mesh : dictionary
        The data about the mesh to deposit the temperature onto. 
        See load.get_mesh for the details of what should be in the dictionary.
    
    Returns
    -------
    den : array-like
        An array of particle number density on the grid.
    temp : array-like
        An array of the temperature on the grid.
    """
    Ndim = mesh['numCells'].shape[0]
    if Ndim == 2:
        Nx = mesh['numCells'][0]
        Ny = mesh['numCells'][1]
        XStart = attrs['lowerBounds'][0]
        YStart = attrs['lowerBounds'][1]
        dx = mesh['cellSize'][0]
        dy = mesh['cellSize'][1]
        x = get_x(data)
        y = get_y(data)
        weights = get_weights(data)
        energy = get_ptc_energy(data, mass) - mass
    temp = np.zeros(mesh['numCells'], dtype='double')
    const = 2 / (3*8.6173303e-5)
    den = np.zeros(mesh['numCells'], dtype='double')
    for i in range(Nx):
        xlo = XStart + i*dx
        xhi = xlo + dx
        xmask = np.ma.masked_outside(x, xlo, xhi)
        for j in range(Ny):
            ylo = YStart + i*dy
            yhi = ylo + dy
            ymask = np.ma.masked_outside(y, ylo, yhi)
            mask = np.logical_and(xmask.mask, ymask.mask)
            temp[i, j] = const*np.average(energy[mask], weights=weights[mask])
            den[i, j] = np.sum(weights[mask]) * attrs['ptsInMacro'] / dx / dy
    return den ,temp


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
    return data[:, -1]
