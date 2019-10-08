#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:05:19 2017

@author: robert
"""

import scipy.constants as const
import numpy as np
import sys

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
    dim = int(data.attrs['numSpatialDims'])
    ux = get_ux(data,dim)
    uy = get_uy(data,dim)
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

def get_gamma_rms(data):
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
    
    average = np.average(gamma, weights=weights)
    variance = np.average((gamma-average)**2, weights=weights)
    
    return np.sqrt(variance)/average*100

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
    """ Calculates the emittance of the beam from the data object (in Vsim y (x)).

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
    dim = int(data.attrs['numSpatialDims'])
    
    y = get_y(data,dim)
    ux = get_ux(data,dim)
    uy = get_uy(data,dim)
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


def get_emittance_y(data):
    """ Calculates the emittance of the beam from the data object (in VSim z).

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
    dim = int(data.attrs['numSpatialDims'])
    
    z = get_z(data,dim)
    ux = get_ux(data,dim)
    uz = get_uz(data,dim)
    weights = get_weights(data)
    zp = uz / ux
    dz = z - np.average(z, weights=weights)
    dzp = zp - np.average(zp, weights=weights)
    # Calculate the RMS sizes and the correlation
    sigmaz2 = np.average(dz**2, weights=weights)
    sigmazp2 = np.average(dzp**2, weights=weights)
    sigmazzp = np.average(dz*dzp, weights=weights)
    # Calculate the emittance
    e = np.sqrt(sigmaz2*sigmazp2 - sigmazzp**2)*1e6
    return e


def get_sigmar(data):
    """ Calculates the sigmar of the beam from the data object.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    sigmar : double
        The sigmar of the beam in m.
    """
    y = get_y(data)
    weights = get_weights(data)
    dy = y - np.average(y, weights=weights)
    sigmar = np.sqrt(np.average(dy**2, weights=weights))
    return sigmar


def get_sigmar_y(data):
    """ Calculates the sigmar of the beam from the data object in y.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    sigmar : double
        The sigmar of the beam in m.
    """
    dim = int(data.attrs['numSpatialDims'])
    
    z = get_z(data, dim)
    weights = get_weights(data)
    dz = z - np.average(z, weights=weights)
    sigmar = np.sqrt(np.average(dz**2, weights=weights))
    return sigmar


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


def get_sigma_xp(data):
    """ Calculates the sigmar of the beam from the data object in xp.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    sigma_xp : double
        The sigmar of the beam in rad.
    """
    dim = int(data.attrs['numSpatialDims'])
    
    ux = get_ux(data,dim)
    uy = get_uy(data,dim)
    weights = get_weights(data)
    yp = uy / ux
    # Calculate the RMS sizes and the correlation
    dyp = yp - np.average(yp, weights=weights)
    sigma_yp = np.sqrt(np.average(dyp**2, weights=weights))
    return sigma_yp


def get_sigma_yp(data):
    """ Calculates the sigmar of the beam from the data object in yp.

    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_species_data to
        load the dataset object from a file.

    Returns
    -------
    sigma_yp : double
        The sigmar of the beam in rad.
    """
    dim = int(data.attrs['numSpatialDims'])
    
    ux = get_ux(data,dim)
    uz = get_uz(data,dim)
    weights = get_weights(data)
    zp = uz / ux
    # Calculate the RMS sizes and the correlation
    dzp = zp - np.average(zp, weights=weights)
    sigma_zp = np.sqrt(np.average(dzp**2, weights=weights))
    return sigma_zp


def find_length(rhoX, x):
    """ Finds the start and end of the wake, the difference is the length.
    
    Parameters
    ----------
    rhoX : array-like
        A 1D array of the plasma density along axis.
    x : array-like
        X positions for each of the elements in rhoX, result will share units.
        
    Returns
    -------
    start : double
        The start of the wake.
    end : double 
        The end of the wake.
    """
    Nx = len(x)
    eps = 1e6
    for i in range(Nx-1, -1, -1):
        if rhoX[i] >= eps and rhoX[i-1] < eps:
            start = x[i-1]
        if rhoX[i] <= eps and rhoX[i-1] > eps:
            end = x[i]
            return start, end
   
     
def find_width(rhoY, y):
    """ Finds the bottom and top of the wake at a given x position.
    
    Parameters
    ----------
    rhoY : array-like
        A 1D array of the plasma density transverse to the propagation.
    y : array-like
        y positions for each of the elements in rhoX, result will share units.
        
    Returns
    -------
    start : double
        The start of the wake.
    end : double 
        The end of the wake.
    """
    Ny = len(y)
    eps = 1e6
    start = None
    for i in range(int(Ny/4), int(3*Ny/4)):
        if rhoY[i] >= eps and rhoY[i+1] < eps:
            start = y[i+1]
        if rhoY[i] <= eps and rhoY[i+1] > eps:
            end = y[i]
            if start is None:
                return None, None
            else:
                return start, end
    return None, None


def find_wake_width(rhoXY, y):
    """ Finds the maximum width of the wake for all x positions.
    
    Parameters
    ----------
    rhoXY : array-like
        A 2D array of the plasma density in a transverse slice.
    y : array-like
        y positions for each of the elements in rhoX, result will share units.
        
    Returns
    -------
    start : double
        The start of the wake.
    end : double 
        The end of the wake.
    """
    Nx, Ny = np.shape(rhoXY)
    beg = 0
    end = Nx
    widthArr = np.zeros(end-beg, dtype='double')
    for j in range(beg, end):
        rhoYj = rhoXY[j, :]
        width = find_width(rhoYj, y)
        if width[0] is not None:
            widthArr[j-beg] = width[1]-width[0]
        else:
            widthArr[j-beg] = 0
    rhoY = rhoXY[np.argmax(widthArr)+beg, :]
    return find_width(rhoY, y)


def get_shape(data):
    """ Returns the grid size of a field data object
    
    Parameters
    ----------
    data : HDF5 dataset
        The data set for the beam of interest, use load.get_field_data or 
        load_field to load the dataset object from a file.
    
    Returns
    -------
    Nx, Ny, Nz : int
        The size of the grid, returns only Nx and Ny if 2D.
    """
    shape = data.shape
    if len(shape) == 4: #3D
        return shape[0], shape[1], shape[1]
    if len(shape) == 3: #2D
        return shape[0], shape[1]


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


def get_x(data, dim=2):
    """ Get the array of x positions from the data object.
    """
    if dim == 2:
        ind = 0
    elif dim == 3:
        ind = 0
    return data[:, ind]


def get_y(data, dim=2):
    """ Get the array of y positions from the data object.
    """
    if dim == 2:
        ind = 1
    elif dim == 3:
        ind = 1
    return data[:, ind]

def get_z(data, dim=3):
    """ Get the array of y positions from the data object.
    """
    if dim == 2:
        print("z Coord. not in 2D!")
    elif dim == 3:
        ind = 2
    return data[:, ind]

def get_ux(data, dim=2):
    """ Get the array of x velocities from the data object.
    """
    if dim == 2:
        ind = 2
    elif dim == 3:
        ind = 3
    return data[:, ind]


def get_uy(data, dim=2):
    """ Get the array of y velocities from the data object.
    """
    if dim == 2:
        ind = 3
    elif dim == 3:
        ind = 4
    return data[:, ind]

def get_uz(data, dim=3):
    """ Get the array of y velocities from the data object.
    """
    if dim == 2:
        ind = 4
    elif dim == 3:
        ind = 5
    return data[:, ind]


def get_weights(data):
    """ Get the array of particle weights from the data object.
    """
    return data[:, -1]
