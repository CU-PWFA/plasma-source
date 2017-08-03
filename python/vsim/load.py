#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 11:01:34 2017

@author: robert
"""

import numpy as np
import h5py
import glob


def get_field_data(fileName, fieldName):
    """ Retrieves all the data for a given field from a field file.

    Loads the value of the field at each grid point from the passed file.

    Parameters
    ----------
    fileName : string
        The full path and filename to load.
    fieldName : string
        The name of the field to be loaded.

    Returns
    -------
    data : array-like
        The numpy array of the data at each point on the grid.
    """
    f = h5py.File(fileName, 'r') # r for read only, don't modify the data
    data = f[fieldName]
    return data


def get_field_attrs(fileName, fieldName):
    """ Retrieves useful attributes from a field file.
    
    Loads the useful attribute metadata from the passed file. The metadata
    includes things like the limits of the grid and the simulation time.

    Parameters
    ----------
    fileName : string
        The full path and filename to load.
    fieldName : string
        The name of the field to be loaded.

    Returns
    -------
    attrs : array-like
        Attribute dictionary containing useful information.
    """
    attrs = {}
    f = h5py.File(fileName, 'r') # r for read only, don't modify the data
    # Grad the grid limits
    limits = f[fieldName].attrs['vsLimits']
    limits = f[limits]
    lowBounds = limits.attrs['vsLowerBounds']
    upBounds = limits.attrs['vsUpperBounds']
    bounds = np.array([lowBounds[0], upBounds[0], lowBounds[1], upBounds[1]])
    attrs['bounds'] = bounds
    # Get the simulation time
    attrs['time'] = f[fieldName].attrs['time']
    return attrs


def load_field(path, simName, fieldName):
    """ Load all the files for a field and retrieve useful attributes.

    Loads all the dump files for a single field and combines them into a
    single numpy array. Creates an attribute dictionary of all the useful
    metadata.

    Parameters
    ----------
    path : string
        Path to the simulation folder.
    simName : string
        The simulation name, all the output files start with this string. It is
        also the name of the in file.
    fieldName : string
        The name of the field to be loaded.

    Returns
    -------
    data : dictionary
        Data array in the form (dump#, field).
    attrs : dictionary
        Attribute dictionary containing useful information.
    """
    attrs = {'time': {},
             'bounds': {}}
    data = {}
    files = glob.glob(path + simName + '_' + fieldName + '_*.h5')
    if len(files) > 0:
        for fileName in files:
            dumpInd = int(fileName.split('.h5')[-2].split('_')[-1])
            # Get the data from the file
            fData = get_field_data(fileName, fieldName)
            fAttrs = get_field_attrs(fileName, fieldName)
            attrs['time'][dumpInd] = fAttrs['time']
            attrs['bounds'][dumpInd] = fAttrs['bounds']
            data[dumpInd] = fData
    else:
        print('No files found matching the passed names.')
    return data, attrs
    

def get_species_data(fileName, speciesName):
    """ Retrieves all the data from a species file.

    Loads data for each particle in a species from the passed file.

    Parameters
    ----------
    fileName : string
        The full path and filename to load.
    speciesName : string
        The name of the species to be loaded.

    Returns
    -------
    data : array-like
        The numpy array of the data for each macroparticle.
    """
    f = h5py.File(fileName, 'r') # r for read only, don't modify the data
    data = f[speciesName]
    return data


def get_species_attrs(fileName, speciesName):
    """ Retrieves useful attributes from a species file.
    
    Loads the useful attribute metadata from the passed file. The metadata
    includes things like the labels of each data column and the time.

    Parameters
    ----------
    fileName : string
        The full path and filename to load.
    speciesName : string
        The name of the species to be loaded.

    Returns
    -------
    attrs : array-like
        Attribute dictionary containing useful information.
    """
    attrs = {}
    f = h5py.File(fileName, 'r') # r for read only, don't modify the data
    pts = f[speciesName].attrs
    # Get particle data
    attrs['charge'] = pts['charge']
    attrs['mass'] = pts['mass']
    attrs['ptsInMacro'] = pts['numPtclsInMacro']
    attrs['time'] = pts['time']
    attrs['dim'] = pts['numSpatialDims']
    # Get the data column headings
    attrs['columns'] = pts['vsLabels'].split(',')
    return attrs


def load_species(path, simName, speciesName):
    """ Load all the files for a species and retrieve useful attributes.

    Loads all the dump files for a single species and combines them into a
    single numpy array. Creates an attribute dictionary of all the useful
    metadata.

    Parameters
    ----------
    path : string
        Path to the simulation folder.
    simName : string
        The simulation name, all the output files start with this string. It is
        also the name of the in file.
    speciesName : string
        The name of the species to be loaded.

    Returns
    -------
    data : dictionary
        Data array in the form (dump#, species).
    attrs : dictionary
        Attribute dictionary containing useful information.
    """
    attrs = {'time': {}}
    data = {}
    files = glob.glob(path + simName + '_' + speciesName + '_*.h5')
    if len(files) > 0:
        for fileName in files:
            dumpInd = int(fileName.split('.h5')[-2].split('_')[-1])
            # Get the data from the file
            pData = get_field_data(fileName, speciesName)
            pAttrs = get_field_attrs(fileName, speciesName)
            attrs['time'][dumpInd] = pAttrs['time']
            data[dumpInd] = pData
        for name in pAttrs:
            if name is not 'time': attrs[name] = pAttrs[name]
    else:
        print('No files found matching the passed names.')
    return data, attrs
