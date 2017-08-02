#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 14:28:31 2017

@author: robert
"""

import sys
# TODO need to make this something that can be set
VSimPath = "/home/robert/VSim-8.1/"
# Path to VSim scripts
sys.path.append(VSimPath + "Contents/engine/share/scripts")
# VSim packages
sys.path.append(VSimPath + "Contents/engine/lib/python2.7/site-packages")

import glob
import VsHdf5
import numpy as np


def load_particle_file(fileName, speciesName):
    """ Loads the data for a particle species from a signle dump file.

    Parameters
    ----------
    fileName : string
        The filename to load the particle data from.
    speciesName : string
        The name of the partical species in VSim.

    Returns
    -------
    pData : array-like
        Array of particle data for all the macroparticles.
    pAttrs : array-like
        Array of tuples specifying properties about the data and the particles.
    """
    particles = VsHdf5.Particles(speciesName)
    pData, pAttrs = particles.readParticles(fileName, speciesName)
    return pData, pAttrs


def get_particle_attrs(data, pAttrs):
    """ Parses the particle attributes and adds useful items to data.

    This function parses the the pAttrs array and formats and adds useful
    information to the data object

    Parameters
    ----------
    data : dictionary
        Dictionary object to add attribute data to.
    pAttrs : array-like
        Attribute array that results from VSim internal scripts.

    Returns
    -------
    data : dictionary
        The updated dictionary with attribute data.
    """
    ind = 0
    data['dataInd'] = {}
    for attr in pAttrs:
        name = attr[0]
        value = attr[1]
        if name == 'charge' or name == 'mass' or name == 'numPtclsInMacro' or name == 'numSpatialDims':
            data[name] = value
        elif name == 'time':
            # Save the index of the time tuple so we can get times later
            data['timeInd'] = ind
        elif name == 'vsLabels':
            # Extract the name and index of each saved variable
            labels = value.split(',')
            for i in range(0, len(labels)):
                data['dataInd'][labels[i]] = i
        ind += 1
    return data


def load_species(path, simName, speciesName):
    """ Loads all the dump data for a given species.

    Loads all the useful data and attributes for one species and creates a
    single object with all the data in it.

    Parameters
    ----------
    path : string
        Path to the VSim folder for the simulation of interest.
    simName : string
        The name of the VSim simulation (name of the input file).
    speciesName : string
        The name of the partical species in VSim.

    Returns
    -------
    data : dictionary
        Dictionary containing all the output data an attributes for 1 species.
    """
    # Grab all the input files
    files = glob.glob(path + simName + '_' + speciesName + '_*.h5')
    data = {}
    data['time'] = {}
    data['dumpInd'] = []
    first = True # We have to grab attriubtes on the first run
    if len(files) > 0:
        for fileName in files:
            print('Reading ' + fileName)
            dumpInd = int(fileName.split('.h5')[-2].split('_')[-1])
            pData, pAttrs = load_particle_file(fileName, speciesName)
            data[dumpInd] = pData
            if first:
                get_particle_attrs(data, pAttrs) # Add attributes to data
                first = False
            data['time'][dumpInd] = pAttrs[data['timeInd']][1]
            data['dumpInd'].append(dumpInd)
        # Add this as an array to loop through
        data['dumpInd'] = np.sort(data['dumpInd'])
    else:
        print('No files found matching the passed names.')
    return data


def save_species(path, simName, speciesName):
    """ Saves all the data about a species from VSim in numpy formay.

    This function saves all the data produced by VSim about a single species in
    the numpy format which can be read by Python 3 code for analysis. The file
    is saved in the VSim directory with the name "path + simName + '_' +
    speciesName + 'Python.npy'".

    Parameters
    ----------
    path : string
        Path to the VSim folder for the simulation of interest.
    simName : string
        The name of the VSim simulation (name of the input file).
    speciesName : string
        The name of the partical species in VSim.
    """
    data = load_species(path, simName, speciesName)
    np.save(path + simName + '_' + speciesName + 'Python', data)
