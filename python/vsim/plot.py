#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:55:08 2017

@author: robert
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from vsim.load import load_field, load_species
from scipy.ndimage.interpolation import zoom


def alpha_colormap_cutoff(cmap, cutoff, flip=True):
    """ Creates a colormap where the end beyond cutoff is transparent.

    Parameters
    ----------
    cmap : cmap
        The matplotlib colormap object to modify.
    cutoff : float
        Between 0-1, the percentage of the beginning of the colormap to change.
    flip : bool, optional
        Which side of the colormap to make transparent.

    Returns
    -------
    cmapt : cmap
        The partially transparent colormap.
    """
    N = cmap.N
    cmapt = cmap(np.arange(N))
    alpha = np.ones(N)
    if flip:
        alpha[:int(cutoff*N)] = 0.0
    else:
        alpha[int((1-cutoff)*N):] = 0.0
    cmapt[:, -1] = alpha
    cmapt = ListedColormap(cmapt)
    return cmapt


def drive_witness_density(params):
    """ Makes a plot of the drive and witness beam desnities with electrons.

    Shows the charge density of both the drive and witness beams in different
    colors. Shows the electrons as small dots to visualize the wake shape.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            drive : string
                The field name for the drive beam charge density.
            witness : string
                The field name for the witness beam charge denisty.
            plasma : string
                The species name of the plasma electrons.
            dumpInd : int
                The dump index to plot the beams at.
            path : string
                The path to the VSim output folder.
            simName : strign
                The simulation name, the first part of every simulation file.
            zoom : float
                Increase in the number of pixels, spline interpolation used.
    """
    # First load in all the necessary data
    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    zf = params['zoom']
    dData, dAttrs = load_field(path, simName, params['drive'])
    wData, wAttrs = load_field(path, simName, params['witness'])
    pData, pAttrs = load_species(path, simName, params['plasma'])
    # Grab the dump we are interested in
    if pAttrs['dim'] == 2:
        drive = zoom(np.flipud(np.transpose(dData[ind][:, :, 0])), zf)
        witness = zoom(np.flipud(np.transpose(wData[ind][:, :, 0])), zf)
        plasmaZ = pData[ind][:, 0] * 1e3
        plasmaX = pData[ind][:, 1] * 1e6

    plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[24, 1])
    # Create the two colormaps
    cmapD = alpha_colormap_cutoff(plt.cm.bone, 0.01, False)
    cmapW = alpha_colormap_cutoff(plt.cm.pink, 0.01, False)

    # Create the axis grid spaces
    colorAx1 = plt.subplot(gs[0, 1])
    colorAx2 = plt.subplot(gs[1, 1])
    # Convert units
    extent = np.array(dAttrs['bounds'][ind])
    extent[:2] *= 1e3
    extent[2:] *= 1e6

    # Create the actual plot
    plt.subplot(gs[:, 0])
    plt.imshow(drive, aspect='auto', extent=extent, cmap=cmapD)
    cb1 = plt.colorbar(cax=colorAx1)
    cb1.set_label(r'Drive beam - Charge density ($C/m^3$)')
    plt.imshow(witness, aspect='auto', extent=extent, cmap=cmapW)
    cb2 = plt.colorbar(cax=colorAx2)
    cb2.set_label(r'Wintess Beam - Charge density ($C/m^3$)')
    # Plasma electrons
    plt.plot(plasmaZ, plasmaX, 'bo', markersize=0.25)
    plt.xlabel(r'z ($mm$)')
    plt.ylabel(r'x ($\mu m$)')
    plt.title('Drive and witness beam charge density with electron wake')

    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'DriveWitnessDensity.pdf', format='pdf')
    plt.savefig(path+'DriveWitnessDensity.png', format='png')
    plt.show()
