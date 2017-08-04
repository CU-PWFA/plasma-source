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
from vsim import load
from scipy.ndimage.interpolation import zoom
# This needs to be before animation is imported
plt.rcParams['animation.ffmpeg_path'] = '/home/robert/anaconda3/envs/CU-PWFA/bin/ffmpeg'
import matplotlib.animation as animation


def get_filename(path, simName, name, dump):
    """ Returns the VSim filename. 
    
    Parameters
    ----------
    path : string
        The path to the VSim output folder.
    simName : string
        The simulation name, the first part of every simulation file.
    name : string
        The name of the field or species of interest.
    dump : int
        The dump number.
    """
    return path + simName + '_' + name + '_' + str(dump) + '.h5'

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
    # Grad items we use alot from the params
    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    zf = params['zoom']
    drive = params['drive']
    witness = params['witness']
    plasma = params['plasma']
    # First load in all the necessary data
    dFile = get_filename(path, simName, drive, ind)
    dData = load.get_field_data(dFile, drive)
    dAttrs = load.get_field_attrs(dFile, drive)

    wFile = get_filename(path, simName, witness, ind)
    wData = load.get_field_data(wFile, witness)

    pFile = get_filename(path, simName, plasma, ind)
    pData = load.get_species_data(pFile, plasma)
    pAttrs = load.get_species_attrs(pFile, plasma)
    
    # Up the resolution and convert units
    if pAttrs['dim'] == 2:
        driveData = zoom(np.flipud(np.transpose(dData[:, :, 0])), zf)
        witnessData = zoom(np.flipud(np.transpose(wData[:, :, 0])), zf)
        plasmaZ = pData[:, 0] * 1e3
        plasmaX = pData[:, 1] * 1e6
    else:
        print('Only 2D simulations are currently supported.')

    plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[24, 1])
    # Create the two colormaps
    cmapD = plt.cm.bone
    cmapW = alpha_colormap_cutoff(plt.cm.pink, 0.01, False)

    # Create the axis grid spaces
    colorAx1 = plt.subplot(gs[0, 1])
    colorAx2 = plt.subplot(gs[1, 1])
    # Convert units
    extent = np.array(dAttrs['bounds'])
    extent[:2] *= 1e3
    extent[2:] *= 1e6

    # Create the actual plot
    plt.subplot(gs[:, 0])
    plt.imshow(driveData, aspect='auto', extent=extent, cmap=cmapD)
    cb1 = plt.colorbar(cax=colorAx1)
    cb1.set_label(r'Drive beam - Charge density ($C/m^3$)')
    plt.imshow(witnessData, aspect='auto', extent=extent, cmap=cmapW)
    cb2 = plt.colorbar(cax=colorAx2)
    cb2.set_label(r'Wintess Beam - Charge density ($C/m^3$)')
    # Plasma electrons
    plt.plot(plasmaZ, plasmaX, 'bo', markersize=0.25)
    plt.xlabel(r'z ($mm$)')
    plt.ylabel(r'x ($\mu m$)')
    plt.title('Drive and witness beam charge density with electron wake')

    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'DriveWitnessDensity_'+str(ind)+'.pdf', format='pdf')
    plt.savefig(path+'DriveWitnessDensity_'+str(ind)+'.png', format='png')
    plt.show()


def phase_space(params):
    """ Plot the macroparticles of a beam in transverse phase space.

    Plots the macro particles in a beam in transverse phase space, each
    particle is color coded according to its weight, particles with weights
    below a certain cutoff are ignored.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            species : string
                The species name for the beam of interest.
            dumpInd : int
                The dump index to plot the beam at.
            path : string
                The path to the VSim output folder.
            simName : strign
                The simulation name, the first part of every simulation file.
            cutoff : float
                The minimum weight cutoff. Note maximum weight seems to be ~1.
    """
    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    species = params['species']
    # Load the particle data
    pFile = get_filename(path, simName, species, ind)
    pData = load.get_species_data(pFile, species)
    pAttrs = load.get_species_attrs(pFile, species)
    # Filter out super light particles, convert units
    if pAttrs['dim'] == 2:
        ptWeight = pData[:, 6]
        sel = ptWeight > params['cutoff']
        ptWeight = ptWeight[sel]
        ptX = pData[:, 1][sel] * 1e6
        ptUx = (pData[:, 3] / pData[:, 2])[sel]*1e3
    else:
        print('Only 2D simulations are currently supported.')

    # Sort the arrays so heavier particles appear on top
    sort = np.argsort(ptWeight)
    ptWeight = ptWeight[sort]
    ptX = ptX[sort]
    ptUx = ptUx[sort]
    # Create the plot
    plt.figure(figsize=(16, 9))
    plt.scatter(ptX, ptUx, c=ptWeight, cmap=plt.cm.hot_r)
    cb = plt.colorbar()
    cb.set_label('Particle weight')
    plt.xlabel(r'x ($\mu m$)')
    plt.ylabel(r"x' ($mrad$)")
    plt.title(species+' transverse phase space distribution')
    plt.grid(True)

    # Save the figure and display it
    plt.savefig(path+species+'PhaseSpace_'+str(ind)+'.pdf', format='pdf')
    plt.savefig(path+species+'PhaseSpace_'+str(ind)+'.png', format='png')
    plt.show()


def phase_space_animation(params):
    """ Animate the macroparticles of a beam in transverse phase space.

    Creates an animation showing the evolution of the beam in transverse phase
    space. Each macroparticle above the cutoff is plotted as the beam evolves.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            species : string
                The species name for the beam of interest.
            dumpInd : int
                The dump index to plot the beam at.
            path : string
                The path to the VSim output folder.
            simName : strign
                The simulation name, the first part of every simulation file.
            cutoff : float
                The minimum weight cutoff. Note maximum weight seems to be ~1.
            Nt : int
                The number of frames to render, set to the last dump number.
            fps : int
                The number of frames per second, determines animation length.
            xlim : array-like
                X plot limits, two element array [xmin, xmax].
            ylim : array-like
                Y plot limits, two element array [ymin, ymax].
    """
    Nt = params['Nt'] = 400
    path = params['path']
    simName = params['simName']
    species = params['species']
    
    pFile = get_filename(path, simName, species, 0)
    pAttrs = load.get_species_attrs(pFile, species)
    
    if pAttrs['dim'] != 2: 
        print('Only 2D simulations are currently supported.')
        return
    
    # Grab the dump we are interested in
    def get_data(ind):
        pFile = get_filename(path, simName, species, ind)
        pData = load.get_species_data(pFile, species)
        ptWeight = pData[:, 6]
        sel = ptWeight > params['cutoff']
        ptWeight = ptWeight[sel]
        ptX = pData[:, 1][sel] * 1e6
        ptUx = (pData[:, 3] / pData[:, 2])[sel]*1e3
        # Sort the arrays so heavier particles appear on top
        sort = np.argsort(ptWeight)
        ptWeight = ptWeight[sort]
        ptX = ptX[sort]
        ptUx = ptUx[sort]
        return ptWeight, ptX, ptUx
    
    # Get the first piece of data
    ptWeight, ptX, ptUx = get_data(0)

    # Create the plot
    fig = plt.figure(figsize=(16, 9))
    sct = plt.scatter(ptX, ptUx, c=ptWeight, cmap=plt.cm.hot_r)
    cb = plt.colorbar()
    cb.set_label('Particle weight')
    plt.xlabel(r'x ($\mu m$)')
    plt.ylabel(r"x' ($mrad$)")
    plt.title(species+' transverse phase space distribution')
    plt.xlim(params['xlim'])
    plt.ylim(params['ylim'])
    plt.grid(True)

    # Update the scatter plot data
    i = 1;
    def updatefig(*args):
        nonlocal i
        ptWeight, ptX, ptUx = get_data(i)
        sct.set_offsets(np.stack((ptX, ptUx), axis=-1))
        sct.set_array(ptWeight)
        i += 1
        # If we run over, loop
        if i == Nt+1:
            i = 0
        return sct,

    ani = animation.FuncAnimation(fig, updatefig, blit=True, frames=Nt)
    ani.save(params['path']+'PhaseSpaceEvolution.mp4', fps=params['fps'])
