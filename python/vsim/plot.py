#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:55:08 2017

@author: robert
"""

from vsim import C
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
from vsim import load
from vsim import analyze
from scipy.ndimage.interpolation import zoom
# This needs to be before animation is imported
plt.rcParams['animation.ffmpeg_path'] = '/home/robert/anaconda3/envs/CU-PWFA/bin/ffmpeg'
import matplotlib.animation as animation
from scipy import constants
import scipy.constants as const
import matplotlib.colors as colors


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
            alphaCutoff : float
                Must be between 0-1, start with 0.05, if you can't see the
                drive beam, increase
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
        dData2D = dData[:,:,int(dData.shape[2]/2)].reshape((dData.shape[0],dData.shape[1]))
        driveData = zoom(np.flipud(np.transpose(dData2D)), zf)
        wData2D = wData[:,:,int(wData.shape[2]/2)].reshape((wData.shape[0],wData.shape[1]))
        witnessData = zoom(np.flipud(np.transpose(wData2D)), zf)
        plasmaZ = pData[:, 0] * 1e3
        plasmaX = pData[:, 1] * 1e6
        
        #Now only plot plasma e's within 2um of the center
        plasmaY = pData[:, 2]
        ycenter = np.array(np.where(np.abs(plasmaY) < 2e-6)[0])
        plasmaZ = plasmaZ[ycenter]
        plasmaX = plasmaX[ycenter]

    plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[24, 1])
    
    # Create the two colormaps
    cmapD = plt.cm.get_cmap('bone')
    cmapW = alpha_colormap_cutoff(plt.cm.get_cmap('pink'),
                                  params['alphaCutoff'], False)
        
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

def single_drive_density(params):
    """ Makes a plot of the drive beam desnity with electrons.

    Shows the charge density of both the drive and witness beams in different
    colors. Shows the electrons as small dots to visualize the wake shape.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            drive : string
                The field name for the drive beam charge density.
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
            alphaCutoff : float
                Must be between 0-1, start with 0.05, if you can't see the
                drive beam, increase
    """
    # Grad items we use alot from the params
    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    zf = params['zoom']
    drive = params['drive']
    plasma = params['plasma']
    # First load in all the necessary data
    dFile = get_filename(path, simName, drive, ind)
    dData = load.get_field_data(dFile, drive)
    dAttrs = load.get_field_attrs(dFile, drive)

    pFile = get_filename(path, simName, plasma, ind)
    pData = load.get_species_data(pFile, plasma)
    pAttrs = load.get_species_attrs(pFile, plasma)
    
    # Up the resolution and convert units
    if pAttrs['dim'] == 2:
        driveData = zoom(np.flipud(np.transpose(dData[:, :, 0])), zf)
        plasmaZ = pData[:, 0] * 1e3
        plasmaX = pData[:, 1] * 1e6
    else:
        dData2D = dData[:,:,int(dData.shape[2]/2)].reshape((dData.shape[0],dData.shape[1]))
        driveData = zoom(np.flipud(np.transpose(dData2D)), zf)
        plasmaZ = pData[:, 0] * 1e3
        plasmaX = pData[:, 1] * 1e6
        
        #Now only plot plasma e's within 2um of the center
        plasmaY = pData[:, 2]
        ycenter = np.array(np.where(np.abs(plasmaY) < 2e-6)[0])
        plasmaZ = plasmaZ[ycenter]
        plasmaX = plasmaX[ycenter]

    plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[24, 1])
    # Create the two colormaps
    cmapD = plt.cm.get_cmap('bone')

    # Create the axis grid spaces
    colorAx1 = plt.subplot(gs[0, 1])
    # Convert units
    extent = np.array(dAttrs['bounds'])
    extent[:2] *= 1e3
    extent[2:] *= 1e6

    # Create the actual plot
    plt.subplot(gs[:, 0])
    plt.imshow(driveData, aspect='auto', extent=extent, cmap=cmapD)
    cb1 = plt.colorbar(cax=colorAx1)
    cb1.set_label(r'Drive beam - Charge density ($C/m^3$)')
    # Plasma electrons
    plt.plot(plasmaZ, plasmaX, 'bo', markersize=0.25)
    plt.xlabel(r'z ($mm$)')
    plt.ylabel(r'x ($\mu m$)')
    plt.title('Drive and witness beam charge density with electron wake')

    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'SingleDriveDensity_'+str(ind)+'.pdf', format='pdf')
    plt.savefig(path+'SingleDriveDensity_'+str(ind)+'.png', format='png')
    plt.show()

def drive_witness_density_new(params):
    e = const.physical_constants['elementary charge'][0]

    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    zf = params['zoom']
    drive = params['drive']
    plasma = params['plasma']
    
    # Load in plasma density
    rho, rhoAttrs = load.load_field(path, simName, 'rhoPlasma')
    Nx, Ny, Nz = analyze.get_shape(rho[ind])
    #rhoXY = -np.transpose(rho[ind][:, :, int(Nz+1)/2, 0]/e/1e6)+2 #for one transverse dim
    rhoXY = -np.transpose(rho[ind][:, int(Nz+1)/2, :, 0]/e/1e6)+2 #for the other tran dim
    
    #Load in drive beam density
    rho, rhoAttrs = load.load_field(path, simName, 'rhoDrive')
    Nx, Ny, Nz = analyze.get_shape(rho[ind])
    rhoBXY = -np.transpose(rho[ind][:, :, int(Nz+1)/2, 0]/e/1e6)
    
    #Load in witness beam density
    rho, rhoAttrs = load.load_field(path, simName, 'rhoWitness')
    Nx, Ny, Nz = analyze.get_shape(rho[ind])
    rhoWXY = -np.transpose(rho[ind][:, :, int(Nz+1)/2, 0]/e/1e6)
    
    def alpha_colormap(cmap, cutoff, flip=True):
        N = cmap.N
        cmapt = cmap(np.arange(N))
        alpha = np.ones(N)
        if flip:
            temp = alpha[:int(cutoff*N)]
            M = len(temp)
            alpha[:int(cutoff*N)] = np.linspace(0, 1, M)
        else:
            alpha[int((1-cutoff)*N):] = 0.0
        cmapt[:, -1] = alpha
        cmapt = colors.ListedColormap(cmapt)
        return cmapt
    
    # Plot the plasma density
    fig = plt.figure(figsize=(10, 7.5), dpi=100, frameon=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    
    # Plot the drive beam
    ax.imshow(rhoBXY, interpolation='gaussian', extent=[-125, 125, -125, 125], cmap='copper')
    
    # Plot the witness beam
    #cmapW = alpha_colormap(plt.cm.get_cmap('nipy_spectral'), 0.1, True)
    cmapW = alpha_colormap(plt.cm.get_cmap('rainbow'), 0.1, True)
    ax.imshow(rhoWXY, interpolation='gaussian', extent=[-125, 125, -125, 125], cmap=cmapW)
    npcase = 2e16
    vmax = 1e18*(npcase/3.7e17)
    vmin = 3e16*(npcase/3.7e17)
    
    # Plot the plasma density
    cmapP = alpha_colormap(plt.cm.get_cmap('inferno'), 0.2, True)
    ax.imshow(rhoXY, interpolation='gaussian', aspect='auto', extent=[-125, 125, -125, 125],
               norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cmapP)#1e16,2e18
    
    plt.savefig(path+'Title_Wake.png')
    plt.show()

def wake_cross_section(params):
    e = const.physical_constants['elementary charge'][0]

    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    vmin = params['vmin']
    vmax = params['vmax']
    central_off = params['centralOff']
    tranExtent = params['tranExtent']
    plot = params['plot']
    
    def alpha_colormap(cmap, cutoff, flip=True):
        N = cmap.N
        cmapt = cmap(np.arange(N))
        alpha = np.ones(N)
        if flip:
            temp = alpha[:int(cutoff*N)]
            M = len(temp)
            alpha[:int(cutoff*N)] = np.linspace(0, 1, M)
        else:
            alpha[int((1-cutoff)*N):] = 0.0
        cmapt[:, -1] = alpha
        cmapt = colors.ListedColormap(cmapt)
        return cmapt
    
    # Load in plasma density
    rho, rhoAttrs = load.load_field(path, simName, 'rhoPlasma')
    #print(rho); print(np.shape(rho))
    #print(rhoAttrs); print(np.shape(rhoAttrs))
    Nx, Ny, Nz = analyze.get_shape(rho[ind])
    #rhoYZ = -np.transpose(rho[ind][int(Nx+1)/2+central_off, :, :, 0]/e/1e6)+2
    rhoYZ = -(rho[ind][int(Nx+1)/2+central_off, :, :, 0]/e/1e6)+2
    x = np.linspace(-tranExtent, tranExtent, Nz)
    y = np.linspace(-tranExtent, tranExtent, Ny)
    """
    rho_grad = rhoYZ[:,int((Ny+1)/2)]
    plt.plot(y,rho_grad)
    plt.title("Linear Density Gradient")
    plt.xlabel(r'$x (\mu m)$')
    plt.ylabel(r'$n_p (cm^{-3})$')
    plt.grid();plt.show()
    """
    if plot is True:
        # Plot the plasma density
        fig = plt.figure(figsize=(7.5, 7.5), dpi=100, frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        #ax.set_axis_off()
        fig.add_axes(ax)
        cmapP = alpha_colormap(plt.cm.get_cmap('inferno'), 0.2, True)
        ax.imshow(rhoYZ, interpolation='gaussian', aspect=1, extent=[-tranExtent, tranExtent, -tranExtent, tranExtent],
                   norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cmapP)#1e16,2e18
        """
        r_circ = 80
        plt.plot(x, np.sqrt(r_circ**2-np.square(x)),  c='b',ls='--')
        plt.plot(x, -np.sqrt(r_circ**2-np.square(x)), c='b',ls='--')
        """
        plt.savefig(path+'Title_Wake.png')
        plt.show()
    
    r = np.zeros(Nz*Ny); theta = np.zeros(Nz*Ny); rhoPol = np.zeros(Nz*Ny)
    z=0
    for i in range(Nz):
        for j in range(Ny):
            r[z] = np.sqrt(np.square(x[i])+np.square(y[j]))
            if x[i] == 0:
                if y[j] >= 0:
                    theta[z] = np.pi/2
                else:
                    theta[z] = -np.pi/2
            else:
                theta[z] = np.arctan(y[j]/x[i])
            if x[i] < 0:
                theta[z] = theta[z] + np.pi
            rhoPol[z] = rhoYZ[i,j]
            z = z+1
            
    threshset = np.array(np.where(r < 100)[0])
    theta = theta[threshset]; r = r[threshset]; rhoPol = rhoPol[threshset]
    sort = np.argsort(theta)
    theta = theta[sort]; r = r[sort]; rhoPol = rhoPol[sort]
    
    thetabins = np.linspace(min(theta), max(theta), 60)
    rhosheath = np.zeros(len(thetabins)-1)
    rsheath = np.zeros(len(thetabins)-1)

    k = 0
    for i in range(len(thetabins)-1):
        while theta[k] < thetabins[i+1]:
            
            if rhoPol[k] < 1e16:
                if r[k] > rsheath[i]:
                    rsheath[i] = r[k]
            """
            if rhosheath[i] < rhoPol[k]:
                rhosheath[i] = rhoPol[k]
                rsheath[i] = r[k]
            """
            k=k+1
            
    return (thetabins[1:]+.5*(thetabins[2]-thetabins[1])), rsheath
    #return r, theta, rhoPol

def wake_cross_section_field(params):
    e = const.physical_constants['elementary charge'][0]

    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    central_off = params['centralOff']
    tranExtent = params['tranExtent']
    plot = params['plot']
    vec = params['vector']
    field = params['field']
    
    # Load in plasma density
    efield, efieldAttrs = load.load_field(path, simName, field)
    Nx, Ny, Nz = analyze.get_shape(efield[ind])
    #rhoYZ = -np.transpose(rho[ind][int(Nx+1)/2+central_off, :, :, 0]/e/1e6)+2
    eYZ = -(efield[ind][int(Nx+1)/2+central_off, :, :, vec])

    if plot is True:
        # Plot the plasma density
        fig = plt.figure(figsize=(7.5, 7.5), dpi=100, frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        #ax.set_axis_off()
        fig.add_axes(ax)
        h = ax.imshow(eYZ, interpolation='gaussian', aspect=1, extent=[-tranExtent, tranExtent, -tranExtent, tranExtent])
        plt.savefig(path+'Title_Wake.png')
        plt.colorbar(h)
        plt.show()
    
    evx = np.flip(eYZ[int((Nz+1)/2),:],0)
    evy = np.flip(eYZ[:,int((Nz+1)/2)],0)
    x = np.linspace(-tranExtent,tranExtent,Nz)
    y = np.linspace(-tranExtent,tranExtent,Ny)

    return x, y, evx, evy, eYZ

def drive_witness_animation(params):
    """ Creates an animation of the electron wake and beam density evolution.

    An animation of drive_witness_density where each dump is a frame in the
    animation.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
    """
    path = params['path']
    simName = params['simName']
    ind = params['dumpInd']
    zf = params['zoom']
    drive = params['drive']
    witness = params['witness']
    plasma = params['plasma']
    Nt = params['Nt']
    
    dFile = get_filename(path, simName, drive, ind)
    dAttrs = load.get_field_attrs(dFile, drive)
    
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[24, 1])
    # Create the two colormaps
    cmapD = plt.cm.get_cmap('bone')
    cmapW = alpha_colormap_cutoff(plt.cm.get_cmap('pink'),
                                  params['alphaCutoff'], False)

    # Create the axis grid spaces
    colorAx1 = plt.subplot(gs[0, 1])
    colorAx2 = plt.subplot(gs[1, 1])
    # Convert units
    extent = np.array(dAttrs['bounds'])
    extent[:2] *= 1e3
    extent[2:] *= 1e6
    
    # Grab the dump we are interested in
    def get_data(ind):
        dFile = get_filename(path, simName, drive, ind)
        wFile = get_filename(path, simName, witness, ind)
        pFile = get_filename(path, simName, plasma, ind)
        dData = load.get_field_data(dFile, drive)
        wData = load.get_field_data(wFile, witness)
        pData = load.get_species_data(pFile, plasma)
        
        driveData = zoom(np.flipud(np.transpose(dData[:, :, 0])), zf)
        witnessData = zoom(np.flipud(np.transpose(wData[:, :, 0])), zf)
        plasmaZ = pData[:, 0] * 1e3
        plasmaX = pData[:, 1] * 1e6
        return driveData, witnessData, plasmaZ, plasmaX
    
    driveData, witnessData, plasmaZ, plasmaX = get_data(ind)
    
    # Create the actual plot
    plt.subplot(gs[:, 0])
    im1 = plt.imshow(driveData, aspect='auto', animated=True,
                     extent=extent, cmap=cmapD)
    cb1 = plt.colorbar(cax=colorAx1)
    cb1.set_label(r'Drive beam - Charge density ($C/m^3$)')
    im2 = plt.imshow(witnessData, aspect='auto', animated=True,
                     extent=extent, cmap=cmapW)
    cb2 = plt.colorbar(cax=colorAx2)
    cb2.set_label(r'Wintess Beam - Charge density ($C/m^3$)')
    # Plasma electrons
    ax1 = plt.plot(plasmaZ, plasmaX, 'bo', markersize=0.25, animated=True)
    plt.xlabel(r'z ($mm$)')
    plt.ylabel(r'x ($\mu m$)')
    plt.title('Drive and witness beam charge density with electron wake')
    plt.tight_layout()
    
    # Update the scatter plot data
    i = ind+1;
    def updatefig(*args):
        nonlocal i
        driveData, witnessData, plasmaZ0, plasmaX = get_data(i)
        im1.set_array(driveData)
        im2.set_array(witnessData)
        print(ax1) #XXX I don't know, this is the correct object
        ax1.get_data(plasmaZ, plasmaX)
        i += 1
        # If we run over, loop
        if i == Nt+1:
            i = ind
        return im1, im2, ax1

    ani = animation.FuncAnimation(fig, updatefig, blit=True, frames=Nt)
    ani.save(path+'WakefieldEvolution.mp4', fps=params['fps'])


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
    plt.scatter(ptX, ptUx, c=ptWeight, cmap=plt.cm.get_cmap('hot_r'))
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
            size : double, optional
                The size of the particles in the plot.
    """
    Nt = params['Nt']
    path = params['path']
    simName = params['simName']
    species = params['species']
    ind = params['dumpInd']
    
    pFile = get_filename(path, simName, species, 0)
    pAttrs = load.get_species_attrs(pFile, species)
    
    if 'size' in params:
        size = params['size']
    else:
        size = 1.0
    
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
    ptWeight, ptX, ptUx = get_data(ind)

    # Create the plot
    fig = plt.figure(figsize=(16, 9))
    sct = plt.scatter(ptX, ptUx, c=ptWeight, cmap=plt.cm.get_cmap('hot_r'), 
                      s=size)
    cb = plt.colorbar()
    cb.set_label('Particle weight')
    plt.xlabel(r'x ($\mu m$)')
    plt.ylabel(r"x' ($mrad$)")
    plt.title(species+' transverse phase space distribution')
    plt.xlim(params['xlim'])
    plt.ylim(params['ylim'])
    plt.grid(True)

    # Update the scatter plot data
    i = ind+1;
    def updatefig(*args):
        nonlocal i
        ptWeight, ptX, ptUx = get_data(i)
        sct.set_offsets(np.stack((ptX, ptUx), axis=-1))
        sct.set_array(ptWeight)
        i += 1
        # If we run over, loop
        if i == Nt+1:
            i = ind
        return sct,

    ani = animation.FuncAnimation(fig, updatefig, blit=True, frames=Nt-ind-2)
    ani.save(params['path']+'PhaseSpaceEvolution_'+species+'.mp4',
             fps=params['fps'])


def phase_space_energy(params):
    """ Plot the macroparticles of a beam in transverse phase space.

    Plots the macro particles in a beam in transverse phase space, each
    particle is color coded according to its energy with an alpha corresponding
    to its weight, particles with weights below a certain cutoff are ignored.

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
    mass = 1e-6*pAttrs['mass']*constants.c**2 / constants.e # convert to MeV
    # Filter out super light particles, convert units
    if pAttrs['dim'] == 2:
        ptWeight = pData[:, 6]
        sel = ptWeight > params['cutoff']
        ptWeight = ptWeight[sel]
        ptX = pData[:, 1][sel] * 1e6
        ptUx = (pData[:, 3] / pData[:, 2])[sel]*1e3
        ptEnergy = analyze.get_ptc_energy(pData, mass)[sel]
    else:
        print('Only 2D simulations are currently supported.')

    if 'cmap' in params: cmap = params['cmap']
    else: cmap = 'gnuplot'
    
    # Sort the arrays so higher energy particles appear on top
    sort = np.argsort(ptEnergy)
    ptWeight = ptWeight[sort]
    ptX = ptX[sort]
    ptUx = ptUx[sort]
    ptEnergy = ptEnergy[sort]
    # Create the colormap
    norm = Normalize(vmin=np.amin(ptEnergy), vmax=np.amax(ptEnergy))
    alpha = ptWeight / np.amax(ptWeight)
    colors = plt.cm.get_cmap(cmap)(norm(ptEnergy))
    colors[:, 3] = alpha
    # Create the plot
    plt.figure(figsize=(16, 9))
    plt.scatter(ptX, ptUx, c=colors)
    #cb = plt.colorbar()
    #cb.set_label('Particle weight')
    plt.xlabel(r'x ($\mu m$)')
    plt.ylabel(r"x' ($mrad$)")
    plt.title(species+' transverse phase space distribution')
    plt.grid(True)

    # Save the figure and display it
    plt.savefig(path+species+'PhaseSpaceEnergy_'+str(ind)+'.pdf', format='pdf')
    plt.savefig(path+species+'PhaseSpaceEnergy_'+str(ind)+'.png', format='png')
    plt.show()


def phase_space_energy_animation(params):
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
            cmap : string, optional
                Colormap name to use for the energy.
    """
    Nt = params['Nt']
    path = params['path']
    simName = params['simName']
    species = params['species']
    ind = params['dumpInd']
    
    pFile = get_filename(path, simName, species, 0)
    pAttrs = load.get_species_attrs(pFile, species)
    
    if pAttrs['dim'] != 2: 
        print('Only 2D simulations are currently supported.')
        return
    if 'cmap' in params: cmap = params['cmap']
    else: cmap = 'gnuplot'
    
    mass = 1e-6*pAttrs['mass']*constants.c**2 / constants.e # convert to MeV
    # Grab the dump we are interested in
    def get_data(ind):
        pFile = get_filename(path, simName, species, ind)
        pData = load.get_species_data(pFile, species)
        ptWeight = pData[:, 6]
        sel = ptWeight > params['cutoff']
        ptWeight = ptWeight[sel]
        ptX = pData[:, 1][sel] * 1e6
        ptUx = (pData[:, 3] / pData[:, 2])[sel]*1e3
        ptEnergy = analyze.get_ptc_energy(pData, mass)[sel]
        # Sort the arrays so heavier particles appear on top
        sort = np.argsort(ptEnergy)
        ptWeight = ptWeight[sort]
        ptX = ptX[sort]
        ptUx = ptUx[sort]
        ptEnergy = ptEnergy[sort]
        # Create the colormap
        norm = Normalize(vmin=np.amin(ptEnergy), vmax=np.amax(ptEnergy))
        alpha = ptWeight / np.amax(ptWeight)
        colors = plt.cm.get_cmap(cmap)(norm(ptEnergy))
        colors[:, 3] = alpha
        return ptWeight, ptX, ptUx, colors
    
    # Get the first piece of data
    ptWeight, ptX, ptUx, colors = get_data(ind)

    # Create the plot
    fig = plt.figure(figsize=(16, 9))
    sct = plt.scatter(ptX, ptUx, c=colors)
    #cb = plt.colorbar()
    #cb.set_label('Particle weight')
    plt.xlabel(r'x ($\mu m$)')
    plt.ylabel(r"x' ($mrad$)")
    plt.title(species+' transverse phase space distribution')
    plt.xlim(params['xlim'])
    plt.ylim(params['ylim'])
    plt.grid(True)

    # Update the scatter plot data
    i = ind+1;
    def updatefig(*args):
        nonlocal i
        ptWeight, ptX, ptUx, colors = get_data(i)
        sct.set_offsets(np.stack((ptX, ptUx), axis=-1))
        sct.set_color(colors)
        i += 1
        # If we run over, loop
        if i == Nt+1:
            i = ind
        return sct,

    ani = animation.FuncAnimation(fig, updatefig, blit=True, frames=Nt-ind-2)
    ani.save(params['path']+'PhaseSpaceEvolutionEnergy_'+species+'.mp4',
             fps=params['fps'])


def emittance_energy(params):
    """ Plots the emittance growth and energy evolution of a given species.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            species : string
                The species name for the beam of interest.
            path : string
                The path to the VSim output folder.
            simName : strign
                The simulation name, the first part of every simulation file.
            Nt : int
                The number of dumps to use, set to the last dump number.
            mass : double
                The mass of the particle in GeV.
                
    Returns
    -------
    d : array-like
        The array of locations for the beam energy.
    energy : array-like
        The energy at each dump.
    """
    Nt = params['Nt']
    path = params['path']
    simName = params['simName']
    species = params['species']
    
    en = np.zeros(Nt+1) # normalized emittance
    energy = np.zeros(Nt+1)
    time = np.zeros(Nt+1)

    for i in range(0, Nt+1):
        pFile = get_filename(path, simName, species, i)
        pData = load.get_species_data(pFile, species)
        pAttrs = load.get_species_attrs(pFile, species)
        en[i] = analyze.get_normemittance(pData)
        energy[i] = analyze.get_energy(pData, params['mass'])
        time[i] = pAttrs['time']
    
    d = time*C*1e2
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    p1 = plt.plot(d, en, 'b-', label='Emittance')
    plt.title('Transverse emittance and energy evolution')
    plt.xlabel('Distance (cm)')
    plt.ylabel(r'Normalized emittance ($mm\cdot mrad$)')
    plt.twinx()
    p2 = plt.plot(d, energy, 'm-', label='Energy')
    plt.ylabel('Energy (GeV)')
    p = p1 + p2
    ax.legend(p, [l.get_label() for l in p])
    
    # Save the figure and display it
    plt.savefig(path+species+'emittanceEnergy.pdf', format='pdf')
    plt.savefig(path+species+'emittanceEnergy.png', format='png')
    plt.show()
    return d, energy


def energy_distribution(params):
    """ Plots the distribution of energy of a given species.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            species : string
                The species name for the beam of interest.
            path : string
                The path to the VSim output folder.
            simName : strign
                The simulation name, the first part of every simulation file.
            dumpInd : int
                The dump index to plot the beam at.
            mass : double
                The mass of the particle in GeV.
            bins : int
                The number of bins in the histogram.
    """
    path = params['path']
    simName = params['simName']
    species = params['species']
    ind = params['dumpInd']
    
    pFile = get_filename(path, simName, species, ind)
    pData = load.get_species_data(pFile, species)
    # pAttrs = load.get_species_attrs(pFile, species)
    
    weights = analyze.get_weights(pData)
    energy = analyze.get_ptc_energy(pData, params['mass'])
    
    plt.figure(figsize=(16, 9))
    plt.hist(energy, params['bins'], weights=weights)
    plt.title('Energy distribution of '+species)
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Counts')
    
    # Save the figure and display it
    plt.savefig(path+'EnergySpread_'+species+'_'+str(ind)+'.pdf', format='pdf')
    plt.savefig(path+'EnergySpread_'+species+'_'+str(ind)+'.png', format='png')
    plt.show()
