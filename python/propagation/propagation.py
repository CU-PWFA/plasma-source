# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 17:05:30 2017

@author: robert
"""

import numpy as np
from propagation import laser
from ionization import ionization
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Useful prep function
def prep_data(data):
    return np.flipud(np.transpose(ionization.intensity_from_field(abs(data))))


def laser_prop(params, Efunc):
    """ Propagate a laser pukse through a region of uniform index.
    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            Nx : int
                Number of grid points in the x direction. Powers of 2 will run
                faster in the FFT.
            Ny : int
                Number of grid points in the y direction. Powers of 2 will run
                faster in the FFT.
            Nz : int
                Number of grid points in the z direction.
            X : double
                Width of the grid in x in um.
            Y : double
                Width of the grid in y in um.
            Z : double
                Length of the grid in z in um.
            z0 : double, optional
                Beginning of the Z grid. 
            E0 : double
                Peak elecric field strength. Acts as a multiplier to Efunc, if
                Efunc isn't normalized this can be set to 1.
            lam : double
                Wavelength of the optical pulse in um.
            n : double
                Index of refraction the light is propagating through.
            path : string
                File loaction for storing simulation results.
    Efunc : function
        Function that returns the initial E field at z=0, Efunc(x,y).
        Must be able to x as a  size (Nx, 1) numpy array and y as a
        size (1, Ny) numpy array and return the field in the form
        (x, y). The grid is centered on x = y = 0. Normally Efunc is
        normallized to 1 and the intensity comes from params[E0].
    """
    # Store some useful parameters as variables
    X = params['X']
    Y = params['Y']
    Z = params['Z']
    Nx = params['Nx']
    Ny = params['Ny']
    Nz = params['Nz']
    if 'z0' in params:
        z0 = params['z0']
    else:
        z0 = 0
    # Initialize the grid
    x = np.linspace(-X/2, X/2, Nx, False)
    y = np.linspace(-Y/2, Y/2, Ny, False)
    z = np.linspace(z0, z0+Z, Nz)
    # Setup an array to store the simulation results
    Efile = np.zeros((Nz, Nx, Ny))
    # Calculate the initial electric field on the boundary
    Ei = Efunc(np.reshape(x, (Nx, 1)), np.reshape(y, (1, Ny))) * params['E0']
    # Propogate the field
    Efile = laser.fourier_prop2(Ei, x, y, z, params['lam'], n=params['n'])
    # Save the data
    path = params['path']
    np.save(path+'electricField', Efile)
    np.save(path+'inputField', Ei)
    np.save(path+'params', params)


def laser_prop_plot(path):
    """ Creates a plot of the results of laser_prop.

    Specify the path to the output files and this function will save an image
    in the results directory with a slice along the x-z plane and a slice along
    the x-y plane at the beginning middle and end of the z grid.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    """
    # Open the data files
    Eplot = np.load(path+'electricField.npy')
    Ei = np.load(path+'inputField.npy')
    params = np.load(path+'params.npy').item()
    # Useful simulation parameters
    X = params['X']
    Y = params['Y']
    Z = params['Z']
    Nx = params['Nx']
    Ny = params['Ny']
    Nz = params['Nz']
    if 'z0' in params:
        z0 = params['z0']
    else:
        z0 = 0
    x = np.linspace(-X/2, X/2, Nx, False)

    # Create the figure
    gridSize = (2, 10)
    plt.figure(figsize=(16, 9))
    gridspec.GridSpec(gridSize[0], gridSize[1])

    # Initial electric field
    plt.subplot2grid(gridSize, (0, 0), colspan=2)
    plt.plot(abs(Ei[:, int(Ny/2)]), x/1e3, 'b-')
    plt.xlabel(r'E (GV/m)')
    plt.ylabel(r'x ($mm$)')
    plt.title('Transverse pulse shape')

    # X-Z intensity profile
    plt.subplot2grid(gridSize, (0, 2), colspan=8)
    plt.imshow(prep_data(Eplot[:, :, int(Ny/2)]),
               aspect='auto',
               extent=[z0/1e3, (z0+Z)/1e3, -X/2e3, X/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14}\,W/cm^3$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'z ($mm$)')
    plt.ylabel(r'x ($mm$)')
    plt.title('Intensity profile in the x-z plane')
    plt.xlim([z0/1e3, (z0+Z)/1e3])
    plt.ylim([-X/8e3, X/8e3])

    # X-Y initial intensity
    plt.subplot2grid(gridSize, (1, 1), colspan=3)
    plt.imshow(prep_data(Eplot[0, :, :]),
               aspect='auto',
               extent=[-X/2e3, X/2e3, -Y/2e3, Y/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14}\,W/cm^3$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'x ($mm$)')
    plt.ylabel(r'y ($mm$)')
    plt.title('Initial transverse intensity')
    plt.xlim([-X/2e3, X/2e3])
    plt.ylim([-Y/2e3, Y/2e3])

    # X_Y center intensity
    plt.subplot2grid(gridSize, (1, 4), colspan=3)
    plt.imshow(prep_data(Eplot[int(Nz/2), :, :]),
               aspect='auto',
               extent=[-X/2e3, X/2e3, -Y/2e3, Y/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14}\,W/cm^3$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'x ($mm$)')
    plt.ylabel(r'y ($mm$)')
    plt.title('Midpoint (Z/2) transverse intensity')
    plt.xlim([-X/4e3, X/4e3])
    plt.ylim([-Y/4e3, Y/4e3])

    # X_Y center intensity
    plt.subplot2grid(gridSize, (1, 7), colspan=3)
    plt.imshow(prep_data(Eplot[Nz-1, :, :]),
               aspect='auto',
               extent=[-X/2e3, X/2e3, -Y/2e3, Y/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14}\,W/cm^3$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'x ($mm$)')
    plt.ylabel(r'y ($mm$)')
    plt.title('Final transverse intensity')
    plt.xlim([-X/4e3, X/4e3])
    plt.ylim([-Y/4e3, Y/4e3])

    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'summaryFig.pdf', format='pdf')
    plt.savefig(path+'summaryFig.png', format='png')
    plt.show()

    # Close the file and clear the memory (fixes a bug in numpy)
    del Eplot
    del Ei
    del params
