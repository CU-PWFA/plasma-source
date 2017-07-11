# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:18:14 2017

@author: rariniello
"""

import numpy as np
from propagation import laser
from ionization import adk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from propagation import propagation


def plasma_refraction(params, Efunc, Tfunc, n0=False):
    """ Propagate a laser pulse through a plasma accounting for refraction.

    Propogates a laser pulse through a region of partially ionized gas. This
    function accounts for refraction from the plasma. It determines the plasma
    density by calculating the ionization that has resulted from each temporal
    piece of the pulse. The results are stored in a file, only the central x-z
    plane is recorded.

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
            Nt : int
                Number of pieces the pulse is discretized into.
            X : double
                Width of the grid in x in um.
            Y : double
                Width of the grid in y in um.
            Z : double
                Length of the grid in z in um.
            T : double
                Duration of the pulse discretization, in fs.
            n0 : double
                Initial gas density in 10^17 cm^-3. Can also be passed as an
                array argument for non-uniform plasmas.
            E0 : double
                Peak elecric field strength. Acts as a multiplier to Efunc, if
                Efunc and Tfunc aren't normalized this can be set to 1.
            alpha : double
                Atomic polarizability of the gas in A^3.
            EI : double
                Ionization energy of the first electron in the gas in eV.
            lam : double
                Wavelength of the optical pulse in um.
            path : string
                File loaction for storing simulation results.
            nFinal : bool
                Save the final 3D density profile.
    Efunc : function
        Function that returns the initial E field at z=0, Efunc(x,y).
        Must be able to x as a  size (Nx, 1) numpy array and y as a
        size (1, Ny) numpy array and return the field in the form
        (x, y). The grid is centered on x = y = 0. Normally Efunc is
        normallized to 1 and the intensity comes from params[E0].
    Tfunc : function
        Function that returns the pulse envelope as a function of time.
        The pulse is centered on t = 0. Normally Tfunc is normalized to 1.
    n0 : array_like, optional
        Can pass an (Nx, Ny, Nz) array of gas density. If absent, a uniform gas
        density of params['n0'] will be assumed. Units of 10^17 cm^-3.
    """
    # Store some useful parameters as variables
    X = params['X']
    Y = params['Y']
    Z = params['Z']
    T = params['T']
    Nx = params['Nx']
    Ny = params['Ny']
    Nz = params['Nz']
    Nt = params['Nt']
    if 'z0' in params:
        z0 = params['z0']
    else:
        z0 = 0
    lam = params['lam']
    # Initialize the grid
    x = create_centered_grid(X, Nx)
    y = create_centered_grid(Y, Ny)
    z = create_z_grid(Z, z0, Nz)
    t = create_centered_grid(T, Nt)
    dt = T/(Nt-1)
    # Setup the index of refraction array and plasma density array
    nih = np.zeros((Nx, Ny, Nz))
    n = np.zeros((Nx, Ny, Nz))
    # Handle non uniform initial densities
    if n0 is False:
        n0 = params['n0']
    # Setup arrays to store the simulation results
    Efile = np.zeros((Nt, Nz, Nx))
    nfile = np.zeros((Nt, Nx, Nz))
    # Calculate the intial field and the pulse
    # TODO add in a valid_params function that tests the params object
    Ei = Efunc(np.reshape(x, (Nx, 1)), np.reshape(y, (1, Ny))) * params['E0']
    Et = Tfunc(t)
    # Calculate constants for gas and plasma index of refraction
    # These are (index of refraction - 1) per 10^17 cm^-3 density
    ngas = params['alpha'] * 5e-8
    nplasma = plasma_index(1, lam) - 1
    nh = 1 + ngas*params['n0']
    # Main loop through the pulse
    for i in range(Nt):
        E = Ei * Et[i]
        # Propagate the beam through the plasma
        Ef = laser.beam_prop2(E, nih, x, y, z, lam, nh)
        Efile[i, :, :] = abs(Ef[:, :, int(Ny/2)])
        # Calculate ionization rate from the beam
        rate = adk.rate_linear(params['EI'], abs(Ef), 1)
        rate = np.moveaxis(rate, 0, -1)
        # Plasma density = Un-ionized gas density * fraction ionized + ionized
        # gas density
        n = (n0 - n)*(1 - np.exp(-rate * dt)) + n
        nfile[i, :, :] = n[:, int(Ny/2), :]
        # Calculate new index of refraction
        nih = n * nplasma + (n0 - n - params['n0']) * ngas
        print('Completed time slice ', i+1, ' of ', Nt)
    # Write params and  results to file
    path = params['path']
    np.save(path+'electricField', Efile)
    np.save(path+'ionizationFrac', nfile)
    np.save(path+'inputField', Ei)
    np.save(path+'inputPulse', Et)
    np.save(path+'params', params)
    if 'nFinal' in params and 'nFinal':
        np.save(path+'finalDensity', n)


def create_centered_grid(X, N):
    """ Create numpy array for the grid of width X, centered on 0.

    Parameters
    ----------
    X : double
        Width of the grid, -X/2 to X/2 not including X/2.
    N : int
        Number of grid points.

    Returns
    -------
    x : array-like
        The numpy array of grid locations in x.
    """
    return np.linspace(-X/2, X/2, N, False)


def create_z_grid(Z, z0, N):
    """ Create numpy array for the grid of width Z, starting at z0.

    Parameters
    ----------
    Z : double
        Width of the grid, z0 to Z inclusive.
    z0 : double
        Starting point of the grid.
    N : int
        Number of grid points.

    Returns
    -------
    z : array-like
        The numpy array of grid locations in z.
    """
    return np.linspace(z0, z0+Z, N)


def plasma_index(n, lam):
    """ Calculates the index of refraction of a plasma.

    Parameters
    ----------
    n : double
        Density of the plasma in 10^17 cm^-3.
    lam : double
        Wavelength of the incident light in um.
    """
    return 1 - n * lam**2 * 4.47869e-5


def open_data(path):
    """ Opens the data files and sets up basic variables.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    """
    # Open the data files
    Eplot = np.load(path+'electricField.npy')
    nplot = np.load(path+'ionizationFrac.npy')
    Ei = np.load(path+'inputField.npy')
    Et = np.load(path+'inputPulse.npy')
    params = np.load(path+'params.npy').item()
    # Useful simulation parameters
    X = params['X']
    Z = params['Z']
    T = params['T']
    Nx = params['Nx']
    Ny = params['Ny']
    Nz = params['Nz']
    Nt = params['Nt']
    if 'z0' in params:
        z0 = params['z0']
    else:
        z0 = 0
    return Eplot, nplot, Ei, Et, params, X, Z, T, Nx, Ny, Nz, Nt, z0


def summary_plot(path):
    """ Creates a summary plot of the results of plasma_refraction.

    Specify the path to the output files and this function will save an image
    in the results directory that summarizes the simulation.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    """
    Eplot, nplot, Ei, Et, params, X, Z, T, Nx, Ny, Nz, Nt, z0 = open_data(path)
    x = create_centered_grid(X, Nx)
    t = create_centered_grid(T, Nt)

    # Create the figure
    gridSize = (2, 5)
    plt.figure(figsize=(16, 9))
    gridspec.GridSpec(gridSize[0], gridSize[1])

    # Initial electric field
    plt.subplot2grid(gridSize, (0, 0))
    plt.plot(abs(Ei[:, int(Ny/2)]), x/1e3, 'b-')
    plt.xlabel(r'E (GV/m)')
    plt.ylabel(r'x ($mm$)')
    plt.title('Transverse pulse shape')

    # Final intensity profile
    plt.subplot2grid(gridSize, (0, 1), colspan=4)
    plt.imshow(propagation.prep_data(Eplot[int(Nt/2), :, :]),
               aspect='auto',
               extent=[z0/1e6, (z0+Z)/1e6, -X/2e3, X/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14}\,W/cm^3$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'z ($m$)')
    plt.ylabel(r'x ($mm$)')
    plt.title('Intensity profile at $t=0$ (pulse center)')
    plt.xlim([z0/1e6, (z0+Z)/1e6])
    plt.ylim([-X/4e3, X/4e3])

    # Temporal pulse shape
    plt.subplot2grid(gridSize, (1, 0))
    plt.plot(t, Et, 'b-')
    plt.xlabel(r't ($fs$)')
    plt.ylabel(r'E ($GV/m$)')
    plt.title('Temporal pulse shape')

    # Final plasma density
    plt.subplot2grid(gridSize, (1, 1), colspan=4)
    plt.imshow(np.flipud(nplot[Nt-1, :, :]),
               aspect='auto',
               extent=[z0/1e6, (z0+Z)/1e6, -X/2e3, X/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Plasma density ($10^{17}\,cm^{-3}$)')
    plt.set_cmap('plasma')
    plt.xlabel(r'z ($m$)')
    plt.ylabel(r'x ($mm$)')
    plt.title('Density of the laser ionized plasma')
    plt.xlim([z0/1e6, (z0+Z)/1e6])
    plt.ylim([-X/4e3, X/4e3])
    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'summaryFig.pdf', format='pdf')
    plt.savefig(path+'summaryFig.png', format='png')
    plt.show()

    # Close the file and clear the memory (fixes a bug in numpy)
    del Eplot
    del nplot
    del Ei
    del Et
    del params
    
def profile_plot(path, xlim=None):
    """ Creates plasma density plots through different profiles.

    Specify the path to the output files and this function will save an image
    in the results directory that shows how the plasma density varies
    longitudinally and transversly.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    xlim : array-like, optional
        Two element array of limits for x.
    """
    Eplot, nplot, Ei, Et, params, X, Z, T, Nx, Ny, Nz, Nt, z0 = open_data(path)
    x = create_centered_grid(X, Nx)
    z = create_z_grid(Z, z0, Nz)
    
    if xlim is None:
        xlim = [-X/4, X/4]

    # Create the figure
    gridSize = (2, 1)
    plt.figure(figsize=(16, 9))
    gridspec.GridSpec(gridSize[0], gridSize[1])
    
    # On-axis plasma density
    plt.subplot2grid(gridSize, (0, 0))
    plt.plot(z/1e6, nplot[Nt-1, int(Nx/2), :], 'b-')
    plt.xlabel(r'z ($m$)')
    plt.ylabel(r'Plasma density ($10^{17}\,cm^{-3}$)')
    plt.title('On-axis plasma density profile')
    # Transverse plasma density
    plt.subplot2grid(gridSize, (1, 0))
    plt.plot(x, nplot[Nt-1, :, int(Nz/2)], 'b-')
    plt.xlabel(r'x ($\mu m$)')
    plt.ylabel(r'Plasma density ($10^{17}\,cm^{-3}$)')
    plt.title('Transverse plasma density profile')
    plt.xlim(xlim)
    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'profileFig.pdf', format='pdf')
    plt.savefig(path+'profileFig.png', format='png')
    plt.show()

    # Close the file and clear the memory (fixes a bug in numpy)
    del Eplot
    del nplot
    del Ei
    del Et
    del params
