#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:14:24 2017

@author: robert
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.fft import fft, ifft, fftfreq, fftshift
from ionization import ionization
from ionization import adk
from propagation import propagation
from propagation import plasma
from scipy.interpolate import interp1d
from ht import intht
import os


def spectrum_from_axis(E, z):
    """ Returns the spatial spectrum in kz from the electric field on-axis.

    Returns the spatial frequencies, kz, and the spatial spectrum in terms of
    kz, of the electric field along the optical axis.

    Parameters
    ----------
    E : array_like
        Array of complex electric field values along the optical axis.
    z : array-like
        Array of z coordinates along the optical axis. Must be evenly spcaed
        for the FFT.

    Returns
    -------
    kz : array-like
        Array of spatial frequencies in z.
    S : array-like
        Spatial spectrum in terms of kz of the electromagnetic field.
    """
    N = np.size(E)
    dz = z[1] - z[0]
    S = fftshift(fft(E)) / N
    kz = 2*np.pi * fftshift(fftfreq(np.size(z), dz))
    return kz, S


def kr_from_kz(kz, lam, S=None):
    """ Returns the spatial frequencies in r from the frequencies in z.
    
    Parameters
    ----------
    kz : array-like
        Spatial frequencies in z.
    lam : double
        Wavelength of the electromagnetic wave in vacuum.
    S : array-like, optional
        The spatial spectrum expressed in terms of kz.

    Returns
    -------
    kr : array-like
        The spatial spectrum in terms of kr. Note may be a different length
        than kz, decaying (imaginary) frequencies are removed.
    S : array-like
        Spatial spectrum in terms of kr with imaginary frequency terms removed.
    """
    k = 2*np.pi / lam
    sel = kz <= k
    kr = np.sqrt(k**2 - kz[sel]**2)
    if S is not None:
        S = S[sel]
        return kr, S
    else:
        return kr


def uniform_bessel(params, Ez, z, n=0):
    """ Calculate the required electric field to create the passed intensity.

    Calculates the electric field necessary to create the passed intensity
    distribution along the optical axis. Unifrom Bessel means that it will
    assign a linear phase to the on axis electric field so that it exists in
    approximatly a single Bessel mode. This means the intensity distribution
    will have a constant width. 

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            N : int
                Number of integration steps in Bessel integrals. Try 1000.
            M : int
                Number of r values to calculate, each requires an integral.
            R : double
                Radius to the first zero of the Bessel function.
            lam : double
                Wavelength of the electromagnetic wave in vacuum.
            rmax : double
                Maximum radius to return the electric field at.
    Ez : array-like
        Desired on axis field in GV/m.
    z : array-like
        Array of z coordinates along the optical axis. Must be evenly spcaed
        for the FFT.
    n : int, optional
        Order of the Bessel function, defaults to order zero.

    Returns
    -------
    r : array-like
        Array of radius coordinates the electric field is given at.
    E : array-like
        Electric field as a function of r on the boundary.
    """
    #TODO figure out if the normalization is just not enough steps or an error
    lam = params['lam']
    k = 2*np.pi/lam
    kz, S = spectrum_from_axis(Ez, z)
    # Add the linear phase based of the width of the Bessel beam
    kr0 = 2.4048 / params['R']
    kz0 = np.sqrt(k**2 - kr0**2)
    kz = kz + kz0
    # Calculate the true spatial spectrum
    S = S / kz
    # Calculate the radial spatial frequencies
    kr, S = kr_from_kz(kz, lam, S)
    # Inverse Hankel transform
    krn = np.linspace(0, np.amax(kr), params['N'])
    Sn = interp1d(kr, S, fill_value=(S[-1], 0.0), bounds_error=False)
    Sn = Sn(krn)
    r = np.linspace(0, params['rmax'], params['M'])
    E = intht.ihtn(Sn, krn, r, n)
    return r, E


def multimode_ionization(params, z, I):
    """ Calculates the ionization fraction from multiple Bessel modes.

    This function calculates the ionization fraction resulting from a series
    of temporally seperated pulses. Each pulse has a a lateral intensity
    profile given by a Bessel function of order n. Modes should be ordered
    temporally.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            L : int
                Number of pulses to include.
            N : int
                Number of integration steps in Bessel integrals. Try 1000.
            M : int
                Number of r values to calculate, each requires an integral.
            R : array-like
                Radius to the first zero of the zero order Bessel function.
                Pass an array of length L with a value for each Bessel pulse.
            lam : double
                Wavelength of the electromagnetic wave in vacuum.
            rmax : array-like
                Maximum radius to return the electric field at. Pass an array
                of length L, one for each Bessel mode. Must be a factor of root
                2 larger than the grid in prop.
            prop : dictionary
                See the params dictionary for laser_prop for details.
            order : array-like
                Array of length L with the order of each Bessel function.
            path : string
                File loaction for storing simulation results.
            atom : dictionary
                See ionization for details and several examples.
            tau : double
                Temporal length of the Gaussian ionizing pulse.
            multi : array-like
                Multiplier for each orders electric field.
            xlim : array-like, optional
                Two element array of limits for the transverse density plot.
    z : array-like
        Array of on axis z values the intensity is specified at.
    I : array-like
        Desired Intensity profile along the optical axis.
    """
    L = params['L']
    order = params['order']
    prop = params['prop']
    atom = params['atom']
    Ez = ionization.field_from_intensity(I)
    besselParams = {'N': params['N'],
                    'M': params['M'],
                    'lam': params['lam']
                    }
    E = {} # Stores electric fields on the boundaries, input to laser_prop
    rmi = {}
    Efield = {} # Stores output electric fields
    frac = {}
    # Define the Efunc for laser_prop with Fourier series phi components
    def Efunc(x, y): # This is complicated because of the atan range
        r = np.sqrt(x**2 + y**2)
        Efield = prop['Efield']
        phi = np.zeros(np.shape(r)) 
        # Handle when x/y -> ininity
        phi[int(prop['Nx']/2), int(prop['Ny']/2):] = np.pi/2
        phi[int(prop['Nx']/2), :int(prop['Ny']/2)] = -np.pi/2
        # Handle the positive x half plane
        sel = np.array(x > 0)
        xp = x[sel]
        xp = np.reshape(xp, (np.size(xp), 1))
        phi[int(prop['Nx']/2+1):, :] = np.arctan(y/xp)
        # Handle the negative x half plane
        sel = np.array(x < 0)
        xn = x[sel]
        xn = np.reshape(xn, (np.size(xn), 1))
        phi[:int(prop['Nx']/2), :] = np.arctan(y/xn) + np.pi
        E0 = Efield(r) * np.exp(1j*prop['order']*phi)
        return E0

    for i in range(0, L):
        # Find the boundary electric field necessary to create the Bessel mode
        besselParams['R'] = params['R'][i]
        besselParams['rmax'] = params['rmax'][i]
        rmi[i], E[i] = uniform_bessel(besselParams, Ez, z, order[i])
        E[i] = 8.15e6*E[i] # Normalization factor I still need to fix
        E[i] *= params['multi'][i]
        # Propagate the modes to find the electric field
        prop['Efield'] = interp1d(rmi[i], E[i])
        prop['order'] = order[i]
        prop['path'] = params['path'] + 'Bessel_' + str(i) +'/'
        if not os.path.exists(prop['path']):
            os.makedirs(prop['path'])
        propagation.laser_prop(prop, Efunc)
        Efield[i] = np.load(prop['path'] + 'electricField.npy')
        # Find the ionization fraction and update the total ionization fraction
        frac[i] = adk.gaussian_frac(atom['EI'], abs(Efield[i]), params['tau'],
                                    atom['Z'], atom['l'], atom['m'])
        if i == 0:
            frac['tot'] = frac[i]
        else:
            frac['tot'] += (1- frac['tot'])*frac[i]
    path = params['path']
    np.save(path+'ionizationFraction', frac['tot'])
    np.save(path+'inputField', E)
    np.save(path+'inputCoordinates', rmi)
    np.save(path+'params', params)


def multimode_refraction(params, Tfunc):
    """ Calculates the ionization fraction including refraction.

    This function calculates the ionization fraction from multiple, time
    delayed, Bessel modes accounting for refractive effects as the plasma
    ionizes. Note that multimode_ionization should be run first to calculate
    the inital fields. The output can also be used to verify the input
    parameters are correct before running this time consuming function.

    Parameters
    ----------
    params : dictionary
        Params should have the following items:
            seq : array-like
                The sequence in which the modes propagate, i.e. [0, 1] means
                0 order first followed by first order.
            path : string
                Path to the output folder from multimode_ionization
            plasma : dictionary
                Additionaly parameters needed by plasma_refraction, see plasma_
                refraction for details about the elements
                    Nt : int
                        Number of temporal grid points.
                    T : double
                        Time duration.
                    n0 : double
                        Initial gas density in 10^17 cm^-3.
                    alpha : double
                        Atomic polarizability of the gas in A^3.
                    EI : double
                        Ionization energy in eV.
    """
    path = params['path']
    plasmaPar = params['plasma']
    directory = 'Density-' + str(plasmaPar['n0']) + '_' + \
                str(params['seq']) + '/' 
    if not os.path.exists(path+directory):
        os.makedirs(path+directory)
    
    # Define the initial electric field loader
    def Efunc(x, y):
        E0 = np.load(modePar['Esource']+'inputField.npy')
        return E0
    
    n = None # No initial ionization
    
    for i in params['seq']:
        modeName = 'Bessel_' + str(i) + '/'
        EsourcePath = path + modeName
        modePar = np.load(EsourcePath+'params.npy').item()
        for name in plasmaPar:
            modePar[name] = plasmaPar[name]
        modePar['Esource'] = EsourcePath
        modePar['nFinal'] = True
        modePar['path'] = path + directory + modeName
        if not os.path.exists(modePar['path']):
            os.makedirs(modePar['path'])
        # Run the refraction calculation
        plasma.plasma_refraction(modePar, Efunc, Tfunc, n=n)
        plasma.summary_plot(modePar['path'])
        n = np.load(modePar['path'] + 'finalDensity.npy')


def open_data(path, data=None):
    """ Opens the data files and sets up basic variables.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    data : string, optional
        Optional path to load the ionization fraction from. Used for plotting
        the output of a refraction calculation.
    """
    # Open the data files
    if data is None:
        frac = np.load(path+'ionizationFraction.npy')
    else:
        frac = np.load(data)
        frac = np.moveaxis(frac, 2, 0)
    E = np.load(path+'inputField.npy').item()
    rmi = np.load(path+'inputCoordinates.npy').item()
    params = np.load(path+'params.npy').item()
    # Useful simulation parameters
    L = params['L']
    X = params['prop']['X']
    Z = params['prop']['Z']
    Nx = params['prop']['Nx']
    Ny = params['prop']['Ny']
    Nz = params['prop']['Nz']
    return frac, E, rmi, params, L, X, Z, Nx, Ny, Nz


def ionization_plot(path, H, dr=None, data=None, suffix=''):
    """ Create plots of the ionization fraction and slice profiles.

    Specify the path to the output files and this function will save an image
    in the results directory that summarizes the simulation.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    H : int
        The number of lineouts to draw.
    dr : double, optional
        Spacing between radial lineouts.
    data : string, optional
        Optional path to load the ionization fraction from. Used for plotting
        the output of a refraction calculation.
    suffix : string, optional
        Suffix to add to the end of the filename.
    """
    frac, E, rmi, params, L, X, Z, Nx, Ny, Nz = open_data(path, data)
    x = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(0, Z, Nz)
    if 'xlim' in params:
        xlim = params['xlim']
    else:
        xlim = [-1e3, 1e3]
    if dr is None:
        dr = X/32
    # Create the figure
    gridSize = (2, 6)
    plt.figure(figsize=(16, 9))
    gridspec.GridSpec(gridSize[0], gridSize[1])

    # Input field plot
    ax1 = plt.subplot2grid(gridSize, (0, 0), colspan=2)
    ax1.set_prop_cycle('color',
                       [plt.cm.Set1(i) for i in np.linspace(0, 1, L)])
    for i in range(0, L):
        plt.plot(rmi[i]/1e3, abs(E[i]))
    plt.title('Input electric field for each order')
    plt.xlabel('r (mm)')
    plt.ylabel(r'$|E(r)|$ (GV/m)')
    plt.legend(['Bessel %d' % i for i in range(0, L)])

    # Ionization fraction plot
    plt.subplot2grid(gridSize, (0, 2), colspan=4)
    plt.imshow(np.flipud(np.transpose(frac[:, :, int(Ny/2)])),
           aspect='auto',
           extent=[0, Z/1e6, -X/2e3, X/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Ionization fraction')
    plt.set_cmap('plasma')
    plt.xlabel(r'z ($m$)')
    plt.ylabel(r'x ($mm$)')
    plt.title('Combined Ionization fraction')
    plt.xlim([0, Z/1e6])
    plt.ylim([-X/8e3, X/8e3])

    # Transverse density plot
    ax3 = plt.subplot2grid(gridSize, (1, 0), colspan=3)
    sliceInd = np.linspace(Nz/3, 2*Nz/3, H, dtype=np.int)
    dz = Z/(Nz-1)
    slicePosition = (sliceInd-1)*dz/1e6
    slicePosition = ['%.2f' % num for num in slicePosition]
    ax3.set_prop_cycle('color',
                      [plt.cm.gist_rainbow(i) for i in np.linspace(0, 1, H)])
    for i in range(0, H):
        plt.plot(x, frac[sliceInd[i], :, int(Ny/2)])
    plt.title('Transverse intensity profile at different distances')
    plt.xlabel('x ($\mu m$)')
    plt.ylabel('Ionization fraction')
    plt.legend(slicePosition, title='z in m')
    plt.xlim(xlim)
    plt.grid(True)

    # Longitudinal density plot
    ax4 = plt.subplot2grid(gridSize, (1, 3), colspan=3)
    dx = X/(Nx-2)
    sliceInd = np.linspace(Nx/2, Nx/2+dr*H/dx, H, dtype=np.int)
    slicePosition = (sliceInd-1)*dx - X/2
    slicePosition = ['%.2f' % num for num in slicePosition]
    ax4.set_prop_cycle('color',
                      [plt.cm.winter(i) for i in np.linspace(0, 1, H)])
    for i in range(0, H):
        plt.plot(z/1e6, frac[:, sliceInd[i], int(Ny/2)])
    plt.title('Longitudinal intensity profiles for different radiuses')
    plt.xlabel('z (m)')
    plt.ylabel('Ionization fraction')
    plt.legend(slicePosition, title='Radius in  $\mu m$')
    plt.grid(True)

    # Save the figure and display it
    plt.tight_layout()
    plt.savefig(path+'ionizationFig'+suffix+'.pdf', format='pdf')
    plt.savefig(path+'ionizationFig'+suffix+'.png', format='png')
    plt.show()
    
    # Close the file and clear the memory (fixes a bug in numpy)
    del frac
    del E
    del rmi
    del params
