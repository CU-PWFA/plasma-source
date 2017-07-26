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
                of length L, one for each Bessel mode.
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
        # Propagate the modes to find the electric field
        prop['Efield'] = interp1d(rmi[i], E[i])
        prop['order'] = order[i]
        prop['path'] = params['path'] + 'Bessel_' + str(i) +'/'
        if not os.path.exists(prop['path']):
            os.makedirs(prop['path'])
        propagation.laser_prop(prop, Efunc)
        Efield[i] = np.load(prop['path'] + 'electricField.npy')
        Efield[i] *= params['multi'][i]
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


def open_data(path):
    """ Opens the data files and sets up basic variables.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    """
    # Open the data files
    frac = np.load(path+'ionizationFraction.npy')
    E = np.load(path+'inputField.npy')
    rmi = np.load(path+'inputCoordinates.npy')
    params = np.load(path+'params.npy')
    print(params)
    # Useful simulation parameters
    L = params['L']
    X = params['prop']['X']
    Z = params['prop']['Z']
    Nx = params['prop']['Nx']
    Ny = params['prop']['Ny']
    Nz = params['prop']['Nz']
    return frac, E, rmi, params, L, X, Z, Nx, Ny, Nz


def ionization_plot(path):
    """ Create plots of the ionization fraction and slice profiles.

    Specify the path to the output files and this function will save an image
    in the results directory that summarizes the simulation.

    Parameters
    ----------
    path : string
        The path specifying the directory with the simulation results.
    """
    frac, E, rmi, params, L, X, Z, Nx, Ny, Nz = open_data(path)
    
    # Create the figure
    gridSize = (2, 6)
    plt.figure(figsize=(16, 9))
    gridspec.GridSpec(gridSize[0], gridSize[1])

    # Input field plot
    plt.subplot2grid(gridSize, (0, 0), colspan=2)
    ax = plt.axes()
    ax.set_prop_cycle('color', [plt.cm.cool(i) for i in np.linspace(0, 1, L)])
    for i in range(0, L):
        plt.plot(rmi[i], E[i])
    plt.title('Input electric field for each order')
    plt.xlabel('r (m)')
    plt.ylabel(r'$|E(r)|$ (GV/m)')

    # Ionization fraction plot
    plt.subplot2grid(gridSize, (0, 1), colspan=4)
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
    
    plt.tight_layout()
    plt.savefig(path+'ionizationFig.pdf', format='pdf')
    plt.savefig(path+'ionizationFig.png', format='png')
    plt.show()
    
    # Close the file and clear the memory (fixes a bug in numpy)
    del frac
    del E
    del rmi
    del params
