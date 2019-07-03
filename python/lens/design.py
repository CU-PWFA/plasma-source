#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:04:21 2019

@author: robert
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq, fftshift
import matplotlib.pyplot as plt
from ionization import ionization
from lens import bessel
from propagation import laser
from beam.beams import laserbeam
from beam import interactions
from beam.beams import laserpulse
from beam.elements import plasma
import matplotlib.colors as colors


def calculate_tran_field(z, I, R, width, lam, dk=None, xlim=None, rlim=None):
    """ Calculate the transverse field from an desired on-axis intensity.
    
    Parameters
    ----------
    z : array of doubles
        Z positions corresponding to the intensity array (must be evenly spaced).
    I : array of doubles
        The desired on axis intensity distribution.
    R : double
        The maximum radius to calculate the field at, in um.
    width : double
        The width of the plasma in um, measured from the zeros of the Bessel function.
    lam : double
        The laser wavelength in um.
    dk : optional, double
        The range of spatial frequencies to plot.
    xlim : optional, tuple or array
        Bounds of z-axis plots.
    rlim : optional, tuple or array
        The limits for the transverse field plot
        
    Returns
    -------
    r : array of doubles
        Radius values for the electric field array.
    E : array of complex128
        The electric field at each radius.
    """
    # Calculate and display the kz spectrum
    #--------------------------------------------------------------------------
    params = {
        'Nr' : 5000,
        'R' : R,
        'w' : width,
        'lam' : lam
    }
    Nz = len(z)
    r = np.linspace(0, params['R'], params['Nr'])
    E = np.zeros(params['Nr'], dtype='complex128')
    Ez = ionization.field_from_intensity(I)
    k = 2*np.pi/lam
    dz = z[1] - z[0]
    # Shift frequencies
    kr0 = 2.4048 / params['w']
    kz0 = np.sqrt(k**2 - kr0**2)
    kz = 2*np.pi * fftshift(fftfreq(Nz, dz)) + kz0
    e = fftshift(fft(Ez)) / Nz
    
    plt.figure(figsize=(8, 2), dpi=150)
    plt.plot(kz, ionization.intensity_from_field(e))
    plt.xlabel(r'$k_z$ ($\mu m^{-1}$)')
    plt.ylabel(r'I ($\mathrm{10^{14}W/cm^2}$)')
    plt.xlim(k-dk, k)
    plt.show()
    
    # Calculate the required transverse field and show the actual longitudinal profile
    #--------------------------------------------------------------------------
    r, E = bessel.bessel_expansion(params, z, I)
    zFres = np.linspace(1e5, 8e6, 1000)
    eFres = laser.fresnel_axis(E, r, zFres, lam)
    IFres = ionization.intensity_from_field(eFres)
    
    plt.figure(figsize=(8, 2), dpi=150)
    plt.plot(z/1e6, I)
    plt.plot((zFres+z[0])/1e6, IFres, 'm')
    plt.legend(['Target intensity', 'Actual intensity'])
    plt.xlabel(r'z (m)')
    plt.ylabel(r'I ($\mathrm{10^{14}W/cm^2}$)')
    if xlim is not None:
        plt.xlim(xlim)
    plt.show()
    # Plot the transverse phase and intensity
    #--------------------------------------------------------------------------
    # Radial dependence of the phase and and intensity after the beam shaping optics
    if rlim is None:
        rlim = [0, R/1e3]
    plt.figure(figsize=(8, 2), dpi=150)
    plt.subplot(121)
    plt.plot(r/1e3, ionization.intensity_from_field(E))
    plt.xlabel(r'r (mm)')
    plt.ylabel(r'I ($\mathrm{10^{14}W/cm^2}$)')
    plt.xlim(rlim)
    plt.subplot(122)
    plt.plot(r/1e3, np.unwrap(np.angle(E)))
    plt.xlabel(r'r (mm)')
    plt.ylabel(r'$\phi$ (radians)')
    plt.xlim(rlim)
    plt.show()
    
    return r, E

def propagate_to_start(r, E, Z, X, Nx, path, lam, tau, threads, xlim=None):
    """ Create a beam from a transverse field and propogate to the plasma.
    
    Parameters
    ----------
    r : array of doubles
        The radius values for the field.
    E : array of complex128
        The complex electric field as a function of r.
    Z : double'
        Distance to the plasma in um.
    X : double
        Transverse domain size in um.
    Nx : int
        Number of cells in the transverse dimension.
    path : string
        Path to save the data at.
    lam : double
        The laser wavelength in um.
    tau : double
        The RMS pulse length in fs.
    xlim : optional, tuple or array
        Bounds for the transverse intensity plot after propagation
    threads : int
        The number of threads to run the calculation on.
        
    Returns
    -------
    beam : laserBeam object
        The beam object with the field information.
    pulseParams : dict
        The params object used to build the beam.
    """
    pulseParams = {
        'Nt' : 2**6,
        'Nx' : Nx,
        'Ny' : Nx,
        'X' : X,
        'Y' : X,
        'T' : 3*tau,
        'lam' : lam,
        'path' : path,
        'load' : False,
        'threads' : 4,
        'cyl' : True,
        'tau' : tau,
        'name' : 'To_Start',
    }
    beam = laserbeam.Laser(pulseParams)
    e = beam.reconstruct_from_cyl(r, E, beam.x, beam.y)
    beam.initialize_field(e)
    
    # Plot the initial transverse intensity
    #--------------------------------------------------------------------------
    e = beam.e
    I = beam.intensity_from_field(e)
    I = beam.prep_data(I)
    plt.figure(figsize=(8, 3), dpi=150)
    plt.subplot(121)
    plt.imshow(I, aspect='auto', extent=[-X/2e3, X/2e3, -X/2e3, X/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14} W/cm^2$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'x (mm)')
    plt.ylabel(r'y (mm)')
    plt.title('Peak transverse intensity at z=0.0cm')
    
    # Plot the initial transverse intensity
    #--------------------------------------------------------------------------
    beam.propagate(Z, 1.0)
    e = beam.e
    I = beam.intensity_from_field(e)
    I = beam.prep_data(I)
    plt.subplot(122)
    plt.imshow(I, aspect='auto', extent=[-X/2e3, X/2e3, -X/2e3, X/2e3])
    cb = plt.colorbar()
    cb.set_label(r'Intensity ($10^{14} W/cm^2$)')
    plt.set_cmap('viridis')
    plt.xlabel(r'x (mm)')
    plt.ylabel(r'y (mm)')
    plt.title('Peak transverse intensity at z=%0.1fcm' % (Z/1e4))
    if xlim is not None:
        plt.xlim(xlim)
        plt.ylim(xlim)
    plt.tight_layout()
    plt.show()
    
    return beam, pulseParams

def domain_test(X, Nx, Z, Nz, beam0, pulseParams, z_target, I_target, start, ylim=None, log=False):
    """ Propagate the beam to see if the domain is large enough.
    
    Parameters
    ----------
    X : double
        Transverse domain size in um.
    Nx : int
        Number of cells in the transverse dimension.
    Z : double
        Distance to propogate.
    Nz : int
        Number of z steps in the propogation.
    beam0 : laserBeam object
        The beam object to pull the field from.
    pulseParams : dict
        The pulse params from the beam0 object.
    z_target : array of doubles
        Z array for the target intensity.
    I_target : array of doubles
        Target intensity to compare the actual against.
    start : double
        The distance to shift the targte intensity pattern.
    ylim : optional, array or tuple
        The y limits of the plot.
    log : optional, bool
        Plot the intensity on a log scale
    """
    z = np.linspace(0, Z, Nz)
    pulseParams['name'] = 'Test_Beam'
    pulseParams['Nx'] = Nx
    pulseParams['Ny'] = Nx
    pulseParams['X'] = X
    pulseParams['Y'] = X
    beam1 = laserbeam.Laser(pulseParams)
    e = beam1.reconstruct_from_cyl(beam0.x, beam0.e[:, int(beam0.Ny/2)], beam1.x, beam1.y)
    beam1.initialize_field(e)
    beam1.propagate(z, 1.0)
    
    e1 = np.zeros((Nz, Nx), dtype='complex128')
    for i in range(Nz):
        e1[i, :] = beam1.load_field(i+1)[0]
    I = ionization.intensity_from_field(e1)
    I_max = np.amax(I)
    
    ext = [0, Z/1e4, -X/2, X/2]
    plt.figure(figsize=(8, 2), dpi=150)
    if log:
        norm = colors.LogNorm(vmin=I_max*1e-4, vmax=I_max)
        plt.imshow(np.flipud(np.transpose(I)), aspect='auto', extent=ext, cmap='viridis', norm=norm)
    else:
        plt.imshow(np.flipud(np.transpose(I)), aspect='auto', extent=ext, cmap='viridis')
    cb = plt.colorbar()
    cb.set_label(r'Laser Intensity ($10^{14} W/cm^2$)')
    plt.xlabel('z (cm)')
    plt.ylabel(r'x ($\mathrm{\mu m}$)')
    if ylim is not None:
        plt.ylim(ylim)
        
    plt.twinx()
    plt.plot((z_target-start)/1e4, I_target, 'w-', label='Target')
    plt.plot(np.array(beam1.z[:-1])/1e4, I[:, int(Nx/2)], 'r--', label='Simulated')
    plt.legend(loc=8)
    plt.xlim(0, Z/1e4)
    plt.show()
    
def plasma_refraction(X, Nx, Z, Nz, beam0, pulseParams, species, n, start, m_e):
    """ Propagate the laser pulse through the gas profile.
    
    Parameters
    ----------
    X : double
        Transverse domain size in um.
    Nx : int
        Number of cells in the transverse dimension.
    Z : double
        Distance to propogate.
    Nz : int
        Number of z steps in the propogation.
    beam0 : laserBeam object
        The beam object to pull the field from.
    pulseParams : dict
        The pulse params from the beam0 object.
    species : dict
        Ionization species dict from ionization.ionization.
    n : func
        Function that returns the gas density as a function of z.
    m_e : double
        Multiplier for the intensity.
    """
    ne0 = 3.4e16/1e17 #This doesn't do anything but we need it or else the plasma
    # class throws an error, plasma should have the parameter removed
    pulseParams['name'] = 'Refracted_Beam'
    pulseParams['Nx'] = Nx
    pulseParams['Ny'] = Nx
    pulseParams['X'] = X
    pulseParams['Y'] = X
    plasmaParams = {
        'Nx' : Nx,
        'Ny' : Nx,
        'Nz' : Nz,
        'X' : X,
        'Y' : X,
        'Z' : Z,
        'atom' : species,
        'path' : pulseParams['path'],
        'load' : False,
        'cyl' : True,
        'name' : 'Plasma_Source',
        'n0' : ne0
    }
    tau = pulseParams['tau']
    pulse = laserpulse.Pulse(pulseParams)
    e = np.sqrt(m_e)*pulse.reconstruct_from_cyl(beam0.x, np.array(beam0.e)[:, int(beam0.Ny/2)], pulse.x, pulse.y)
    e = e[None, :, :]*np.exp(-pulse.t[:, None, None]**2*np.pi/(2*tau**2))
    pulse.initialize_field(e)
    plasma_source = plasma.Plasma(plasmaParams)

    # Initialize gas density
    n_gas = np.zeros((Nx, Nx, Nz), dtype='double')
    ne = np.zeros((Nx, Nx, Nz), dtype='double')
    for i in range(Nz):
        n_gas[:, :, i] = n(plasma_source.z[i]+start)
    plasma_source.initialize_plasma(n_gas, ne)
    # Propagate the laser beam through the oven
    interactions.pulse_plasma(pulse, plasma_source)
    e = np.zeros((Nz, Nx), dtype='complex128')
    ne = np.zeros((Nz, Nx))
    for i in range(0, Nz-1):
        ne[i, :] = plasma_source.load_plasma_density(i)[0]
    for i in range(Nz):
        e[i, :] = pulse.load_field(i)[0][int(pulseParams['Nt']/2), :]
    I = ionization.intensity_from_field(e)
    ne = ne*1e17
    return pulse, I, ne

def plot_laser_plasma(I, ne, ext):
    plt.figure(figsize=(16, 4), dpi=150)
    plt.subplot(121)
    plt.imshow(np.flipud(np.transpose(I)), aspect='auto', extent=ext, cmap='viridis')
    cb = plt.colorbar()
    cb.set_label(r'Laser Intensity ($10^{14} W/cm^2$)')
    plt.xlabel('z (cm)')
    plt.ylabel(r'x ($\mathrm{\mu m}$)')
    plt.ylim(-500, 500)

    plt.subplot(122)
    plt.imshow(np.flipud(np.transpose(ne)), aspect='auto', extent=ext, cmap='plasma')
    cb = plt.colorbar()
    cb.set_label(r'$n_e$ ($\mathrm{cm^{-3}}$)')
    plt.xlabel('$z$ (cm)')
    plt.ylabel(r'$x$ ($\mathrm{\mu m}$)')
    plt.ylim(-500, 500)
    plt.tight_layout()
    plt.show()
    
def pulse_evolution():
    pass
    