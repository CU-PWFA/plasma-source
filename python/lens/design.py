#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:04:21 2019

@author: robert
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq, fftshift
from scipy.integrate import simps
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from ionization import ionization
from lens import bessel
from lens import ray
from propagation import laser
from beam.beams import laserbeam
from beam import interactions
from beam.beams import laserpulse
from beam.elements import plasma
from beam.elements import optic
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
plt.rcParams['animation.ffmpeg_path'] = '/home/robert/anaconda3/envs/CU-PWFA/bin/ffmpeg'
import matplotlib.animation as animation


def calculate_tran_field(z, I, R, width, lam, path, dk=None, xlim=None, rlim=None):
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
    plt.tight_layout()
    plt.show()
    
    np.save(path+'r.npy', r)
    np.save(path+'e.npy', E)
    return r, E

def propagate_to_start(r, E, Z, X, Nx, path, lam, tau, threads, xlim=None, plot=True):
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
    plot : bool
        Show plots or not.
        
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
        'threads' : threads,
        'cyl' : True,
        'tau' : tau,
        'name' : 'To_Start',
    }
    beam = laserbeam.Laser(pulseParams)
    e = beam.reconstruct_from_cyl(r, E, beam.x, beam.y)
    beam.initialize_field(e)
    
    # Plot the initial transverse intensity
    #--------------------------------------------------------------------------
    if plot:
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
    if plot:
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

def domain_test(X, Nx, Z, Nz, beam0, pulseParams, z_target, I_target, start, ylim=None, log=False, plot=True, legend=True):
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
        Plot the intensity on a log scale.
    plot : bool
        Show plots or not.
        
    Returns
    -------
    I : array of doubles
        The x-z intensity from the simulation.
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
    if plot:
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

        dz = z[1]-z[0]
        plt.twinx()
        plt.plot((z_target-start-dz)/1e4, I_target, 'w-', label='Target')
        plt.plot(np.array(beam1.z[:-1])/1e4, I[:, int(Nx/2)], 'c--', label='Simulated')
        if legend:
            plt.legend(loc=8)
        plt.xlim(0, Z/1e4)
        plt.show()
    
    return I
    
def plasma_refraction(X, Nx, Z, Nz, beam0, pulseParams, species, n, start, m_e, ne0, name=None, t=0.0, n2=0.0, 
                      ionization_type='adk', multispecies= False, species2= None):
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
    name : string, optional
        Name for the pulse, defaults to Refracted_Beam.
    t : double, optional
        Plasma heating energy to remove from the field, in eV.
    n2 : double, optional
        The nonlinear index of refraction at atmospheric pressure. In cm^2/W.
    ionization_type : string, optional
        Function to use for the ionization model.
    multispecies : boolean
        turn multispecies on or not
    species2: dict
        Ionization species dict from ionization.ionization. 
        If multispecies is on, second sepcies is specified here
    """
    if multispecies == False:        
        if name is None:
            name = 'Refracted_Beam'
        pulseParams['name'] = name
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
        print('Initial pulse energy %0.2fmJ' % (pulse.pulse_energy()*1e3))
        plasma_source = plasma.Plasma(plasmaParams)
    
        # Initialize gas density
        n_gas = np.zeros((Nx, Nx, Nz), dtype='double')
        ne = np.zeros((Nx, Nx, Nz), dtype='double')
        for i in range(Nz):
            n_gas[:, :, i] = n(plasma_source.z[i]+start)
        plasma_source.initialize_plasma(n_gas, ne)
        # Propagate the laser beam through the oven
        interactions.pulse_plasma_energy(pulse, plasma_source, temp=t, n2=n2, ionization=ionization_type)
        print('Final pulse energy %0.2fmJ' % (pulse.pulse_energy()*1e3))
        e = np.zeros((Nz, Nx), dtype='complex128')
        ne = np.zeros((Nz, Nx))
        for i in range(0, Nz-1):
            ne[i, :] = plasma_source.load_plasma_density(i)[0]
        for i in range(Nz):
            e[i, :] = pulse.load_field(i)[0][int(pulseParams['Nt']/2), :]
        I = ionization.intensity_from_field(e)
        ne = ne*1e17
        return pulse, I, ne
    
    else:
#TODO check if species2 is none, if so, through error
        if name is None:
            name = 'Refracted_Beam'
        pulseParams['name'] = name
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
        plasmaParams2 = {
            'Nx' : Nx,
            'Ny' : Nx,
            'Nz' : Nz,
            'X' : X,
            'Y' : X,
            'Z' : Z,
            'atom' : species2,
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
        print('Initial pulse energy %0.2fmJ' % (pulse.pulse_energy()*1e3))
        plasma_source = plasma.Plasma(plasmaParams)
        plasma_source2 = plasma.Plasma(plasmaParams2)
    
        # Initialize gas density
        n_gas = np.zeros((Nx, Nx, Nz), dtype='double')
        ne = np.zeros((Nx, Nx, Nz), dtype='double')
        for i in range(Nz):
            n_gas[:, :, i] = n(plasma_source.z[i]+start)
        plasma_source.initialize_plasma(n_gas, ne)
        # Propagate the laser beam through the oven
        interactions.pulse_plasma_energy_second(pulse, plasma_source, plasma_source2, temp=t, n2=n2, ionization=ionization_type)
        print('Final pulse energy %0.2fmJ' % (pulse.pulse_energy()*1e3))
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
    """ Plot the laser intensity at z=0 and the final plasma density.
    
    Parameters
    ----------
    I : array of doubles
        The x-z slice of the intensity array.
    ne : array of doubles
        The x-z slice of the plasma density array.
    ext : array of doubles
        The extent of the image for imshow extent.
    """
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

def plot_plasma_density(pulse, ne, ne0, ext, lines=[20, 40, 60], name=None, xlim=None, ylim=None, xlim2=None, yticks=None):
    """ Plot the plasma desnity with line outs.

    Parameters
    ----------
    pusle : Pulse Object
        The laser pulse object that created the plasma (used for bounds).
    ne : array of doubles
        The x-z slice of the plasma density array.
    ne0 : double
        Flattop plasma density measured in 1e17 cm^-3.
    ext : array of doubles
        The extent of the image for imshow extent.
    lines : array of doubles
        The locations of the transverse lineouts in cm.
    name : string, optional
        Name and path for the image to be saved to.
    """
    if xlim is None:
        xlim = [ext[0], ext[1]]
    if ylim is None:
        ylim = [-500, 500]
    if xlim2 is None:
        xlim2 = [-350, 350]
    if yticks is None:
        yticks = [-200, 0, 200]
    orange = '#EE7733'
    fig = plt.figure(figsize=(10, 7), dpi=150)
    gs = gridspec.GridSpec(4, 2, height_ratios=(1.5, 0.5, 3, 1), width_ratios=(40, 1),
                           hspace=0)
    ax1 = plt.subplot(gs[2, 0])
    ne_im = np.flipud(np.transpose(ne/1e16))
    im = plt.imshow(ne_im, aspect='auto', extent=ext, cmap='plasma',
                    interpolation='Spline16')
    plt.ylabel(r'$x$ ($\mathrm{\mu m}$)')
    plt.ylim(ylim)
    grey2 = '#AAAAAA'
    linewidth = 0.8
    for i in lines:
        plt.plot([i, i], [ylim[0], ylim[1]], '-', c=grey2, linewidth=0.8)

    ax2 = plt.subplot(gs[2, 1])
    cb = plt.colorbar(im, cax=ax2)
    cb.set_label(r'$n_e$ ($10^{16}\mathrm{cm^{-3}}$)')

    ax3 = plt.subplot(gs[3, 0], sharex=ax1)
    plt.plot(np.array(pulse.z)/1e4, ne[:, int(pulse.Nx/2)]/1e16, c=orange)
    plt.xlim(xlim)
    plt.ylim(0, 1.1*ne0*10)
    plt.xlabel(r'$z$ (cm)')
    plt.ylabel(r'$n_e$')
    plt.grid(True)

    plt.tight_layout()
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    Nz = np.shape(ne)[0]
    Z = ext[1]
    gs2 = gridspec.GridSpecFromSubplotSpec(1, 3, wspace=0.0, subplot_spec=gs[0, 0])
    ax01 = plt.subplot(gs2[0, 2])
    plt.plot(ne[int(Nz*lines[2]/Z), :]/1e16, pulse.x, c=orange)
    plt.plot([ne0*10, ne0*10], xlim2, '-', c=grey2, linewidth=linewidth)
    plt.grid(True, axis='y')
    plt.xlim(0, 1.1*ne0*10)
    plt.xlabel(r'$n_e$')

    ax02 = plt.subplot(gs2[0, 1], sharey=ax01)
    plt.plot(ne[int(Nz*lines[1]/Z), :]/1e16, pulse.x, c=orange)
    plt.plot([ne0*10, ne0*10], xlim2, '-', c=grey2, linewidth=linewidth)
    plt.grid(True, axis='y')
    plt.xlim(0, 1.1*ne0*10)
    plt.xlabel(r'$n_e$')

    ax03 = plt.subplot(gs2[0, 0], sharey=ax02)
    plt.plot(ne[int(Nz*lines[0]/Z), :]/1e16, pulse.x, c=orange)
    plt.plot([ne0*10, ne0*10], xlim2, '-', c=grey2, linewidth=linewidth)
    plt.grid(True, axis='y')
    plt.xlim(0, 1.1*ne0*10)
    plt.xlabel(r'$n_e$')
    plt.ylim(xlim2)
    plt.ylabel(r'$x$ ($\mathrm{\mu m}$)')
    plt.yticks(yticks)

    plt.setp([a.get_yticklabels() for a in fig.axes[3:-1]], visible=False)
    if name is not None:
        plt.savefig(name+'.png')
    plt.show()
    
def plot_pulse(pulse, ind, ylim=None, log=False, smooth=False):
    """ Plot the pulse intensity in t-x space.
    
    Parameters
    ----------
    pulse : Pulse object
        The laser pulse to plot the evolution for.
    ind : int
        The z index to plot the pulse at.
    ylim : optional, array or tuple
        The y limits of the plot.
    log : optional, bool
        Plot the intensity on a log scale.
    """
    Nt = pulse.Nt
    T = pulse.T
    Nx = pulse.Nx
    X = pulse.X
    e = np.zeros((Nt, Nx), dtype='complex128')
    e[:, :] = pulse.load_field(ind)[0]
    I = ionization.intensity_from_field(e)
    I_max = np.amax(I)
    
    ext = [-T/2, T/2, -X/2, X/2]
    plt.figure(figsize=(4, 2), dpi=150)
    if smooth == True:
        interp = 'Spline16'
    else:
        interp = None
    if log:
        norm = colors.LogNorm(vmin=I_max*1e-4, vmax=I_max)
        plt.imshow(np.fliplr(np.flipud(np.transpose(I))), aspect='auto', extent=ext, cmap='viridis', norm=norm, interpolation=interp)
    else:
        plt.imshow(np.fliplr(np.flipud(np.transpose(I))), aspect='auto', extent=ext, cmap='viridis', interpolation=interp)
    cb = plt.colorbar(format="%0.2f")
    cb.set_label(r'Laser Intensity ($10^{14} W/cm^2$)')
    plt.xlabel('t (fs)')
    plt.ylabel(r'x ($\mathrm{\mu m}$)')
    if ylim is not None:
        plt.ylim(ylim)
    plt.show()
    
def pulse_evolution(pulse, name, ylim=None, log=False, smooth=False):
    """ Create an animation of the pulse evolution.
    
    Parameters
    ----------
    pulse : Pulse object
        The laser pulse to plot the evolution for.
    name : string
        Filename of the animation.
    ylim : optional, array or tuple
        The y limits of the plot.
    log : optional, bool
        Plot the intensity on a log scale.
    smooth : optional, bool
        Whether to use a Spline 16 interpolation during rendering.
    """
    Nt = pulse.Nt
    T = pulse.T
    Nx = pulse.Nx
    X = pulse.X
    Nz = len(pulse.z)
    e = np.zeros((Nt, Nx), dtype='complex128')
    e[:, :], z = pulse.load_field(0)
    I = ionization.intensity_from_field(e)
    I_max = np.amax(I)
    
    ext = [-T/2, T/2, -X/2, X/2]
    fig = plt.figure(figsize=(4, 2), dpi=300)
    if smooth == True:
        interp = 'Spline16'
    else:
        interp = None
    if log:
        norm = colors.LogNorm(vmin=I_max*1e-4, vmax=I_max)
        im = plt.imshow(np.fliplr(np.flipud(np.transpose(I))), aspect='auto', extent=ext, cmap='viridis', norm=norm, interpolation=interp)
    else:
        im = plt.imshow(np.fliplr(np.flipud(np.transpose(I))), aspect='auto', extent=ext, cmap='viridis', vmin=0.0, interpolation=interp)
    cb = plt.colorbar(format="%0.2f")
    cb.set_label(r'Laser Intensity ($10^{14} W/cm^2$)')
    plt.xlabel('t (fs)')
    plt.ylabel(r'x ($\mathrm{\mu m}$)')
    if ylim is not None:
        plt.ylim(ylim)
    plt.tight_layout()
    i = 1
    
    def updatefig(*args):
        nonlocal i
        e[:, :] = pulse.load_field(i)[0]
        I = ionization.intensity_from_field(e)
        I = np.fliplr(np.flipud(np.transpose(I)))
        im.set_array(I)
        if log==False:
            im.autoscale()
        i += 1
        # If we run over, loop
        if i == Nz+1:
            i = 0
        if i % 30 == 0:
            print("Frame", i, "completed")
        return im,
    ani = animation.FuncAnimation(fig, updatefig, blit=True, frames=Nz-3)
    ani.save(pulse.path+name+'.mp4', fps=30)
    plt.clf()
    
def load_plasma_design(path, name=None):
    """ Load the required information from the plasma design notebook.
    
    Parameters
    ----------
    path : string
        Path to save the data at.
    name : string, optional
        Name of the pulse from the plasma source, defaults to 'Refracted_Beam'.
    
    Returns
    -------
    plasma : array of doubles
        Plasma density in the x-z plane.
    I : array of doubles
        Target laser intensity.
    z : array of doubles
        Z array corresponding to the I array.
    sim_start : double
        Z location the refraction simulation starts at.
    sim_length : double
        Length of the refraction sim.
    pulse : Pulse object
        The laser pulse object from the refraction simulation.
    """
    plasma = np.load(path+'plasma.npy', allow_pickle=True)
    I = np.load(path+'intensity.npy', allow_pickle=True)
    z = np.load(path+'z.npy', allow_pickle=True)
    sim_start, sim_length = np.load(path+'sim_size.npy', allow_pickle=True)
    
    if name is None:
        name = 'Refracted_Beam'
    
    pulseParams = {
        'Nt' : 2**6,
        'Nx' : 2**10,
        'Ny' : 2**10,
        'X' : 10e3,
        'Y' : 10e3,
        'T' : 90,
        'lam' : 0.8,
        'path' : path,
        'load' : True,
        'threads' : 4,
        'cyl' : True,
        'tau' : 30,
        'name' : name,
    }
    pulse = laserpulse.Pulse(pulseParams)
    return plasma, I, z, sim_start, sim_length, pulse

def extend_zI(z0, loc, z, I, sim_start, sim_length):
    """ Insert zeros on the front of the z and I arrays to reach optic B.
    
    Parameters
    ----------
    z0 : double
        Distance from optic B to either the start or center of the simulation
        domain. Start or center is defined by loc.
    loc : string ['start', 'center']
        If the distane to the optic is measured from the start or the center.
    z : array of doubles
        Original z array for the intensity.
    I : array of doubles
        Original intensity array. 
    sim_start : double
        Z location the refraction simulation starts at.
    sim_length : double
        Length of the refraction sim.
    
    Returns
    -------
    z_optic : array of doubles
        Extended z array.
    I_optic : array of doubles
        Extended I array.
    """
    if loc=='start':
        z_start = z0-sim_start
    elif loc=='center':
        z_start = z0-sim_start-sim_length/2
    dz = z[1]-z[0]
    Np = int(z_start/dz)
    z_optic = np.linspace(0, z[-1]+z_start, len(z)+Np)
    I_optic = np.append(np.zeros(Np), I)
    return z_optic, I_optic

def create_lens_A(ri, Ei, r, E, L, path, lam, X, Nx):
    """ Create the phase profile for lens A.

    Parameters
    ----------
    ri : array of doubles
        Radius array for the initial laser field.
    Ei : array of doubles or complex
        Electric field (magnitude or complex) of the incoming laser that will 
        pass through the lens, will be scaled to have the correct energy.
    r : array of doubles
        Radius array for the target laser field.
    E : array of doubles or complex
        Target electric field at the second lens.
    L : double
        The distance between optic A and optic B.
    path : string
        Path to save the data at.
    lam : double
        The laser wavelength in um.
    X : double
        Transverse domain size in um for lensA object.
    Nx : int
        Number of cells in the transverse dimension for lensA object.
    
    Returns
    -------
    rA : array of doubles
        Radius for each phase point.
    phiA : array of doubles
        Phase delay in radians.
    lensA : Phase object
        Lens object for the beam propagation code.
    multi : double
        The multiplier for the input field strength.
    """
    Ii = ionization.intensity_from_field(Ei)
    Io = ionization.intensity_from_field(E)
    Pi = 2*np.pi*simps(ri*Ii*1e-4, ri*1e-4)*100
    Po = 2*np.pi*simps(r*Io*1e-4, r*1e-4)*100
    multi = np.sqrt(Po/Pi)
    Ei *= multi
    Ii = ionization.intensity_from_field(Ei)
    
    # The first lens shapes the intensity on the second one
    rA, phiA = ray.lens_design(Ii, ri, Io, r, L)
    phiA *= 2*np.pi/lam
    
    # Create the first lens
    lensParams = {'Nx' : Nx,
                  'Ny' : Nx,
                  'X' : X,
                  'Y' : X,
                  'path' : path,
                  'name' : 'LensA',
                  'lam' : lam,
                  'load' : False}
    
    lensA = optic.Phase(lensParams)
    phi = lensA.reconstruct_from_cyl(rA, phiA, lensA.x, lensA.y)
    lensA.initialize_phase(phi)
    
    # The phase difference between neighboring cells must be less than pi
    dphi = np.sort(abs(phi[1:, int(Nx/2+1)]-phi[:-1, int(Nx/2+1)]))[-3]
    dphi = dphi*(Nx-1)/X
    print('Maximum phase gradient %0.4f rad/um' % dphi)
    return rA, phiA, lensA, multi

def propagate_to_lens_B(r0, E0, L, path, lam, lensA, tau, threads, plot=True):
    """ Create the phase profile for lens A.

    Parameters
    ----------
    r0 : array of doubles
        Radius array for the initial laser field.
    E0 : array of doubles or complex
        Electric field (magnitude or complex) of the incoming laser.
    L : double
        The distance between optic A and optic B.
    path : string
        Path to save the data at.
    lam : double
        The laser wavelength in um.
    lensA : Phase object
        Lens object for the beam propagation code.
    tau : double
        The RMS pulse length in fs.
    threads : int
        The number of threads to run the calculation on.
    plot : bool
        Show plots or not.
    
    Returns
    -------
    beam0 : LaserBeam object
        The laser beam object that propagates through the system.
    """
    # Create the initial beam to pass through the lens
    XA = lensA.X
    NxA = lensA.Nx
    beamParams = {'Nx' : NxA,
                  'Ny' : NxA,
                  'X' : XA,
                  'Y' : XA,
                  'lam' : lam,
                  'path' : path,
                  'name' : 'Beam0_A_to_B',
                  'threads' : threads,
                  'cyl' : True,
                  'load' : False}
    
    # Super Gaussian for simulation
    beam = laserbeam.Laser(beamParams)
    e = beam.reconstruct_from_cyl(r0, E0, beam.x, beam.y)
    beam.initialize_field(e)
    I0 = ionization.intensity_from_field(E0)
    print('Total input energy %0.2fmJ' % (beam.total_cyl_power(r0, I0)*tau))
    
    # Plot the initial transverse intensity
    #--------------------------------------------------------------------------
    if plot:
        e = beam.e
        I = beam.intensity_from_field(e)
        I = beam.prep_data(I)*1e4
        plt.figure(figsize=(8, 3), dpi=150)
        plt.subplot(121)
        plt.imshow(I, aspect='auto', extent=[-XA/2e3, XA/2e3, -XA/2e3, XA/2e3])
        cb = plt.colorbar()
        cb.set_label(r'Intensity ($10^{10} W/cm^2$)')
        plt.set_cmap('viridis')
        plt.xlabel(r'x (mm)')
        plt.ylabel(r'y (mm)')
        plt.title('Peak transverse intensity at Optic A')
    
    # Propagate to the second lens and plot the intensity
    #--------------------------------------------------------------------------
    interactions.beam_phase(beam, lensA)
    beam.propagate(L, 1.0)
    if plot:
        e = beam.e
        I = beam.intensity_from_field(e)
        I = beam.prep_data(I)*1e4
        plt.subplot(122)
        plt.imshow(I, aspect='auto', extent=[-XA/2e3, XA/2e3, -XA/2e3, XA/2e3])
        cb = plt.colorbar()
        cb.set_label(r'Intensity ($10^{10} W/cm^2$)')
        plt.set_cmap('viridis')
        plt.xlabel(r'x (mm)')
        plt.ylabel(r'y (mm)')
        plt.title('Peak transverse intensity at Optic B')
        plt.tight_layout()
        plt.show()
    
    return beam

def create_lens_B(beam0, r, E, path, lam, X, Nx):
    """ Create the phase profile for lens B.

    Parameters
    ----------
    beam0 : LaserBeam object
        The laser beam object that propagated through the system.
    r : array of doubles
        Radius array for the target laser field.
    E : array of doubles or complex
        Target electric field at the second lens.
    path : string
        Path to save the data at.
    lam : double
        The laser wavelength in um.
    X : double
        Transverse domain size in um for lensA object.
    Nx : int
        Number of cells in the transverse dimension for lensA object.
    
    Returns
    -------
    rB : array of doubles
        Radius for each phase point.
    phiB : array of doubles
        Phase delay in radians.
    lensB : Phase object
        Lens object for the beam propagation code.
    """
    # The second lens (B) removes the phase from the first lens and adds an axicon like phase
    r0 = -beam0.x[:int(beam0.Nx/2+1)]
    e0 = beam0.e[:int(beam0.Nx/2+1), int(beam0.Ny/2+1)]
    phi0 = np.unwrap(np.angle(e0))
    phi0 = phi0 - phi0[-1]
    phiB = np.unwrap(np.angle(E)) - beam0.reconstruct_from_cyl(r0, phi0, r, np.zeros(1))[:, 0]
    
    # Create the second lens, the domain is smaller for this one
    lensParams = {'Nx' : Nx,
                  'Ny' : Nx,
                  'X' : X,
                  'Y' : X,
                  'path' : path,
                  'name' : 'LensB',
                  'lam' : lam,
                  'load' : False}
    
    lensB = optic.Phase(lensParams)
    phi = lensB.reconstruct_from_cyl(r, phiB, lensB.x, lensB.y)
    lensB.initialize_phase(phi)
    dphi = np.sort(abs(phi[1:, int(Nx/2+1)]-phi[:-1, int(Nx/2+1)]))[-3]
    dphi = dphi*(Nx-1)/X
    print('Maximum phase gradient %0.2f rad/um' % dphi)
    
    return r, phiB, lensB

def plot_phase(rA, phiA, rB, phiB, xlimB):
    """ Plot the phases of the two lenses creating the beam.
    
    Parameters
    ----------
    rA : array of doubles
        Radius array for lens A.
    phiA : array of doubles
        Phase of lens A in radians.
    rB : array of doubles
        Radius array for lensB
    phiB : array of doubles
        Phase of lens B in radians.
    xlimB : array or tuple
        X limits of the optic B plot.
    """
    plt.figure(figsize=(8, 2), dpi=150)
    plt.subplot(121)
    plt.plot(rA/1e3, phiA)
    plt.xlabel(r'r (mm)')
    plt.ylabel(r'$\phi$ (radians)')
    plt.xlim(0, np.amax(rA)/1e3)
    plt.subplot(122)
    plt.plot(rB/1e3, phiB)
    plt.xlabel(r'r (mm)')
    plt.ylabel(r'$\phi$ (radians)')
    plt.xlim(xlimB)
    plt.show()
    
def field_after_lens_B(beam0, rB, phiB, r, E, rlim=None, plot=True):
    """ Calculate the radial field after the beam shaping optics.
    
    Parameters
    ----------
    beam0 : LaserBeam object
        The laser beam object that propagates through the system.
    rB : array of doubles
        Radius for each phase point.
    phiB : array of doubles
        Phase delay in radians.
    r : array of doubles
        Radius array for the target laser field.
    E : array of doubles or complex
        Target electric field at the second lens.
    rlim : optional, tuple or array
        The limits for the transverse field plot.
    plot : bool
        Show plots or not.
    
    Returns
    -------
    r1 : array of doubles
        Radius vector for the electric field.
    e1 : array of complex
        Electric field after the beam shaping optics.
    """
    r0 = -beam0.x[:int(beam0.Nx/2+1)]
    e0 = beam0.e[:int(beam0.Nx/2+1), int(beam0.Ny/2+1)]
    e = interp1d(r0, e0)
    phi = interp1d(rB, phiB)
    e1 = e(r)*np.exp(1j*phi(r))
    
    if rlim is None:
        rlim = [0, np.amax(r)/1e3]
    if plot:
        plt.figure(figsize=(8, 2), dpi=150)
        plt.subplot(121)
        plt.plot(r/1e3, ionization.intensity_from_field(E)/1e-4)
        plt.plot(r/1e3, ionization.intensity_from_field(e1)/1e-4, 'm--')
        plt.xlabel(r'r (mm)')
        plt.ylabel(r'I ($\mathrm{10^{10}W/cm^2}$)')
        plt.xlim(rlim)
        plt.subplot(122)
        plt.plot(r/1e3, np.unwrap(np.angle(E)))
        plt.plot(r/1e3, np.unwrap(np.angle(e1))-np.unwrap(np.angle(e1))[0], 'm--')
        plt.xlabel(r'r (mm)')
        plt.ylabel(r'$\phi$ (radians)')
        plt.xlim(rlim)
        plt.show()
    return r, e1

def propagate_down_beampipe(pulse, apertures, diameters, X, Nx, Nz, ret_pulse=None):
    """ Propagate the pulse through the apertures in the beamline.
    
    Parameters
    ----------
    pulse : LaserPulse object
        The laser pulse used in calculating the refraction.
    apertures : array of doubles
        Distance from the end of the plasma to each aperture in the system.    
    diameters : array of doubles
        The diameter of each aperture in the beam pipe.
    X : double
        Transverse domain size in um for the pulse.
    Nx : int
        Number of cells in the transverse dimension for the pulse.
    Nz : int
        Total number of cells to use down the beampipe
    ret_pulse : bool
        Whether to return the pulse or not.
    
    Returns
    -------
    beam_hm : LaserPulse object
        The pulse that propagates through the beamline.
    """
    pulseParams = pulse.params.copy()
    pulseParams['name'] = 'Beamline_Pulse'
    pulseParams['Nx'] = Nx
    pulseParams['Ny'] = Nx
    pulseParams['X'] = X
    pulseParams['Y'] = X
    pulseParams['load'] = False
    pulse_bl = laserpulse.Pulse(pulseParams)
    e0 = np.array(pulse.e)
    e = pulse_bl.e
    for i in range(pulse.Nt):
        e[i, :, :] = pulse_bl.reconstruct_from_cyl(pulse.x, e0[i, :, int(pulse.Ny/2)], pulse_bl.x, pulse_bl.y)
    pulse_bl.initialize_field(e)
    
    holeParams = {
        'Nx' : Nx,
        'Ny' : Nx,
        'X' : X,
        'Y' : X,
        'path' : pulseParams['path'], 
        'lam' : pulseParams['lam'],
        'load' : False
    }
    
    Z = apertures[-1]
    dz = Z/(Nz-1)
    N = len(apertures)
    
    prev = 0
    m_tot = 0
    for i in range(N):
        l = apertures[i]-prev
        prev = apertures[i]
        m = int(l/dz)
        m_tot += m
        z = np.linspace(0, l, m)
        pulse_bl.propagate(z, 1.0)
        holeParams['name'] = 'Aperture_'+str(i)
        holeParams['r'] = diameters[i]/2
        hole = optic.Aperture(holeParams)
        interactions.beam_intensity(pulse_bl, hole)
    I_bl = np.zeros((Nz, Nx), dtype='double')
    for i in range(m_tot):
        I_bl[i, :] = np.amax(pulse_bl.intensity_from_field(pulse_bl.load_field(i+1)[0]), axis=0)
    I_bl = pulse_bl.prep_data(I_bl)
    if ret_pulse is None:
        return I_bl
    else:
        return I_bl, pulse_bl
    
    
        
