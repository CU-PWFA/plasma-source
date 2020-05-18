#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing functions for propagating two THz pulses through a crystal.
Created on Fri May 31 10:20:09 2019

@author: keenan
"""
import copy;
from cycler import cycler;
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
import scipy.fftpack as sfft;

# Crystal properties module
import eocry as ec;

# Constants
eps0 = 8.854187817e-12;

# Colors for plotting.
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77','#CC6677', '#882255',  '#AA4499'];
cy = cycler('color', plot_colors);


def thz_field_raw(params):
    '''
    Function to calculate the raw THz electric field of an electron 
    bunch or bunches and its FFT.
    
    Parameters:
    -----------
    params : dictionary
             The params dictionary should contain
                t      : array_like
                         Time array in s.
                beam1   : dictionary
                          A dictionary containing 'Q':charge (C), 
                          'x' : horizontal offset from crystal (m), 
                          'y' : vertical 'sig_t': temporal rms length (s)
                          for the primary bunch.
                beam2  : dictionary
                         A dictionary containing the same parameters as beam1, 
                         with values for the trailing beam.
                t_off  : float
                         The longitudinal temporal offset of the beams in s.
                plot   : bool, optional
                         Whether or not to plot the field, its ft and its ift
    
    Returns:
    --------
    Er1  : array_like
           The radial electric field of the primary beam
    Er2  : array_lole
           The radial electric field of the trailing beam
    Er   : array_like
           The radial electric field of the combined bunches
    t    : array_like
           The time domain of the field
    FEr1 : array_like
          The FT of the THz field from beam 1.
    FEr2 : array_like
           The FT of the THz field from beam2. 
    FEr  : The FT of the combined THz field
    f    : array_like
           The frequency domain of the FT
    '''
    
    # Extract variables from params
    t     = params['t'];
    beam1 = params['beam1'];
    beam2 = params['beam2'];
    t_off = params['t_off'];
    #print("field raw t_off", t_off);
    # Optional parameters
    if bool(params.get('plot')):
        plot = params['plot'];
    else:
        plot = False;
    
    # Compute Electric field amplitude for both beams
    r1 = np.sqrt(beam1['x'] ** 2 + beam1['y']**2)
    E01 = beam1['Q'] / (2 * np.pi * eps0 * r1 * np.sqrt(2 * np.pi) \
         * c * beam1['sig_t']);
    Er1 = E01 * np.exp(-t**2 / (2 * beam1['sig_t']**2));
    r2  = np.sqrt(beam2['x']**2 + beam2['y']**2)
    E02 = beam2['Q'] / (2 * np.pi * eps0 * r2 * np.sqrt(2 * np.pi) \
         * c * beam2['sig_t']);
    
    t2  = t + t_off;
    Er2 = E02 * np.exp(-t2**2 / (2 * beam2['sig_t']**2));
    Er  = Er1 + Er2;
    
    # Calculate FFTs
    Ts = abs(t[1] - t[0]) # sampling time interval
    Fs = 1/Ts; # Sampling frequency
    L  = len(t); # Number of samples
    # FFT is significantly faster if the array length is a power of 2, pad array
    # with zeros to achieve this
    def nextpow2(i):
        n = 1;
        while 2**n < i : n += 1
        return n
    NFFT   = 2**(nextpow2(L) + 1)
    f      = Fs * np.linspace(0, 1, NFFT)*1e-12; # THz frequencies
    f_plot = (Fs/2) * np.linspace(0, 1, int(NFFT / 2 + 1)) * 1e-12; # Nyquist
    FEr    = sfft.fft(Er, NFFT) / L;
    FEr1   = sfft.fft(Er1, NFFT) / L;
    FEr2   = sfft.fft(Er2, NFFT) / L;
    # As a general check calculate the IFFT
    Nsamp = 2 * (len(f_plot) - 1);
    fmax  = f_plot[-1]*1e12; # Nyquist freqency
    
    tt    = np.arange(0, Nsamp - 1) * (1 / (2 * fmax)); # s, time axis
    IFER  = sfft.ifft(FEr) * L;
    IFER1 = sfft.ifft(FEr1) * L;
    IFER2 = sfft.ifft(FEr2) * L;
    # Remove extra zeros
    
    tdiff = len(tt) - len(t);
    tt    = tt[0:-1 - tdiff + 1];
    IFER  = IFER[0:-1 - tdiff];
    IFER1 = IFER1[0:-1 - tdiff];
    IFER2 = IFER2[0:-1 - tdiff];
    # recenter time axis on '0', only works if input t is symmetric
    tt    = tt - tt[-1]/2;
    tt    = tt * 1e15;
    
    if plot:
        
        xlims = np.array([-5 * beam2['sig_t'], 5 * beam1['sig_t']]) * 1e15;
        fig1 = plt.figure(figsize = (4,4), dpi = 200);
        ax1  = fig1.gca();
        ax1.set_prop_cycle(cy)
        ax1.plot(t * 1e15, Er, '-', label = 'Combined');
        ax1.plot(t * 1e15, Er1, '--', label = 'Drive');
        ax1.plot(t * 1e15, Er2, '-.', label = 'Witness');
        plt.legend()
        ax1.set_xlabel('t [fs]')
        ax1.set_ylabel(r'$E_r$(t) [V/m]');
        #ax1.set_xlim(xlims)
        ax1.set_title('Input pulse');
       
        fig2 = plt.figure(figsize = (4,4), dpi = 200);
        ax2  = fig2.gca();
        ax2.set_prop_cycle(cy)
        ax2.plot(f_plot, 2 * abs(FEr[1:int(NFFT / 2 + 2)]), '-', \
                 label = 'Combined')
        ax2.plot(f_plot, 2 * abs(FEr1[1:int(NFFT / 2 + 2)]), '-', \
                 label = 'Drive');
        ax2.plot(f_plot, 2 * abs(FEr2[1:int(NFFT / 2 + 2)]), '-.', \
                 label = 'Witness');
        plt.legend();
        ax2.set_xlim([0, 10]);
        ax2.set_xlabel('f [THz]');
        ax2.set_ylabel(r'FT[$E_r$](f)')
        ax2.set_title('FFT');
        
        fig3 = plt.figure(figsize = (4,4), dpi = 200);
        ax3  = fig3.gca();
        ax3.set_prop_cycle(cy)
        ax3.plot(tt, IFER, '-', label = 'Combined');
        ax3.plot(tt, IFER1, '--', label = 'Drive');
        ax3.plot(tt, IFER2, '-.', label = 'Witness');
        plt.legend();
        peak_ind = np.argmax(IFER);
        t0       = tt[peak_ind]
        xlims = [t0 - 5 * beam2['sig_t'] * 1e15, \
                 t0 + 5 * beam1['sig_t'] * 1e15];
        #ax3.set_xlim(xlims);
        ax3.set_xlabel('t [fs]');
        ax3.set_ylabel(r'$E_r$(t) [V/m]');
        ax3.set_title('IFFT');
        
        plt.show();
    return Er1, Er2, Er, t, FEr1, FEr2, FEr, f;

def thz_field_cry(params):
    '''
    Function to calcualte the THz electric field in the crystal as a function 
    of time and its FFT accounting for transmission and the eo response. Also 
    plots a delta-like probe pulse propagating through the crystal if wanted. 
    
    Parameters:
    -----------
    params : dictionary
             The params dictionary should contain
                ctype   : str
                          Specifies the EO crystal, either GaP or ZnTe
                t       : array_like
                          Time array in s.
                d       : float
                          Crystal thickness in m.
                nslice  : int
                          The number of slices to calculate the E-field for.
                beam1   : dictionary
                          A dictionary containing 'Q':charge (C), 
                          'x' : horizontal offset from crystal (m), 
                          'y' : vertical 'sig_t': temporal rms length (s)
                          for the primary bunch.
                beam2   : dictionary
                          A dictionary containing the same keys as beam1, with 
                          values for the trailing bunch.
                t_off   : float
                          The longitudinal temporal offset of the bunches. 
                plot    : bool, optional
                          Whether or not to plot the field.
                verbose : bool, optional
                          Whether or not to include print statements. 
                pockel  : bool, otpional
                          Whether or not to account for the Pockel effect when
                          propagating the THz pulse. default = True.
                ref     : bool, optional
                          Whether or not to include reflection of the THz pulse
                          default = False.
                probe   : dictionary, optional
                          A dictionary containing y0: the central wavelenth, 
                          z_match : the depth at which the probe is aligned with
                          the peak of the THz pulse, a_laser : the angle the 
                          laser makes with the EO crystal in radians. 

    Returns:
    --------
    Ec    : array_like
            Re[Ec] is the electric field in the crystal at depth d
    tt    : array_like
            The time axis of Ec, in ps
    FEc   : array_like
            The Fourier components of the THz field in the crystal at depth d
    sf    : array_like
            The frequency axis of FEc in THz
    d_arr : array_like
            The depths into the crystal at which Ec is calculated, in m
    '''
    # Extract variables from params
    ctype  = params['ctype'];
    d      = params['d'];
    t      = params['t'];
    nslice = params['nslice']
    #print("field cry t_off", params['t_off'])
    # Optional parameters
    if bool(params.get('plot')):
        plot = params['plot'];
    else:
        plot = False;
    if bool(params.get('verbose')):
        verbose = params['verbose'];
    else:
        verbose = False;
    if bool(params.get('pockel')):
        pockel = params['pockel'];
    else:
        pockel = True;
    if bool(params.get('ref')):
        ref = params['ref'];
    else:
        ref = False;
    if bool(params.get('probe')):
        probe = params['probe'];
        plot_probe = True;
    else:
        plot_probe = False
    
    # Calculate the raw fields and their FFT first (don't plot)
    raw_params = copy.deepcopy(params); raw_params['plot'] = False;
    Er1, Er2, Er, t, FEr1, FEr2, FEr, f = thz_field_raw(raw_params);
    
    # Only take first half (forward propagation)
    sFEr1  = FEr1[0:int(len(FEr1) / 2 + 1)];
    sFEr2  = FEr2[0:int(len(FEr2) / 2 + 1)];
    sFEr   = FEr[0:int(len(FEr) / 2 + 1)];
    sf     = f[0:int(len(f) / 2 + 1)];
    f_Hz   = sf * 1e12;
    
    # Calculate and print average phase velocity of frequency range
    
    if verbose:
        v_ph, v_g, dum = ec.velocities({'ctype' : ctype, 'f' : f_Hz});
        v_ph_c         = np.nanmean(v_ph)/c; 
        v_g_c          = np.nanmean(v_g) / c;
        print(f'THz average phase velocity = {v_ph_c}c');
        print(f'THz average group velocity = {v_g_c}c')

    # Compute transmission, attenuation, & pockel coeffs and ind. of ref. \
    eo_params     = {'ctype' : ctype, 'f' : f_Hz}
    A             = ec.transcoeff(eo_params);
    eps, n, kappa = ec.dielec(eo_params);
    if pockel:
        r41           = ec.eocoeff(eo_params);
    else:
        r41           = 1;
    
    # Fourier parameters
    L = len(t);
    
    # Create array of depths
    j      = np.arange(1, nslice, 1);
    dz     = d / nslice;
    d_arr  = (j - 0.5) * dz;
    # Preallocate for loop
    FEc_drive = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    FEc_wit   = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    FEc       = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    
    Ec_drive  = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    Ec_wit    = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    Ec        = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    
    # Propagate through the crystal
    
    for i in range(len(d_arr)):
        FEc_drive[:,i] = A * sFEr1 * r41 * np.exp((1j * 2 * np.pi * f_Hz / c) \
                   * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                   * d_arr[i]);
                 
        FEc_wit[:, i]  = A * sFEr2 * r41 * np.exp((1j * 2 * np.pi * f_Hz / c) \
                   * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                   * d_arr[i]);
               
        FEc[:,i] = A * sFEr * r41 * np.exp((1j * 2 * np.pi * f_Hz / c) \
                   * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                   * d_arr[i]);
           
        Ec_drive[:,i]  = sfft.ifft(FEc_drive[:,i]) * L;
        Ec_drive[:, i] = np.flip(Ec_drive[:, i], axis = 0)
           
        Ec_wit[:,i]    = sfft.ifft(FEc_wit[:,i]) * L;
        Ec_wit[:, i]   = np.flip(Ec_wit[:, i], axis = 0)
        
        Ec[:,i]        = sfft.ifft(FEc[:,i]) * L;
        Ec[:, i]       = np.flip(Ec[:, i], axis = 0)
        #Ec[:,i]  = sfft.ifftshift(Ec[:,i]);    
    # Add reflection if wanted:
    if ref:
        A_ref     = (1 - n - 1j * kappa) / (1 + n + 1j * kappa);
        k         = 2 * np.pi * f_Hz * n / c; 
        alpha     = 2 * np.pi * f_Hz * kappa / c;
        ref_phase = np.exp(1j * 2 * k * d) * np.exp(-2 * alpha * d);
        
        F2_drive = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
        F2_wit   = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
        F2       = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
        
        Eref_drive = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
        Eref_wit   = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
        Eref       = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
        for i in range(len(d_arr)):
            F2_drive[:, i] = A * sFEr1 * r41 * A_ref**2 * ref_phase  \
                               * np.exp((1j * 2 * np.pi * f_Hz / c) \
                               * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                               * d_arr[i]);
                                        
            F2_wit[:, i]   = A * sFEr2 * r41 * A_ref**2 * ref_phase  \
                               * np.exp((1j * 2 * np.pi * f_Hz / c) \
                               * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                               * d_arr[i]);
                                        
            F2[:, i]       = A * sFEr * r41 * A_ref**2 * ref_phase  \
                               * np.exp((1j * 2 * np.pi * f_Hz / c) \
                               * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                               * d_arr[i]);
                                        
            Eref_drive[:, i] = sfft.ifft(F2_drive[:, i]) * L;
            Eref_drive[:, i] = np.flip(Eref_drive[:, i], axis = 0);
            Eref_wit[:, i]   = sfft.ifft(F2_wit[:, i]) * L;
            Eref_wit[:, i]   = np.flip(Eref_wit[:, i], axis = 0)                            
            Eref[:, i]       = sfft.ifft(F2[:, i]) * L;
            Eref[:, i]       = np.flip(Eref[:, i], axis = 0);
    
    # Reconstruct time array:
    Nsamp = len(sf);
    fmax  = sf[-1] * 1e12; 
    tt    = np.arange(0, Nsamp - 1) * (1 / (fmax));
    
    tdiff = len(tt) - len(t);
    
    Ec_drive = Ec_drive[0:int(len(tt) - tdiff + 1), :];
    Ec_wit   = Ec_wit[0:int(len(tt) - tdiff + 1), :];
    Ec       = Ec[0:int(len(tt) - tdiff + 1), :];
    
    if ref:
        Eref_drive = Eref_drive[0:int(len(tt) - tdiff + 1), :]
        Ec_drive   = Ec_drive + Eref_drive;
        
        Eref_wit = Eref_wit[0:int(len(tt) - tdiff + 1), :]
        Ec_wit   = Ec_wit + Eref_wit;
        
        Eref = Eref[0:int(len(tt) - tdiff + 1), :]
        Ec = Ec + Eref;
    tt    = tt[0:int(len(tt) - tdiff + 1)];
    tt    = tt * 1e12; # ps 
    # Shift time axis accordingly
    
    tt    = tt - 3 * tt[-1] / 4;
    
    
    # Values for plotting probe pulse if requested. 
    if plot_probe:
        f_opt       = np.array([c / probe['y0']]);
        dum, vg_opt, dum = \
                          ec.velocities({'ctype' : ctype, 'f' : f_opt});
        vg_eff = vg_opt * np.cos(probe['a_laser']);
        # Add shift so probe is aligned to z_match
        if probe['z_match'] != 0:
            match_ind = np.argmin(abs(d_arr - probe['z_match']));
            
        else:
            match_ind = 0;
        t_peak  = tt[np.argmax(np.real(Ec[:, match_ind]))];
        t_shift = t_peak - (d_arr[match_ind] / vg_eff * 1e12);
        probe['t_shift'] = t_shift;
        probe['vg_eff']  = vg_eff;
        if verbose:
            print("Effective probe velocity =", \
                   str(vg_eff / c), "c"); 
    if plot:
        height = 20;
        fig1 = plt.figure(figsize = (6, height), dpi = 200);
        
        # Get indices to plot (5 points in the crystal);
        z_step = int(d *1e6 / 5); # microns
        z_start = (z_step / 2);
        z_end   = (z_step * (5 - .5)); 
        z_plot = np.arange(z_start, z_end + 1, z_step) * 1e-6; # m
        
        
        # Plot, bottom to top
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[4]));
        ax1           = fig1.add_subplot(515);
        ax1.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind]);
        Ec_plot       = np.real(Ec[:, ind]);
        print(max(Ec_plot));
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax1.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');
        ax1.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax1.plot(tt, Ec_wit_plot, '-.', label = 'Witness');
        ax1.set_xlabel('t [ps]');
        # Peak location
        if verbose:
            print(np.round(d_arr[ind] * 1e6, 2), \
                  np.round(tt[np.argmax(Ec_plot)], 2), \
                  np.round(tt[np.argmax(Ec_drive_plot)], 2), \
                  np.round(tt[np.argmax(Ec_wit_plot)], 2));
        if plot_probe:
            t_probe = (z_plot[4] / vg_eff) * 1e12
            ax1.axvline(x = t_probe + t_shift, color = 'r');
            
        # Pretty plot stuff
        ax1.set_yticks([]);
        ax1.text(0.05, 0.85, depth, ha='center', va='center', \
                 transform=ax1.transAxes)
        ind2 = np.argmax(Ec_plot);
        
        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[3]));
        ax2 = fig1.add_subplot(514, sharex = ax1);
        ax2.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind]);
        Ec_plot       = np.real(Ec[:, ind]);
        print(max(Ec_plot))
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax2.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');
        ax2.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax2.plot(tt, Ec_wit_plot, '-.', label = 'Witness');
        # Peak location
        if verbose:
            print(np.round(d_arr[ind] * 1e6, 2), \
                  np.round(tt[np.argmax(Ec_plot)], 2), \
                  np.round(tt[np.argmax(Ec_drive_plot)], 2), \
                  np.round(tt[np.argmax(Ec_wit_plot)], 2));
        if plot_probe:
            t_probe = (z_plot[3] / vg_eff) * 1e12
            ax2.axvline(x = t_probe + t_shift, color = 'r');
        # Pretty plot stuff
        ax2.text(0.05, 0.85, depth, ha='center', va='center', \
                 transform=ax2.transAxes)
        ax2.set_yticks([]);
        
        
        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[2]));
        ax3 = fig1.add_subplot(513, sharex = ax1);
        ax3.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind]);
        Ec_plot       = np.real(Ec[:, ind]);
        print(max(Ec_plot))
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax3.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');
        ax3.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax3.plot(tt , Ec_wit_plot, '-.', label = 'Witness');
        # Peak location
        if verbose:
            print(np.round(d_arr[ind] * 1e6, 2), \
                  np.round(tt[np.argmax(Ec_plot)], 2), \
                  np.round(tt[np.argmax(Ec_drive_plot)], 2), \
                  np.round(tt[np.argmax(Ec_wit_plot)], 2));
        if plot_probe:
            t_probe = (z_plot[2] / vg_eff) * 1e12
            ax3.axvline(x = t_probe + t_shift, color = 'r');
        # Pretty plot stuff
        ax3.text(0.05, 0.85, depth, ha='center', va='center', \
                 transform=ax3.transAxes)
        ax3.set_yticks([]);

        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[1]));
        ax4 = fig1.add_subplot(512, sharex = ax1);
        ax4.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind]);
        Ec_plot       = np.real(Ec[:, ind]);
        print(max(Ec_plot))

        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax4.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');
        ax4.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax4.plot(tt, Ec_wit_plot, '-.', label = 'Witness');
        # Peak location
        if verbose:
            print(np.round(d_arr[ind] * 1e6, 2), \
                  np.round(tt[np.argmax(Ec_plot)], 2), \
                  np.round(tt[np.argmax(Ec_drive_plot)], 2), \
                  np.round(tt[np.argmax(Ec_wit_plot)], 2));
        if plot_probe:
            t_probe = (z_plot[1] / vg_eff) * 1e12
            ax4.axvline(x = t_probe + t_shift, color = 'r');
        # Pretty plot stuff
        ax4.text(0.05, 0.85, depth, ha='center', va='center', \
                 transform=ax4.transAxes)
        ax4.set_yticks([]);

        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[0]));
        ax5 = fig1.add_subplot(511, sharex = ax1);
        ax5.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind]);
        Ec_plot       = np.real(Ec[:, ind]);
        print(max(Ec_plot))

        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax5.plot(tt, Ec_plot, linestyle = '-',\
                 label = r'$E_{eff}$');
        ax5.plot(tt, Ec_drive_plot, linestyle = '--', \
                 label = 'Drive');
        ax5.plot(tt, Ec_wit_plot, linestyle = '-.', \
                 label = 'Witness');
        
        ax5.legend(ncol = 1, loc = 'upper left', bbox_to_anchor = (-0.15, 1.0), \
                   frameon = False);
        
        
        # Peak location
        if verbose:
            print(np.round(d_arr[ind] * 1e6, 2), \
                  np.round(tt[np.argmax(Ec_plot)], 2), \
                  np.round(tt[np.argmax(Ec_drive_plot)], 2), \
                  np.round(tt[np.argmax(Ec_wit_plot)], 2));
        if plot_probe:
            t_probe = (z_plot[0] / vg_eff) * 1e12
            ax5.axvline(x = t_probe + t_shift, color = 'r');
            
        # Pretty plot stuff
        ax5.set_yticks([]);
        ind = np.argmax(Ec_wit_plot)
        ax5.text(0.05, 0.85, depth, ha='center', va='center', \
                 transform=ax5.transAxes)

        ########################################################################
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        if tt[ind2] + 2 > tt[-1]:
            xlims = [0, tt[-1]]
        else:
            xlims = [0, tt[ind2] + 1.0];
        ax1.set_xlim(xlims);
        plt.subplots_adjust(hspace = 0.0);
    return Ec_drive, Ec_wit, Ec, tt, FEc_drive, FEc_wit, FEc, sf, d_arr;
