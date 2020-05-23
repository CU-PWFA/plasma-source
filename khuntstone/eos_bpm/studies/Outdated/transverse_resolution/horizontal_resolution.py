#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for studying the detector response at different horizontal beam offsets
Created on Thu Jun  6 12:47:24 2019

@author: keenan
"""

import cycler
import matplotlib.pyplot as plt;
import numpy as np;
from scipy import optimize;
import scipy.interpolate as interp;
import time;


import phase_retard as pr;

# Colors for plotting.
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77','#CC6677', '#882255',  '#AA4499'];

cy = cycler.cycler('color', plot_colors);


def offset_signal(params, x, x_off, sname, plot = False):
    '''
    Function to compute the phase retardation parameter of the drive and 
    witness beam for a range of horizontal offsets from the crystal.
    
    Parameters:
    -----------
    params   : dictionary, see phase_retard.phase_retard() parameters
    x        : float, 
               crystals horrizontal offset from beamline (m)
    x_off    : array_like
               Beam offsets from the beamline center (m)
    sname    : str
               Since this function takes a while it is written to save the 
               signals to a .np file, sname is the filename. 
    plot     : bool, optional
               Whether or not to plot peak signal vs. offset
    '''
    # Preallocate for loop
    drive_sigs = np.zeros((len(x_off), len(params['tau'])));
    wit_sigs   = np.zeros((len(x_off), len(params['tau'])));
    
    beam1      = params['beam1'];
    beam2      = params['beam2'];
    
    start = time.time();
    for i in range(len(x_off)):
        beam1['x'] = x + x_off[i];
        beam2['x'] = x + x_off[i];
        params['beam1'] = beam1;
        params['beam2'] = beam2;
        drive_sigs[i, :], wit_sigs[i, :], dum, t_sig = \
                                                      pr.phase_retard(params);
    print(time.time() - start)
    if plot:
        fig1  = plt.figure(figsize = (6, 6), dpi = 200);
        fig2  = plt.figure(figsize = (6, 6), dpi = 200);
        ax1   = fig1.gca(); ax1.set_prop_cycle(cy);
        ax2   = fig2.gca(); ax2.set_prop_cycle(cy);
        bal_drive = np.sin(drive_sigs);
        bal_drive = np.max(bal_drive, axis = 1);
        x_drive   = np.sin(drive_sigs / 2)**2;
        x_drive   = np.max(x_drive, axis = 1);
        
        bal_wit = np.sin(wit_sigs);
        bal_wit = np.max(bal_wit, axis = 1);
        x_wit   = np.sin(wit_sigs / 2)**2;
        x_wit   = np.max(x_wit, axis = 1);
        # Balanced setup
        ax1.plot(x_off * 1e6, bal_drive, label = 'Drive');
        ax1.plot(x_off * 1e6, bal_wit, label = 'Witness');
        ax1.set_ylabel(r'Peak EO signal')
        ax1.set_xlabel(r'$\tau$ [fs]');
        ax1.set_title("Balanced detector")
        ax1.legend();
        # Crossed polarizer
        ax2.plot(x_off * 1e6, x_drive, label = 'Drive');
        ax2.plot(x_off * 1e6, x_wit, label = 'Witness');
        ax2.set_ylabel(r'Peak EO signal')
        ax2.set_xlabel(r'$\tau$ [fs]');
        ax2.set_title("Crossed polarizer")
        ax2.legend();
    fpath = "/home/keenan/eos_bpm/python/outputs/horizontal_signal/"
    np.save(fpath + sname, (drive_sigs, wit_sigs, t_sig));
def plot_sigs(drive_sigs, wit_sigs, x_off, t_sigs, x_plot, r, det):
    '''
    Function to plot the drive beam and witness beam signals for different 
    detector setups. 
    Parameters:
    -----------
    drive_sigs : array_like
                2D array of drive phase retardation axis 1 = x offset, 
                axis 2 = probe delay
    wit_sigs   : array_like
                 Similar array for witness signal
    x_off      : array_like
                 Array of the horizontal offsets corresponding to drive_sigs 
                 and wit sigs.
    t_sigs     : array_like
                 Time axis for plotting signals
    x_plot     : array_like
                 Horizontal offsets at which to plot the signal  
    r          : float
                 The crystal to beamline distance in mm (for labeling)
    det        : str
                 The detection scheme either bal or cross
    '''
    
    # First create interpolation functions 
    f_drive = interp.interp1d(x_off, drive_sigs, axis = 0);
    f_wit   = interp.interp1d(x_off, wit_sigs, axis = 0);
    
    # Plot drive, witness, and combined
    fig1 = plt.figure(figsize = (8, 6), dpi = 200);
    fig2 = plt.figure(figsize = (8, 6), dpi = 200);
    fig3 = plt.figure(figsize = (8,6), dpi  = 200);
    
    # Axes for drive, witness, combined
    ax1  = fig1.gca(); 
    ax2  = fig2.gca();
    ax3  = fig3.gca();
    
    # Line colors
    ax1.set_prop_cycle(cy);
    ax2.set_prop_cycle(cy);
    ax3.set_prop_cycle(cy);
    
    fs = 16;
    ts = 36;
    # x label time
    ax1.set_xlabel(r'$\tau$ [fs]', fontsize = fs);
    ax2.set_xlabel(r'$\tau$ [fs]', fontsize = fs);
    ax3.set_xlabel(r'$\tau$ [fs]', fontsize = fs);
    
    # y label and title (varies by setup)
    if det == 'bal':
        ax1.set_ylabel(r'sin$\Gamma$', fontsize = fs);
        ax2.set_ylabel(r'sin$\Gamma$', fontsize = fs);
        ax3.set_ylabel(r'sin$\Gamma$', fontsize = fs);
        
        ax1.set_title("Balanced detector drive signal " \
                      + r'($r_0$ = ' + str(r) + 'mm)', fontsize = ts, \
                      fontweight = 'bold');
        ax2.set_title("Balanced detector witness signal " \
                      + r'($r_0$ = ' + str(r) + 'mm)', fontsize = ts, \
                      fontweight = 'bold');
        ax3.set_title("Balanced detector signal " \
                      + r'($r_0$ = ' + str(r) + 'mm)', fontsize = ts, \
                      fontweight = 'bold');
    elif det == 'cross':
        ax1.set_ylabel(r'sin$^2 \Gamma / 2$', fontsize = fs);
        ax2.set_ylabel(r'sin$^2 \Gamma / 2$', fontsize = fs);
        ax3.set_ylabel(r'sin$^2 \Gamma / 2$', fontsize = fs);
        
        ax1.set_title(r'Drive signal',\
                      fontsize = ts, \
                      fontweight = 'bold');
        ax2.set_title(r'Witness signal', \
                      fontsize = ts, \
                      fontweight = 'bold');
        ax3.set_title(r'sin$^2(\Gamma/2)$', \
                      fontsize = ts, \
                      fontweight = 'bold');
    else:
        ax1.set_ylabel(r'$\Gamma$', fontsize = fs)
        ax2.set_ylabel(r'$\Gamma$', fontsize = fs)
        ax3.set_ylabel(r'$\Gamma$', fontsize = fs)
        
        ax1.set_title('Phase $\Gamma$', fontsize = ts, \
                      fontweight = 'bold');
                      
        ax2.set_title('Phase $\Gamma$', fontsize = ts, \
                      fontweight = 'bold');
                      
        ax3.set_title('Phase ($\Gamma$)', fontsize = ts, \
                      fontweight = 'bold');
    ax1.tick_params(labelsize = 'large');
    ax2.tick_params(labelsize = 'large');
    ax3.tick_params(labelsize = 'large');         
    for i in range(len(x_plot)):
        
        xlab        = x_plot[i] * 1e6; # For labeling
        xlab        = np.round(xlab);
        g_drive     = f_drive(x_plot[i]);
        g_wit       = f_wit(x_plot[i]); 
        g_total     = g_drive + g_wit;
        if det == 'bal':
            drive_plot   = np.sin(g_drive);
            wit_plot     = np.sin(g_wit);
            both_plot    = np.sin(g_total);
        elif det == 'cross':
            drive_plot = np.sin(g_drive / 2)**2
            wit_plot   = np.sin(g_wit / 2)**2
            both_plot  = np.sin(g_total / 2)**2
        else:
            drive_plot = g_drive;
            wit_plot   = g_wit;
            both_plot  = g_total;
        
        ax1.plot(t_sigs, drive_plot, label = r'$\Delta$ x = ' \
                                                     + str(xlab) + r'$\mu$m');
        ax2.plot(t_sigs, wit_plot, label = r'$\Delta$ x = ' \
                                                     + str(xlab) + r'$\mu$m');
        ax3.plot(t_sigs, both_plot, label = r'$\Delta$ x = ' \
                                                     + str(xlab) + r'$\mu$m');
                 
        
    #box = ax3.get_position();
    #ax3.set_position([box.x0, box.y0, box.width * 0.8, box.height]);
    
    ax1.legend();
    ax2.legend();
    ax3.legend();
    #ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5));

def signal_difference(drive_sigs, wit_sigs, tsigs, x_off, x_interp, del_t, 
                      det, plot):
    '''
    Function to calculate and plot the signals detected by two EO crystals and
    the signal difference at t = 0 and t = del_t
    
    Parameters:
    -----------
    drive_sigs : array_like
                 2D array of the phase retardation parameter due to the drive
                 bunch. axis 1 = horizontal offset, axis 2 = probe delay time
    wit_sigs  : array_like
                similar array for the phase retardation parameter due to the 
                witness bunch
    tsigs     : array_like
                the probe delay time associated with drive_sigs and wit_sigs
    x_off     : array_like 
                The horizontal offset associated with drive_sigs and wit_sigs
    x_interp  : array_like
                Array of horizontal offsets used in interpolation of 
                drive_sigs and wit_sigs
    del_t     : the longitudinal offset between the drive and witness bunches
    
    det       : str
                The detector scheme 'bal' for balanced detectors, 'cross' for 
                cross polarized, '' for just calculating gamma
    plot      : bool, 
                whether or not to plot
    Returns:
    --------
    r_drive_sig : array_like
                  The combined signal of the right crystal at tsigs = 0 at 
                  each value of x_interp
    r_wit_sig   : array_like
                  The combined signal of the right crystal at tsigs = del_t at
                  each value of x_interp
    l_drive_sig : array_like
                  The combined signal of the left crystal at tsigs = 0 at 
                  each value of x_interp
    l_wit_sig  : array_like
                  The combined signal of the right crystal at tsigs = del_t at 
                  each value of x_interp
    '''
    
    # First create interpolation functions of the signals, compute for 
    # x_interp
    f_drive = interp.interp1d(x_off, drive_sigs, axis = 0);
    f_wit   = interp.interp1d(x_off, wit_sigs, axis = 0);
    
    drive_interp = f_drive(x_interp);
    wit_interp   = f_wit(x_interp);
    
    # preallocate left signals for loop (right signals are left signals 
    # flipped)
    l_drive_sig = np.zeros(len(x_interp));
    l_wit_sig   = np.zeros(len(x_interp));
    for i in range(len(x_interp)):
        g_total = drive_interp[i, :] + wit_interp[i, :];
        
        if det == 'bal':
            fg = interp.interp1d(tsigs, g_total);
            l_drive_sig[i] = np.sin(fg(0));
            l_wit_sig[i]   = np.sin(fg(del_t));
        elif det == 'cross':
            fg = interp.interp1d(tsigs, g_total);
            l_drive_sig[i] = np.sin(fg(0) / 2)**2;
            l_wit_sig[i]   = np.sin(fg(del_t) / 2)**2;
        else:
           fg = interp.interp1d(tsigs, g_total);
           l_drive_sig[i] = fg(0);
           l_wit_sig[i]   = fg(del_t);
    r_drive_sig = np.flip(l_drive_sig, axis = 0);
    r_wit_sig   = np.flip(l_wit_sig, axis = 0);
    # Plot results
    fig1 = plt.figure(figsize = (6,14), dpi = 200);
    ax1  = fig1.add_subplot(211);
    ax2  = fig1.add_subplot(212);

    
    ax1.set_prop_cycle(cy);
    ax1.set_xlabel(r'Horizontal offset [$\mu$m]');
    ax1.set_ylabel('EO Signal [AU]');
    
    ax2.set_prop_cycle(cy);
    ax2.set_xlabel(r'Horizontal offset [$\mu$m]');
    ax2.set_ylabel('EO Signal [AU]');
    

    
    ax1.plot(x_interp * 1e6, r_drive_sig, label = "Right crystal");
    ax1.plot(x_interp * 1e6, l_drive_sig, label = "Left crystal");
    
    ax1.plot(x_interp * 1e6, abs(l_drive_sig - r_drive_sig), '--', \
                       label = "|Difference|");
      
    ax2.plot(x_interp * 1e6, r_wit_sig, label = "Right crystal");
    ax2.plot(x_interp * 1e6, l_wit_sig, label = "Left crystal");
    
    ax2.plot(x_interp * 1e6, abs(l_wit_sig - r_wit_sig), '--',\
                                label = r'|Difference|');
    
    # Titles vary based on detector
    if det == 'bal':
        ax1.set_title(r'sin[$\Gamma$(0)]');
        ax2.set_title(r'sin[$\Gamma$($t_{wit}$)]');
    elif det == 'cross':
        ax1.set_title(r'sin$\frac{\Gamma (0)}{2}$');
        ax2.set_title(r'sin$\frac{\Gamma (t_{wit})}{2}$');
    else:
        ax1.set_title(r'$\Gamma$(0)');
        ax2.set_title(r'$\Gamma$($t_{wit}$)');
    ax1.legend();
    ax2.legend();
    fig1.subplots_adjust(hspace = 0.2);
    return r_drive_sig, r_wit_sig, l_drive_sig, l_wit_sig;

def fwhm(drive_sigs, wit_sigs, x_off, tsigs, x_interp, plot = False):
    '''
    Function to calculate and plot the FWHM of the phase retardation curves 
    vs. the horizontal offset of the electron bunches. 
    
    Parameters:
    -----------
    drive_sigs : array_like
                 2D array of phase retardation due to the drive beam, axis 1 
                 is horizontal offset, axis 2 is time. 
    wit_sigs   : array_like
                 Corresponding 2D array of witness phase retardation. 
    x_off      : array_like
                 Array of horizontal offsets corresponding to the phase 
                 retardation arrays (m). 
    tsigs      : array_like
                 time array (fs) corresponding to the phase retardation arrays
    x_interp   : array_like
                 Array of horizontal offsets to interpolate for. 
    plot       : bool, otpional
                 Whether or not to plot FWHM vs x_interp
    Returns:
    --------
    fwhm_drive : array_like
                 full width half maxes of the phase retardation due to the
                 drive bunch
    fwhm_wit   : array_like
                 full width half maxes of the phase retardation due to the 
                 witness bunch
    '''
        
    # First interpolate drive and witness gammas
    f_drive      = interp.interp1d(x_off, drive_sigs, axis = 0)
    drive_interp = f_drive(x_interp);
    f_wit        = interp.interp1d(x_off, wit_sigs, axis = 0)
    wit_interp   = f_wit(x_interp);
    # Preallocate for loop
    fwhm_drive   = np.zeros(len(x_interp));
    fwhm_wit     = np.zeros(len(x_interp));
    for i in range(len(x_interp)):
        g_drive       = drive_interp[i, :];
        hm_drive      = max(g_drive) / 2;
        fd_interp     = interp.interp1d(tsigs, g_drive);
        fd_interp2    = lambda tsigs: fd_interp(tsigs) - hm_drive;
        root1d        = optimize.newton(fd_interp2, -70), 
        root2d        = optimize.newton(fd_interp2, 30);  
        fwhm_drive[i] = abs(root1d - root2d);
        # Repeat for witness
        g_wit = wit_interp[i, :];
        hm_wit = max(g_wit) / 2;
        fw_interp = interp.interp1d(tsigs, g_wit);
        fw_interp2 = lambda tsigs: fw_interp(tsigs) - hm_wit;
        root1w     = optimize.newton(fw_interp2, 510);
        root2w     = optimize.newton(fw_interp2, 610);
        fwhm_wit[i] = abs(root1w - root2w)
    if plot:
        # 1st plot FWHms for each crystal
        fig = plt.figure(figsize = (6,12), dpi = 200);
        ax1  = fig.add_subplot(211);
        ax1.set_prop_cycle(cy);
        drive_plot = abs(fwhm_drive - np.flip(fwhm_drive, axis = 0))
        ax1.plot(x_interp * 1e6, drive_plot, label = r'|$\Delta$FWHM|');
        ax1.set_title("Drive dominant curve");
        ax1.legend();
        
        ax2 = fig.add_subplot(212);
        ax2.set_prop_cycle(cy);
        wit_plot  = abs(fwhm_wit - np.flip(fwhm_wit, axis = 0));
        ax2.plot(x_interp * 1e6, wit_plot, label = r'|$\Delta$FWHM|');

        ax2.set_title(r'Witness dominant curve');
        ax2.legend();
        
def plot_difference(drive_sigs, wit_sigs, x_off, tsigs, x_interp, det):
    '''
    Function to plot the signal difference in a range around the "peak" 
    (in phase) for horizontal offsets of 100, 75, 50, and 25 microns
    Parameters:
    -----------
    drive_sigs : array_like
                 2D array of drive phase. Axis 1 is horizontal offset, axis 2
                 is time
    wit_sigs   : array_like
                 corresponding 2D array of witness phase.
    x_off      : array_like
                 the horizontal offset corresponding to the phases
    tsigs      : array_like
                 the time array corresponding to the phases
    x_interp   : array_like
                 Array of horizontal offsets to interpolate the phase 
    det        : str
                 Detection scheme, either 'bal' or 'cross', for just phase ''
    Returns:
    --------
    drive_diff : array_like
                 Array of the difference in signal from the drive beam
    wit_diff   : array_like 
                 Array of the difference in signal from the witness beam
    '''
    # Create interpolated arrays
    f_drive      = interp.interp1d(x_off, drive_sigs, axis = 0);
    f_wit        = interp.interp1d(x_off, wit_sigs, axis = 0);
    drive_interp = f_drive(x_interp);
    wit_interp   = f_wit(x_interp);  
    
    # Calculate the region of interest (0 ind is strongest signal);
    i_drive = np.argmax(drive_interp[0, :]);
    i_wit   = np.argmax(wit_interp[0, :]);
    d_roi   = np.argmin(abs(drive_interp[0, i_drive:-1] - drive_sigs[0, 0]));
    d_roi   = d_roi + i_drive;
    # Witness is same width as drive shifted. 
    w_roi1  = int(i_wit - d_roi / 2);
    w_roi2  = int(i_wit + d_roi / 2);
    
    t_drive = tsigs[0:d_roi];
    t_wit   = tsigs[w_roi1:w_roi2];
    
    # preallocate for loop
    drive_diff = np.zeros((4, len(t_drive)));
    wit_diff   = np.zeros((4, len(t_wit)));
    
    x_calc = np.array([-100, -75, -50, -25]) * 1e-6;
    
    fig1 = plt.figure(figsize = (8, 6), dpi = 200);
    fig2 = plt.figure(figsize = (8,6), dpi = 200);
    fs   = 16;
    ts   = 20;
    ax1  = fig1.gca();
    ax2  = fig2.gca();
    ax1.set_prop_cycle(cy);
    ax2.set_prop_cycle(cy);
    ax1.tick_params(labelsize = 'large');
    ax2.tick_params(labelsize = 'large');
    ax1.set_xlabel(r'$\tau$ [fs]', fontsize = fs);
    ax2.set_xlabel(r'$\tau$ [fs]', fontsize = fs);
    ax1.set_ylabel('Signal difference [AU]', fontsize = fs);
    ax2.set_ylabel('Signal difference [AU]', fontsize = fs);
    if det == 'bal' or det == 'cross':
        ax1.set_title("Drive beam signal difference", fontsize = ts);
        ax2.set_title("Witness beam signal difference", fontsize = ts);
    else:
        ax1.set_title("Drive beam phase difference");
    
    
    for i in range(len(x_calc)):
        drive_left  = f_drive(x_calc[i])[0:d_roi];
        drive_right = f_drive(-1 * x_calc[i])[0:d_roi];

        wit_left  = f_wit(x_calc[i])[w_roi1:w_roi2];
        wit_right = f_wit(-1 * x_calc[i])[w_roi1:w_roi2];
                
        if det == 'bal':
            drive_left  = np.sin(drive_left);
            drive_right = np.sin(drive_right);
            wit_left    = np.sin(wit_left);
            wit_right   = np.sin(wit_right);
        elif det == 'cross':
            drive_left  = np.sin(drive_left / 2)**2;
            drive_right = np.sin(drive_right / 2)**2;
            wit_left    = np.sin(wit_left / 2)**2;
            wit_right   = np.sin(wit_right / 2)**2;    
        
            
        drive_diff[i, :] = drive_right - drive_left;
        wit_diff[i, :]   = wit_right - wit_left;
    

    # Normalize signal, adjust subtracton:
    max_drive = np.amax(drive_diff);
    
    drive_diff = drive_diff / np.amax(abs(drive_diff));
    wit_diff   = wit_diff / np.amax(abs(wit_diff));
    for i in range(len(x_calc)):
        xlab = str(abs(np.round(x_calc[i] * 1e6)));
        ax1.plot(t_drive, drive_diff[i, :],\
                 label = r'|$\Delta$x| = ' + xlab + r'$\mu$m');
        ax2.plot(t_wit, wit_diff[i, :],\
                 label = r'|$\Delta$x| = ' + xlab + r'$\mu$m');
    

    ax1.legend();
    ax2.legend();
    
    return drive_diff, wit_diff;
    
    
         
        
        
    
