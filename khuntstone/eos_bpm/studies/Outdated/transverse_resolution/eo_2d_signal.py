#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelling the 2D EO signal to look at detector response to vertical offsets.
Created on Tue Sep 17 11:08:50 2019

@author: keenan
"""

# Python modules
import numpy as np;
import matplotlib.pyplot as plt;
from scipy.constants import c, epsilon_0;
from time import time;
eps0 = epsilon_0;
#from scipy.interpolate import interp1d;
import sys;
sys.path.insert(0, "../../python");
# Custom modules
from crystal import crystal;
from ebeam import ebeam;
from laser import laser;
import phase_retard as pr;
from plotting import makefig;
import thz;

def get_2D_gamma(drive, wit, probe, angle, cry, d, tau, y, x, ny, \
                 nslice = 100):
    """
    Compute the 2D phase retardation of the EO crystal
    Parameters:
    -----------
    drive  : dictionary
             Dictionary of keys 'Q', 't', 'del_z', 'sigz' with corresponding
             values for drive bunch charge, time domain, spatial offset, and
             rms length (SI base units)
    wit    : dictionary
             Dictionary of keys 'Q', 't', 'del_z', 'sigz' with corresponding
             values for witness bunch charge, time domain, spatial offset, and
             rms length (SI base units)
    probe  : object
             Instance of laser class w/ probe pulse parameters
    angle  : float
             Angle of incidence of the probe (rad)
    cry    : object
             Instance of the EO crystal class
    d      : float
             The crystal thickness (m)
    tau    : array_like
             Signal time array (s)
    y      : float
             vertical height of the crystal
    x      : float
             The horizontal crystal offset from the beamline
    ny     : number of heights to compute signal for
    nslice : int, optional
             The number of slices to divide the crystal into
    
    Returns:
    --------
    gamma : array_like
            The 2D phase retardation
    
    """
    # Slice crystal
    j     = np.arange(1, nslice, 1);
    dz    = d / nslice;
    d_arr = (j - 0.5) * dz;
    # Create vertical array (signal is symettric, save on computation by 
    # only taking +y)
    y_arr = np.linspace(0, y, ny);
    # Initialize gamma
    gamma = np.zeros((len(tau), len(y_arr)));
    
    # Propagate the pulses through the crystal (only the temporal part)
    t_off = wit['del_z'] / c;
    wsigt = wit['sigz'] / c;
    dsigt = drive['sigz'] / c;
    wit_t = np.exp(-(wit['t'] + t_off)**2 / (2 * wsigt**2));
    drive_t = np.exp(-(drive['t'])**2 / (2 * dsigt**2));
    
    FEr_drive, f_drive = thz.raw_field(drive_t, drive['t']);
    FEr_wit, f_wit     = thz.raw_field(wit_t, wit['t']); 
    
    print("Propagating")
    Ec_drive, tt_drive = thz.cry_field(drive['t'], FEr_drive, f_drive, d, \
                                       probe, cry, nslice = nslice);                               
    Ec_wit, tt_wit     = thz.cry_field(wit['t'], FEr_wit, f_wit, d, \
                                       probe, cry, nslice = nslice);   
    print("Looping");
    for i in range(len(tau)):
        if (i+1)%10 == 0:
            print(np.round((i + 1) / len(tau) * 100), "%");
        xp = x + tau[i] * c /  np.tan(angle); # Horizontal location of pulse
        for j in range(len(y_arr)):
            # Calculate amplitudes
            r = np.sqrt(xp**2 + y_arr[j]**2);
            drive_E0 = drive['Q'] / (2 * np.pi * eps0 * np.sqrt(2 * np.pi) \
                             * c * r * dsigt);
            wit_E0   = wit['Q'] / (2 * np.pi * eps0 * np.sqrt(2 * np.pi) \
                             * c * r * wsigt);
            
            # Birefringence at the spatial location
            bir_drive = drive_E0 * Ec_drive;
            bir_wit   = wit_E0 * Ec_wit;
            
            g_drive, t_drive = pr.phase_retard(bir_drive, tt_drive, d_arr, \
                                               tau, probe, cry, 'spatial', \
                                               psi = angle);
            g_wit, t_wit = pr.phase_retard(bir_wit, tt_drive, d_arr, tau, \
                                               probe,cry, 'spatial', \
                                               psi = angle);
            gamma[:, j] = g_drive + g_wit;
            
    
    return gamma;

def plot_2d_gamma(gamma, tau, angle, x, y):
    fig, ax = makefig(xlab = 'x [mm]', ylab = 'y [mm]');
    xlim = np.array([tau[0], tau[-1]]) * c / np.tan(angle) + x;
    xlim = xlim * 1e3;
    ylim = [0, y * 1e3];
    ext  = [xlim[0], xlim[1], ylim[0], ylim[1]];
    img  = ax.imshow(np.transpose(gamma), cmap = 'CMRmap', extent = ext, 
                     aspect = 'auto')
    cb   = plt.colorbar(mappable = img);
    plt.show();
# Run below
dev = True;
#dev = False;
if dev:
    Q_drive = 1.5e-9;
    Q_wit   = 0.5e-9;
    sigz    = 5.2e-6;
    sigt    = sigz/c;
    del_z   = 150e-6;
    N       = 8000;
    tb      = np.linspace(-N * sigt / 2, N * sigt / 2, N);
    drive   = {'Q' : Q_drive,
               'sigz' : sigz,
               't' : tb, 
               'del_z' : 0};
    wit     = {'Q' : Q_wit,
               'sigz' : sigz,
               't' : tb, 
               'del_z' : del_z};
    #r = 5e-3;
    #beam1 = ebeam(drive);
    #beam2 = ebeam(wit);
    #beam1.get_Er(r);
    #beam2.get_Er(r);
    #Er = beam1.Er + beam2.Er;
    #Er_norm = Er / max(Er);
    #FEr, f = thz.raw_field(Er, beam1.t);
    #FErn, fn = thz.raw_field(Er_norm, beam1.t);
    
    # Crystal
    cry     = crystal('GaP');
    d       = 100e-6;
    y       = 5e-2;
    ny      = 10;
    x       = 1e-3;
    # Probe
    y0    = 800e-9;
    dy    = 27e-9;
    tp    = 30e-15;
    probe = laser({'y0' : y0, 'dy' : dy, 'tp' : tp});
    
    #Ec, tt = thz.cry_field(tb, FEr, f, d, probe, cry);
    #Ecn, ttn = thz.cry_field(tb, FErn, fn, d, probe, cry);
    
    angle = 15 * np.pi / 180;
    tau   = np.linspace(0, 2000, 100) * 1e-15; 
    start = time();
    gamma = get_2D_gamma(drive, wit, probe, angle, cry, d, tau, y, x, ny);
    print(time() - start);
    plot_2d_gamma(gamma, tau, angle, x, y);
    
    
    

    