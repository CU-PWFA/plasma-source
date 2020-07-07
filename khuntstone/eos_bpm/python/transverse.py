#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 15:51:13 2020

Contains useful functions for optimizing transverse EOS-BPM
@author: keenan
"""

import matplotlib.pyplot as plt 
#try:
#    plt.style.use('huntstone')
#except:
#    plt.style.use('dark_background')
import numpy as np
import sys
sys.path.insert(0,"../../python/")
import currents as cp
import eo_signal as eos
from plotting import makefig

def cry_sigs(offsets, setup):
    """
    Function to compute the peak drive signal from the two EO crystals of given 
    setup over multiple transverse offsets.

    Parameters:
    -----------
    offsets : array_like
              Array of beam's transverse offsets
    setup   : dictionary
              Dictionary containing keys 'ind', 'ctype', 'd', 'y0', 'tp', 
              'angle', 'r0', and 'method'. Corresponding to the current profile
              index, the crystal type, crystal thickness, probe central 
              wavelength, the probe pulse duration, the crossing angle, 
              the crystal-beamline offset, and the detection method. 
    Returns:
    --------
    sigA : array_like
           Max signal in crystal A corresponding to each offset
    sigB : array_like
           Max signal in crystal B corresponding to each offset
    """
    # Extract parameters
    ind = setup["ind"]
    r0  = setup["r0"]
    # Preallocate for loop
    sigA = np.zeros(len(offsets))
    sigB = np.zeros(len(offsets))
    gA   = np.zeros(len(offsets))
    gB   = np.zeros(len(offsets))
    # Loop through offsets and compute signal
    for i in range(len(offsets)):
        setup["r0"] = r0 + offsets[i]
        simA = eos.get_signal(ind, setup)
        setup["r0"] = r0 - offsets[i]
        simB = eos.get_signal(ind, setup)
        
        sigA[i] = max(simA[3])
        gA[i]   = max(simA[5])
        sigB[i] = max(simB[3])
        gB[i]   = max(simB[5])
    return sigA, sigB, gA, gB

def comp_delta(sigA, sigB, g0):
    """
    Function to compute delta (percent offset) from the two crystal signals
    and average peak 

    Parameters:
    -----------
    sigA : array_like or float
           Crystal A signal
    sigB : array_like or float
           Crystal B signal
    g0   : array_like or float
           Nominal G0

    Returns:
    --------
    order1  : array_like or float
             delta computed by 1st order taylor approximation
    order2  : array_like or float
              delta as computed by 2nd order taylor approximation
    taylor2 : array_like
              delta computed as a taylor expansion of order 2
    """
    # Compute signal sums and differences
    s_diff = sigA - sigB
    s_sum  = sigA + sigB
    # First order
    num    = s_diff * np.tan(0.5 * g0)
    den    = s_sum * g0
    order1 = -num / den
    # Second order
    a      = (s_diff / s_sum) * (g0**2 / 4) / np.tan(g0)
    b      = g0
    c      = (s_diff / s_sum) * 0.5 * np.tan(0.5 * g0)
    order2 = -b + np.sqrt(b**2 - 4 * a * c)
    order2 = order2/(2*a)
    order2[np.isnan(order2)] = 0
    # Taylor of second order
    taylor2 = -(np.tan(0.5*g0) * (s_diff / s_sum)) / (2 * g0)
    return order1, order2, taylor2

def split_offsets(d_offsets, w_offsets, setup):
    """
    Function to compute the drive and witness peak signal in each EO crystal
    when the bunches have a transverse offset between them. 

    Parameters:
    -----------
    d_offsets : array_like
                Array of drive beam offsets from the beamline
    w_offsets : array_like
                Array of witness beam offsets from the beamline
    setup     : dictionary
                Dictionary object as described above
    
    Returns:
    --------
    dsigA : array_like
            Array of peak drive signal in A for every offset combination
    dsigB : array_like
            Array of peak drive signal in B for every offset combination
    wsigA : array_like
            Array of peak witness signal in A for every offset combination
    wsigB : array_like
            Array of peak witness signal in B for every offset combination 
    """
    # Extract parameters
    ind    = setup["ind"]
    r0     = setup["r0"]
    # Get current
    I, ti, p2p  = cp.get_current(ind, setup["fpath"])
    # Preallocate for loop
    dsigA = np.zeros((len(d_offsets), len(w_offsets)))
    dsigB = np.zeros((len(d_offsets), len(w_offsets)))
    wsigA = np.zeros((len(d_offsets), len(w_offsets)))
    wsigB = np.zeros((len(d_offsets), len(w_offsets)))
    for i in range(len(d_offsets)):
        ddx = d_offsets[i]
        for j in range(len(w_offsets)):
            wdx = w_offsets[j]
            # Get individual electric fields
            EA_drive, te = cp.get_E(I, ti, r0 + ddx)
            EA_wit, te   = cp.get_E(I, ti, r0 + wdx)
            EB_drive, te = cp.get_E(I, ti, r0 - ddx)
            EB_wit, te   = cp.get_E(I, ti, r0 - wdx)
            # SPlit
            ind          = np.argmin(abs(te - p2p/2))
            EA           = np.zeros(len(EA_drive))
            EB           = np.zeros(len(EA_drive))
            EA[0:ind]     = EA_drive[0:ind]
            EA[ind:-1]    = EA_wit[ind:-1] 
            EB[0:ind]     = EB_drive[0:ind]
            EB[ind:-1]    = EB_wit[ind:-1] 
            
            # Get signal
            simA = eos.E_signal(EA, te, setup)
            simB = eos.E_signal(EB, te, setup)
            
            dsigA[i, j] = max(simA[0])
            dsigB[i,j]  = max(simB[0])
            ind = np.argmin(abs(simA[1] - p2p))
            wsigA[i,j] = max(simA[0][ind:-1])
            wsigB[i,j] = max(simB[0][ind:-1])
    return dsigA, dsigB, wsigA, wsigB

def plot_deltas(r0, offsets, d_comp1, d_comp2):
    """
    Function to plot the comptued delta (converted to distance offset)

    Parameters:
    -----------
    r0      : float,
              Nominal crystal-beamline distance
    offsets : array_like,
              Array of transverse offsets 
    d_comp1 : array_like
              1st order computed delta
    d_comp2 : array_like
             2nd order computed delta
    """

    fig, ax = makefig(xlab = r'Drive offset $[\mu$m]',\
                      ylab = r'Computed offset [$\mu$m]')
    ax.plot(offsets * 1e6, -d_comp1 * r0 * 1e6, label = "1st order")
    ax.plot(offsets* 1e6, d_comp2 * r0 * 1e6, label = "2nd order")
    plt.legend()
    plt.show()