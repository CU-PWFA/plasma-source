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

def comp_delta(sigA, sigB, g0, method):
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
           Nominal phase retardation with no offset
    method : str
             Method of detection, either "cross" or "bal"
    cutoff : float
             Cutoff percentage of signal to compute delta
    Returns:
    --------
    delta : array_like
          percent offset of the beam 
    """
    # Compute signal sums and differences
    s_diff = sigA - sigB
    s_sum  = sigA + sigB
    R      = s_diff / s_sum
    # Compute delta based on detection method
    if method == "cross":
        delta = R * np.tan(0.5*g0) / g0
        return delta
    elif method == "bal":
        delta = R * np.tan(g0) / g0
        return delta
    else:
        print("Error: unknown detection method")
        return

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

def optimize_r0(I, ti, r0s, setup, verbose = False):
    """
    Function to optimize the crystal beamline spacing for maximum response
    to a beams transverse offset.
    
    Parameters:
    -----------
    I       : array_like, 
              Input current profile, kA
    ti      : array_like
              Time array corresponding to I, s
    r0s     : array_like
              Array of crystal beamline spacing, m 
    setup   : dictionary_object
              EOS-BPM setup dictionary as defined in eo_signal
    verbose : bool, optional
              Whether or not to print loop status
    Returns:
    --------
    g0        : array_like
                Array of maximum phase retardation for each r0
    bal_res   : array_like
                Array of the balanced detectors response to transverse
                offsets
    cross_res : array_like
                Array of the crossed polarizers response to transverse
                offsets
    """
    # Preallocate
    g0 = np.zeros(len(r0s))
    # Compute nominal phase retardation for each setup
    for i in range(len(r0s)):
        if verbose:
            print(i+1, "of", len(r0s))
        r0 = r0s[i]
        setup["r0"] = r0
        E, te = cp.get_E(I, ti, r0)
        # Set probe timing for the setup
        setup["tau"] = te
        sig, t_sig, gamma, t_gamma = eos.E_signal(E, te, setup)
        g0[i] = np.amax(gamma)
    # Compute response
    bal_res = g0 * np.cos(g0)
    cross_res = g0 * np.sin(g0)
    # Set response to nan where phase wrapping would occur
    bal_res[g0 > 0.5*np.pi] = np.nan
    cross_res[g0 > np.pi] = np.nan
    return g0, bal_res, cross_res
