#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for performing EOS-BPM optimization parameter scans. 

@author: keenan
"""
# Standard python imports
from functools import partial
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
from scipy.interpolate import interp1d
import sys
import time
sys.path.insert(0, "/home/keenan/eos_bpm/python")
# EOS-BPM modules
import eo_signal as eos
from plotting import makefig

def get_err(ind, d, th, fpath, plot = False, nmax = 4):
    """
    Function to get the longitudinal peak to peak error of a given EOS-BPM 
    setup of crystal thickness and angle (other parameters are fixed).
    
    Parameters:
    -----------
    ind   : int,
            Index of the current profile to use
    d     : float,
            Crystal thickness (m)
    th    : float,
            Probe crossing angle (deg.)
	fpath : str, 
	        Full path to the current profile and dz .mat files
    Returns:
    --------
    error_t : float
              Longitudinal error in s
    error_p : float
              Longitudinal error as a %
    """
    # Get current and signal 
    setup = {"ctype"  : "GaP", 
             "d"      : d, 
             "y0"     : 800e-9,
             "tp"     : 30e-15,
             "angle"  : th, 
             "r0"     : 2.5e-3,
             "method" : "cross",
             "fpath"  : fpath,
             "tilt"   : 0,
             "th"     : 0,
             "nslice" : 100,
             "plot"   : False, 
             }
    I, ti, Idt, sig, t_sig, gamma, t_gamma = eos.get_signal(ind, setup)
    Sdt = eos.peak2peak(sig, t_sig)
    error_t = abs(Sdt - Idt)
    error_p = error_t/ Idt
    return error_t, error_p


def scan2D(ds, ths, fpath, N = 3134, N_p = 4):
    """
    Function to scan over crystal thicknesses and angles and compute the
    average error. 
    
    Parameters:
    -----------
    ds    : array_like
            Array of crystal thicknesses (m)
    ths   : array_like
            Array of probe crossing angles (deg.)
    fpath : string, 
            full path to the current profile and dz .mat files
    N     : int, optional  
            Number of current profiles to use, default = 3134 (maximum)
    N_p   : int, optional
            Number of processors to use, default = 4
          
    Returns:
    --------
    errors_t : array_like
               2D array of longitudinal errors as a time
    errors_p : array_like
               2D array of longitudinal errors as a %
    """
    
    # Create array of indices to scan across
    inds     = np.arange(0, N, 1)
    # Preallocate for loops
    errors_p = np.zeros((len(ds), len(ths)))
    errors_t = np.zeros((len(ds), len(ths)))
    start    = time.time()
    for i in range(len(ds)):
        print(i+1, "of", len(ds))
        for j in range(len(ths)):
            pool  = Pool(N_p)
            errf  = partial(get_err, d = ds[i], th = ths[j], fpath = fpath)
            errt, errp   = zip(*pool.map(errf, inds))
            pool.close()
            pool.join()
            errors_t[i,j] = np.nanmean(errt)
            errors_p[i,j] = np.nanmean(errp)      
    print("Completed in", time.time() - start, "seconds")
    
    return errors_t, errors_p
    
def get_peak(ind, d, th, r0, dx, fpath):
    """
    Function to get the peak signal of a given EOS-BPM setup of crystal 
    thickness, probe crossing angle, and crystal-beamline distance. 
    
    Parameters:
    -----------
    ind : int
          Index of the current profile to use
    d   : float
          Crystal thickness (m)
    th  : float
          Probe crossing angle (deg.)
    r   : float
          Crystal-beam distance (m)
    dx  : float
          Beam offset from propagation axis
          
    Returns:
    --------
    S_peak : float
             The maximum signal of the EOS-BPM setup
    """
    setup = {"ctype"  : "GaP", 
             "d"      : d, 
             "y0"     : 800e-9,
             "tp"     : 30e-15,
             "angle"  : th, 
             "r0"     : r0 + dx,
             "method" : "cross",
             "fpath"  : fpath,
             "tilt"   : 0,
             "th"     : 0,
             "nslice" : 100,
             "plot"   : False, 
             }
    I, ti, Idz, sig1, t_sig1, gamma, t_gamma = eos.get_signal(ind, setup)
    setup["r0"] = r0 - dx
    I, ti, Idz, sig2, t_sig2, gamma, t_gamma = eos.get_signal(ind, setup)
    sig = sig1 - sig2
    maxS = np.nanmax(sig)
    minS = -1.0 * np.nanmin(sig)
    if minS > maxS:
        return -1.0*minS
    else:
        return maxS

def scan1D(r0s, dx, d, th, fpath, N = 3134, N_p = 4):
    """
    Function to perform a 1D parameter scan of crystal-beamline distance and 
    compute the average peak signal vs. transverse offset
    
    Parameters:
    -----------
    r0s : array_like
          Array of crystal-beamline distances to scan (m)
    dx  : array_like
          Array of electron bunch transverse offsets to scan for each r0 (m)
    d   : float
          Crystal thickness (m)
    th  : float
          Probe crossing angle (deg.)
    N   : int, optional
          Number of current profiles to use, default = 3134 (max)
    N_p : int, optional
          Number of processors to use, default = 4
    
    Returns:
    --------
    S_peaks : array_like
              2D array of peak signal vs transverse offset for each r0
    """
    
    # Create array of indices to scan
    inds = np.arange(0, N, 1)
    # Preallocate for loo[]
    S_peaks  = np.zeros([len(r0s), len(dx)])
    start = time.time()
    for i in range(len(r0s)):
        print(i+1, "of", len(r0s))
        for j in range(len(dx)):
            smax = partial(get_peak, d = d, th = th, r0 = r0s[i], dx = dx[j], \
            	           fpath = fpath)
            pool  = Pool(N_p)
            peaks = pool.map(smax, inds)
            pool.close()
            pool.join()
            S_peaks[i,j] = np.nanmean(peaks)
    print("Completed in", time.time()-start, "seconds")
    return S_peaks

# Plotting functions

def plot_2D(errors_t, errors_p, ds, ths, cm = "CMRmap"):
    """
    Function to plot the results of scan2D. 
    
    Parameters:
    -----------
    errors_t : array_like
               Array as returned from scan2D
    errors_p : array_like
               Array as returned from scan2D
    ds       : array_like
               Array of crystal thicknesses input into scan2D
    ths      : array_like
               Array of angles input into scan2D
    cm       : str, optional
               Colormap to use in plotting
    """
    
    # Compute extent (needs evenly spaced ds and ths)
    dd    = (ds[1]-ds[0])*1e6
    dmin  = ds[0]*1e6 - dd/2
    dmax  = ds[-1]*1e6 + dd/2
    dt    = ths[1] - ths[0]
    tmin  = ths[0] - dt/2
    tmax  = ths[-1] + dt/2
    ext       = [dmin, dmax, tmin, tmax]
    
    # Transpose error arrays and flip for plotting (convert errors_t to fs)
    errp_plot = np.flipud(np.transpose(errors_p))
    errt_plot = np.flipud(np.transpose(errors_t*1e15))
    
    # Plot percentage
    fig1, ax1 = makefig(xlab = r'd [$\mu$m]', ylab = "Angle [deg.]")
    ax1.set_xticks(ds*1e6)
    img1      = ax1.imshow(errp_plot*100, cmap = cm, \
                           extent = ext, aspect = 'auto')
    plt.colorbar(mappable = img1, label = 'Avg. Error [%]')
    # Plot distance
    fig2, ax2 = makefig(xlab = r'd [$\mu$m]', ylab = "Angle [deg.]")
    ax2.set_xticks(ds*1e6)
    img2      = ax2.imshow(errt_plot, cmap = cm, \
                           extent = ext, aspect = 'auto')
    plt.colorbar(mappable = img2, label = r'Avg. Error [fs]') 
    
    plt.show()

def polyfit(x, y, degree):
    """
    Fits a polynomial to a data set, useful for analyzing scan1D results. 
    
    x      : array_like
             Input data 
    y      : array_like
             Input data
    degree : int
             Degree of the polynomial
             
    Returns:
    --------
    results : dict object
              Polynomial coefficient and determination of the fit
    """
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                        
    ybar = np.sum(y)/len(y)          
    ssreg = np.sum((yhat-ybar)**2)   
    sstot = np.sum((y - ybar)**2)
    results['determination'] = ssreg / sstot
    return results

def plot_1D(r0s, dx, S_peaks, interp = False, norm = "relative"):
    """
    Function to plot the signal vs transverse offset for each r0 in scan_1D. 
    Also computes and plots the overall drop in signal and the linearity of 
    S_peaks vs dx. 
    
    Parameters:
    -----------
    r0s     : array_like
              Array of crystal beam distances input into scan1D
    dx      : array_like
              Array of transverse offsets input into scan1D
    S_peaks : array_like
              Array of peak signals returned from scan1D
    interp  : bool, optional
              Whether or not interpolate input
              
    Returns:
    --------
    r2     : array_like
             Determination of linear fit for S_peaks vs dx for each r0 value
    s_drop : array_like
             Drop in signal over dx for each r0 value
    """
    if interp:
        x_int  = np.linspace(dx[0], dx[-1], 1000)
        Sp_int = np.zeros((len(r0s), len(x_int)))
        for i in range(len(r0s)):
            fSp = interp1d(dx, S_peaks[i, :])
            Sp_int[i,:] = fSp(x_int)
        Sp_use = Sp_int
        x_use  = x_int
    else:
        Sp_use = S_peaks
        x_use  = dx
    # Signal vs offset
    if norm == "absolute":
        scale  = np.nanmax(Sp_use)
    fig1, ax1 = makefig(xlab = r'Offset [$\mu$m]', ylab = 'Peak Signal [AU]')
    for i in range(len(r0s)):
        lab = r'$r_0$ = ' + str(np.round(r0s[i]*1e3,2)) + "mm"
        if norm == "relative":
            scale = np.nanmax(Sp_use[i,:])
        Sp_use[i, :] = Sp_use[i,:] / scale
        ax1.plot(x_use*1e6, Sp_use[i,:], label = lab)
    ax1.legend()
    
    # Determination and Signal drop
    r2     = np.zeros(len(r0s))
    s_drop = np.zeros(len(r0s))
    for i in range(len(r0s)):
        results   = polyfit(x_use, Sp_use[i,:], 1)
        r2[i]     = results['determination']
        s_drop[i] = results['polynomial'][0]*(-100e-6) # magnitude / micron
        
    fig2 = plt.figure(figsize = (8, 12), dpi = 200)
    ax2  = fig2.add_subplot(211)
    ax2.set_ylabel("Determination")
    ax2.plot(r0s*1e3, r2, 'o', markersize=5)
    ax3 = fig2.add_subplot(212)
    ax3.set_xlabel(r'$r_0$ [mm]')
    ax3.set_ylabel(r'Signal Drop / $\mu$m  [AU]')
    ax3.plot(r0s*1e3, abs(s_drop), 'o', markersize=5)
    fig2.align_ylabels()
    plt.subplots_adjust(hspace=0)
    plt.show()
    
    return r2, s_drop

    
