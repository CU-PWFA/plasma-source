from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, epsilon_0
import scipy.interpolate as interpolate
import scipy.io as io
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import interp1d
import sys;
eps0 = epsilon_0;
gamma = 19569.4716; # Nominal FACET-II (E = 10 GeV)
# Simulation module
#import eosim as sim;
# Constants
eps0 = 8.854187817e-12;
# Colors for plotting.
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77','#CC6677', '#882255',  '#AA4499'];
cy = cycler('color', plot_colors)

def get_peaks(signal, t, dt = 0.25e-12):
    '''
    Function to find the drive and witness peaks of a given signal

    Parameters:
    -----------
    signal : array_like
             Array of some beam parameter (current, E-field) that is dual
             peaked
    dt     : float
             Minimal temporal separation of the bunches (in steps)
    Returns:
    --------
    peaks : array_like
            Array containing the indices of the two peaks
    '''
    max_sig = max(signal)
    # Start cutoff at 5 percent the max signal
    cutoff = 0.05
    # Find distance
    t1   = t[0] + dt
    sep  = np.argmin(abs(t-t1))

    # Start recursion
    peaks  = find_peaks(signal, height = cutoff * max_sig, distance = sep)[0]
    npeaks = len(peaks)
    while len(peaks) > 2:
        cutoff = cutoff + 0.025
        # sometimes the 2nd peak repeats
        peaks = find_peaks(signal, height = cutoff * max_sig, distance = sep)[0]
        if len(peaks) == 1: 
            cutoff = cutoff - 0.025
            peaks = find_peaks(signal, height = cutoff * max_sig, \
                               distance = sep)[0]
            peaks = [peaks[0], peaks[1]]
            break
        npeaks == len(peaks)
    return peaks

def get_current(ind, comp = "Desktop"):
    '''
    Loads in a current profile and calculate longitudinal offset

    Parameters:
    -----------
    ind     : int
              The data set index (0-9) *ignore 1, bad set*
    comp    : str, optional
              Computer being used, Desktop or Vela

    Returns:

    I_ka : array_like
           The beam current profile in kA
    ti   : array_like
           Time array corresponding to I_kA in s
    p2p  : float
           Peak to peak spacing of the bunches in seconds
    --------
    '''
    
    # Create file path based on compute
    if comp == "Desktop":
        fpath = "/home/keenan/plasma-source/khuntstone/eos_bpm/studies/" + \
                "current_profiles/"
    elif comp == "Vela":
        fpath = "/home/cu-pwfa/CU-PWFA/plasma-source/khuntstone/eos_bpm/" + \
                "studies/current_profiles/"
    if ind < 10:
        ipdz  = io.loadmat(fpath + "ipdz_for_keenan.mat")["ipdz_for_keenan"]
        shots = io.loadmat(fpath + "shots_for_keenan.mat")["shots_for_keenan"]
    else:
        ind   = ind - 10
        ipdz  = io.loadmat(fpath + "ipdz_forKeenan.mat")["ipdz"]
        shots = io.loadmat(fpath + "ipcurrent_forKeenan.mat")["ipcurrent"]
    dz    = ipdz[ind][0]
    I_ka  = shots[ind]; # kA
    z     = np.array([0 + (i * dz) for i in range(len(I_ka))])
    dzp   = abs(z[0] - z[1])
    ti    = (z * 1e-6) / c
    dt    = (dz * 1e-6) / c

    tsep   = 0.3e-12
    t1     = np.argmin(abs(ti))
    t2     = np.argmin(abs(ti - tsep))
    tsteps = t2 - t1
    # center drive on 0
    peaks = get_peaks(I_ka, ti)
    ti    = ti - ti[peaks[0]]
    p2p   = ti[peaks[1]] - ti[peaks[0]] 

    # Now compute electric field

    return I_ka, ti, p2p
    

def get_E(I, ti, r):
    '''
    Computes the E-field from a current profile

    Parameters:
    -----------
    ind    : int
             The data set index
    r      : float
             Distance from the beam to the crystal (m)
    I      : array_like
             The beam current profile in kA
    ti     : array_like
             Time array corresponding to I in s
    Returns:
    --------
    E    : array_like
           The beam electric field in V/m
    ze   : array_like
           Longitudinal distance array corresponding to E in s
    te   : array_like
           Time array corresponding to E in s
    --------
    '''

    dti = abs(ti[0] - ti[1])

    Er = (e * gamma) / (4 * np.pi * eps0)
    Er = Er * r / ((r**2 + gamma**2 * c**2 * ti**2)**(3/2))
    N  = (I * 1e3) * dti / e # Longitudinal electron distribution
    E  = np.convolve(Er, N)
    te = np.array([0 + (i * dti) for i in range(len(E))])
    # Ensure peak alignment
    peaks = get_peaks(E, te)
    te    = te - te[peaks[0]]
    return E, te

def extend_I(I, ti, t_split, sigma = 0, scaling = 0, method = 'Gauss'):
    """
    Function to increase the longitudinal rms of a single bunch by convolving
    it with a Gaussian.
    
    Parmeters:
    ----------
    I       : array_like
              The current profile
    ti      : array_like
              The time array corresponding to I
    t_split : float
              Time at which to separate the currents
    sigma   : float, optional
              The rms of the Gaussian (used for Gauss method)
    scaling : float, optional
              The scaling factor to extend time by (used for Scaling method)
    Returns:
    --------
    I_conv : array_like
             New current profile with larger rms, corresponding to ti
    """
    t_shift = 50e-15
    if method == 'Gauss':
        def gaussian(x, A, sigma):
            return A * np.exp(-(x/sigma)**2)
            
        I_drive, I_wit = split_I(I, ti, t_split)
        fg      = gaussian(ti, max(I_drive), sigma)
        I_conv  = np.convolve(I_drive, fg)
        I_conv  = max(I)/max(I_conv) * I_conv
        t_conv  = np.linspace(ti[0], ti[-1], len(I_conv))
        f_wit     = interp1d(ti + t_shift, I_wit, \
                             bounds_error = False, fill_value = 0)
        I_wit_int = f_wit(t_conv)
        fI        = interp1d(t_conv, I_conv + I_wit_int)
        I_conv    = fI(ti)
        return I_conv
    elif method == 'scale':
        I_drive, I_wit = split_I(I, ti, t_split)
        t_extend  = ti * scaling
        f_wit     = interp1d(ti + t_shift, I_wit, \
                             bounds_error = False, fill_value = 0)
        I_wit_int = f_wit(t_extend)
        fI        = interp1d(t_extend, I_drive + I_wit_int)
        I_extend  = fI(ti)
        return I_extend

def smooth_I(I, window = 51, poly = 3, double = False):
    '''
    Function to smooth beam current (somewhat like utilizing the 
    laser heater). Smoothing algorithm is Savitzky-Golay filter which fits 
    succesive subsets of adjacent data with a low-degree polynomial.
    Parameters:
    -----------
    I      : array_like
             Input current
    window : int, optional
             The window size (how many adjacent points), must be odd integer
    poly   : int, optional
             The polynomial order. 
    Returns:
    --------
    I_smooth : array_like
               The smoothed current profile
    '''
    
    
    I_smooth = savgol_filter(I, window, poly);
    # No negative current
    I_smooth[I_smooth < 0] = 0;
    if double:
        I_smooth = savgol_filter(I_smooth, window, polyorder = 1);
        I_smooth[I_smooth < 0] = 0;
    return I_smooth;
