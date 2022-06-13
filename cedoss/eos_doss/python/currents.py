from cycler import cycler
import numpy as np
from scipy.constants import c, e, epsilon_0
import scipy.io as sio
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import interp1d
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
        peaks = find_peaks(signal, height = cutoff * max_sig, \
                           distance = sep)[0]
        if len(peaks) == 1: 
            cutoff = cutoff - 0.025
            peaks = find_peaks(signal, height = cutoff * max_sig, \
                               distance = sep)[0]
            peaks = [peaks[0], peaks[1]]
            break
        npeaks == len(peaks)
    return peaks

def get_current(ind, fpath):
    '''
    Loads in a current profile and calculate longitudinal offset

    Parameters:
    -----------
    ind     : int
              The data set index (0-9)
    fpath   : string
              The full path to the current profile and ipdz .mat files

    Returns:

    I_ka : array_like
           The beam current profile in kA
    ti   : array_like
           Time array corresponding to I_kA in s
    p2p  : float
           Peak to peak spacing of the bunches in seconds
    --------
    '''
    
    # Set file path
    #fpath = "/home/keenan/eos_bpm/current_data/"
    #fpath = "/home/keenan/eos_bpm/khuntstone/" + "current_profiles/"
    if ind < 10:
        ipdz  = sio.loadmat(fpath + "dz1.mat")["ipdz_for_keenan"]
        shots = sio.loadmat(fpath + "current1.mat")["shots_for_keenan"]
    else:
        ind   = ind - 10
        ipdz  = sio.loadmat(fpath + "dz2.mat")["ipdz"]
        shots = sio.loadmat(fpath + "current2.mat")["ipcurrent"]
    dz    = ipdz[ind][0]
    I_ka  = shots[ind]; # kA
    #print(I_ka)
    #print(np.shape(I_ka))
    #import matplotlib.pyplot as plt
    #plt.plot(I_ka)
    z     = np.array([0 + (i * dz) for i in range(len(I_ka))])
    ti    = (z * 1e-6) / c


    # center drive on 0
    peaks = get_peaks(I_ka, ti)
    ti    = ti - ti[peaks[0]]
    p2p   = ti[peaks[1]] - ti[peaks[0]] 

    # Now compute electric field

    return I_ka, ti, p2p
    

def get_E(I, ti, r0):
    '''
    Computes the E-field from a current profile

    Parameters:
    -----------
    I      : array_like
             The beam current profile in kA
    ti     : array_like
             Time array corresponding to I in s
    r0     : float or array_like
             Distance from the beam to the crystal or temporal
             profile of distance from the beam to the crystal (m)
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

    #dti = abs(ti[0] - ti[1])

    #Er = (e * gamma) / (4 * np.pi * eps0)
    #dx = c * (ti - ti[0]) * np.tan(tilt)
    #r  = r0 + dx
    #Er = Er * r / ((r**2 + gamma**2 * c**2 * ti**2)**(3/2))
    #N  = (I * 1e3) * dti / e # Longitudinal electron distribution
    #E_int  = np.convolve(Er, N)
    #te_int = np.array([0 + (i * dti) for i in range(len(E_int))])
    # Ensure peak alignment
    #peaks  = get_peaks(E_int, te_int)
    #te_int = te_int - te_int[peaks[0]]
    
    dti  = ti[1]-ti[0]
    dzi  = c * dti
    ne   = I * 1e3 * dti / e
    Ne   = np.sum(ne) # number of particles in the beam
    Qe   = Ne * e # Total charge of the beam
    lamz = I * 1e3 * dti / (Qe * c) # Distribution (m^-1)
    lamz = lamz / np.sum(lamz * dzi)
    E    = 2 * Qe * lamz / (4 * np.pi * epsilon_0 * r0)
    fE     = interp1d(ti, E, bounds_error = False, fill_value = 0)
    te_int = np.linspace(-3, 3, 1000)*1e-12
    E_int  = fE(te_int)
    return E_int, te_int


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

def extractBeam(ind, fpath, onlyI = False):
    """
    Function to extract the dictionaries "data" and "scanvals" 
    describing the SFQED electron bunches as provided by SLAC.

    Parameters:
    -----------
    ind   : int
            Current profile index, 0 to 50
    fpath : string
            Full file path to the .mat workspace file
    onlyI : bool, optional
            Whether or not to only extract the current profile, 
            defualt = False
    Returns:
    --------
    data     : dict
               dictionary object describing the electron beam
    scanvals : dict
               dictionary object describing the linac values to create the 
               beam
    """

    if onlyI:
        mat_contents = sio.loadmat(fpath + "scanDataKeenan_072020.mat")
        oct_struct   = mat_contents["scanData"]
        val          = oct_struct[0, ind]
        I            = val["data"]["I"][0][0]
        dzi          = val["data"]["dz"][0][0][0][0]
        z_arr        = np.linspace(0, dzi*len(I), len(I))
        t_arr        = z_arr / c
        ti           = t_arr - t_arr[np.argmax(I)]
        return np.squeeze(I), ti
    else:
        mat_contents = sio.loadmat(fpath + "scanDataKeenan_072020.mat")
        oct_struct = mat_contents["scanData"]
        # Create dictionaries for beam values
        data       = {}
        scanvals   = {}
        beamparams = {}
        beam       = {}
        bunch      = {}
        beam_keys  = ["sigy", "sigx", "rmsx", "rmsy", "sigE", "sigz", "pkI", \
                      "rmsz", "rmsE", "nx", "ny", "centroidx"]
        # Extract data
        val = oct_struct[0, ind]
        data["scanparams"] = val["data"]["scanparams"][0][0][0]
        data["I"]          = val["data"]["I"][0][0]
        data["Eprof"]      = val["data"]["Eprof"][0][0][0]
        data["dz"]         = val["data"]["dz"][0][0][0][0]
        data["XEdges"]     = val["data"]["XEdges"][0][0][0]
        data["YEdges"]     = val["data"]["YEdges"][0][0][0]
        data["N"]          = val["data"]["N"][0][0]
        #Extract beamparams
        for i in range(12):
            beamparams[beam_keys[i]] = \
                                val["data"]["beamparams"][0][0][0][0][i][0][0]
        # Extract beam (Lucretia values)
        BunchInterval         = val['data']["beam"][0][0][0][0][0][0][0]
        bunch["x"]            = \
                       np.squeeze(val['data']["beam"][0][0][0][0][1][0][0][0])
        bunch["Q"]            = \
                       np.squeeze(val['data']["beam"][0][0][0][0][1][0][0][1])
        bunch["stop"]         = \
                       np.squeeze(val['data']["beam"][0][0][0][0][1][0][0][2])
        beam["BunchInterval"] = BunchInterval
        beam["Bunch"]         = bunch
        # Add to data
        data["beam"]          = beam
        data["beamparams"]    = beamparams
            # Extract scanvals from matlab data
        scanvals["P1"] = val['scanvals']['P1'][0][0][0][0]
        scanvals["P2"] = val['scanvals']['P2'][0][0][0][0]
        scanvals["V1"] = val['scanvals']['V1'][0][0][0][0]
        scanvals["V2"] = val['scanvals']['V2'][0][0][0][0]
        scanvals["qi"] = val['scanvals']['qi'][0][0][0][0]
        scanvals["dx"] = val['scanvals']['dx'][0][0][0][0]
        scanvals["dy"] = val['scanvals']['dy'][0][0][0][0]
        return data, scanvals