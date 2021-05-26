"""
Module for computing EOS-BPM signal and performing analysis
"""
import matplotlib.pyplot as plt 
#try:
#    plt.style.use('huntstone')
#except:
#    plt.style.use('dark_background')
import numpy as np
from scipy.constants import c
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import currents as cp
from crystal import crystal
from laser import laser
import phase_retard as pr
from plotting import makefig
import thz


def get_signal(ind, setup):
    """
    Function to run a full eos simulation and compute a signal for a given 
    current profile (or E-field) and detector setup.  

    Parameters:
    -----------
    ind : int,
          The index of the current profile
    setup : dictionary object
            dictionary containing the following keys and values
            "ctype"  : str, the crystal kind being used
            "d"      : float, the crystal thickness in m
            "y0"     : float, the probe central wavelength in m
            "tp"     : float, the probe FWHM duration in s
            "angle"  : float, the probe crossing angle in deg.
            "r0"     : float, the crystal beamline offset in m
            "method" : str, the detection scheme being used
            "fpath"  : str, the full path to the current and dz .mat files
            "tilt"   : foat, tilt angle of beam
            "th"     : float, the waveplate angle in near-crossed pol. setup
            "nslice" : int, the number of crystal slices to use (rec 100)
            "tau"    : array_like, array of probe delay time
            "plot"   : bool, whether or not to plot the resultant signal
    
    Returns:
    --------
    I     : array_like
            The input current profile (kA)
    ti    : array_like
            The time array corresponding to I
    sig   : array_like
            The EOS signal
    t_sig : array_like
            Time array corresponding to sig
    """
    # Extract variables (setup MUST have the following)
    ctype  = setup["ctype"]
    d      = setup["d"]
    y0     = setup["y0"]
    tp     = setup["tp"]
    angle  = setup["angle"]
    r0     = setup["r0"]
    method = setup["method"]
    fpath  = setup["fpath"]
    th     = setup["th"]
    nslice = setup["nslice"]
    tau    = setup["tau"]
    plot   = setup["plot"]
    # Initialize crystal and parameters
    cry = crystal(ctype)

    # Initialize probe
    dy = 27e-9;
    probe = laser({'y0' : y0, 'dy' : dy, 'tp' : tp});

    I, ti, p2p = cp.get_current(ind, fpath)
    # Compute electric field
    E, te = cp.get_E(I, ti, r0)
    # Make symmetric
    N      = 1000
    fE     = interp1d(te, E)
    t_use  = min([abs(te[0]), te[-1]])
    te_int = np.linspace(-t_use, t_use, N)
    E_int  = np.flip(fE(te_int))

    # Compute Effective THz pulse
    FEr, f = thz.raw_field(E_int, te_int);
    Ec, tt = thz.cry_field(te_int, FEr, f, d, probe, cry, nslice = nslice)
    # Compute probe timing window (camera res. of 3.5 microns assumed)
    #dtau   = (3.5e-6 / c) * np.tan(angle * np.pi / 180)
    #tau    = np.arange(-500e-15, 1150e-15, dtau)
    
    if method == "eosd":
        gamma, t_gamma = pr.eosd(Ec, te, setup)

    # Compute signal
    if method == "cross":
        sig = np.sin(gamma / 2)**2
    elif method == "bal":
        sig = np.sin(gamma)
    elif method == "near":
        sig = 1 - np.cos(gamma + 4 * th)
    # Re-center time
    t_sig = t_gamma - t_gamma[np.argmax(gamma)]

    if plot:
        plot_signal(I, sig, ti, t_sig)
    return I, ti, p2p, sig, t_sig, gamma, t_gamma
def E_signal(E, te, setup):
    """
    Function to compute spatially encoded EOS signal from an electric field.
    
    Parameters:
    -----------
    E      : array_like,
             The input electric field
    te     : array_like,
             Time array corresponding to E (s)
    setup : dictionary object
            dictionary containing the following keys and values
            "ctype"  : str, the crystal kind being used
            "d"      : float, the crystal thickness in m
            "y0"     : float, the probe central wavelength in m
            "tp"     : float, the probe FWHM duration in s
            "angle"  : float, the probe crossing angle in deg.
            "r0"     : float, the crystal beamline offset in m
            "method" : str, the detection scheme being used
            "fpath"  : str, the full path to the current and dz .mat files
            "tilt"   : foat, tilt angle of beam
            "th"     : float, the waveplate angle in near-crossed pol. setup
            "nslice" : int, the number of crystal slices to use (rec 100)
            "plot"   : bool, whether or not to plot the resultant signal
    Returns:
    --------
    """
    # Extract parameters
    ctype  = setup["ctype"]
    d      = setup["d"]
    y0     = setup["y0"]
    tp     = setup["tp"]
    angle  = setup["angle"]
    method = setup["method"]
    th     = setup["th"]
    nslice = setup["nslice"]
    tau    = setup["tau"]
    # Initialize crystal and parameters
    cry = crystal(ctype)
    nslice = 100

    # Initialize probe
    dy = 27e-9
    probe = laser({'y0' : y0, 'dy' : dy, 'tp' : tp})

    # Make symmetric
    N      = len(te)
    fE     = interp1d(te, E)
    t_use  = min([abs(te[0]), te[-1]])
    te_int = np.linspace(-t_use, t_use, N)
    E_int  = np.flip(fE(te_int))

    # Compute Signal
    FEr, f = thz.raw_field(E_int, te_int);
    Ec, tt = thz.cry_field(te_int, FEr, f, d, probe, cry, nslice = nslice)
    
    if method != "eotd":
        # Not using temporal decoding, can just grab the phase
        gamma, t_gamma = pr.phase_retard(Ec, tt*1e-12, d, tau, probe, cry,\
                                         psi = angle)
        # Compute signal
        if method == "cross":
            sig = np.sin(gamma / 2)**2
        elif method == "bal":
            sig = np.sin(gamma)
        elif method == "near":
            sig = 1 - np.cos(gamma + 4 * th)
        # Re-center time
        t_sig = t_gamma - t_gamma[np.argmax(gamma)]
    
    elif method == "eotd":
        # Temporal decoding requires some extra steps
        an = setup["analyzer"] # Method for probe beam 
        I_gate = setup["I_gate"] # Gate pulse intensity (AU)
        # Interpolating function for I_gate
        f_gate = interp1d(t_las, I_gate, bounds_error = False, fill_value = 0)

    return sig, t_sig, gamma, t_gamma


def plot_signal(I, sig, ti, t_sig, save = False, sname = ""):
    """
    Function to plot the input current profile and corresponding EOS signal
    Parameters:
    -----------
    I     : array_like
            Input current profile
    sig   : array_like
            EOS signal
    ti    : array_like
            Time corresponding to I
    t_sig : array_like
            Time corresponding to sign
    save  : bool, optional
            Whether or not to save the figure
    sname : str, optional 
            the name to save the figure
    """
    fig, ax1 = makefig(5, 3.75)
    # Input current
    ax1.spines['left'].set_color('r')
    ax1.tick_params(axis = 'y', color = 'r')
    ax1.tick_params(axis = 'y', labelcolor = 'r')
    ax1.yaxis.label.set_color('r')
    ax1.plot(ti * 1e12, I, '-r')
    ax1.set_ylabel('I [kA]', color = 'r')
    ax1.set_xlabel('t [ps]')
    ax1.set_xlim([min(t_sig*1e12), -1*min(t_sig*1e12)])

    # Signal
    ax2 = ax1.twinx()
    ax2.spines['right'].set_color('b')
    ax2.tick_params(axis = 'y', color = 'b')
    ax2.tick_params(axis = 'y', labelcolor = 'b')
    ax2.yaxis.label.set_color('b')
    ax2.plot(t_sig * 1e12, sig, '--b')
    ax2.set_ylabel('Signal [AU]')
    plt.show()

def peak2peak(signal, t, dt = 0.25e-12):
    """
    Function to compute (if possible) the peak to peak measurement of the 
    signal. 

    Parameters:
    -----------
    sig   : array_like
            The EOS signal
    t_sig : array_like
            The time corresponding to sig

    Returns:
    --------
    p2p : float
          The peak to peak measurement (s), returns nan if noise is strong
    """

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

    try:
        p2p = t[peaks[1]] - t[peaks[0]]
    except:
        #ind1 = np.argmin(abs(t - 0.3e-12))
        #ind2 = np.argmax(signal[ind1:-1])
        #wind = ind1 + ind2
        #drid = np.argmax(signal)
        #p2p  = t[wind] - t[drid]
        p2p = np.nan
    return p2p

def get_wit_stn(sig, t_sig, noise_t1, noise_t2, height, width):
    """
    Function to compute the witness signal to noise ratio.

    Parameters:
    -----------
    sig : array_like
          EOS signal
    t_sig : array_like
            Time corresponding to sig
    noise_t1 : float
               Index of the initial time at which to consider noise
    noise_t2 : float
               Index of the final time at which to consider noise
    height   : float
               Cutoff for finding the witness peak, used in find_peaks
    width    : int
               Integer number of time steps between drive and witness signal
               used for find_peaks
    """

    # Compute noise
    ind1  = np.argmin(abs(t_sig - noise_t1))
    ind2  = np.argmin(abs(t_sig - noise_t2))
    noise = np.mean(sig[ind1:ind2])

    # Compute witness peak
    inds, sigs = find_peaks(sig, height = height, distance = width)
    try:
        wit_sig = sigs['peak_heights'][1]
        return wit_sig, noise 
    except:
        return np.nan, noise

def plot_shot(I, ti, sigA, sigB, tsig, delta, dxd, dxw, angle, r0):
    """
    Function to plot the single shot EOS-BPM signal

    Parameters:
    -----------
    I : array_like
        The input current profile
    ti : array_like
         The time profile of I
    sigA : array_like
           Crystal A signal
    sigB : array_like
           Crystal B signal
    tsig : array_like
           Time array corresponding to sigA/sigB
    dx_comp : arrray_like
              Computed bunch offsets
    dxd     : float
              The actual drive offset
    dxw     : float
              The actual witness offset
    angle   : float
              Probe crossing angle deg. 
    r0      : flaot
              Crystal beamline distance
    """
    S_bpm = sigA - sigB
    # EOS-BPM example shot figure
    fig = plt.figure(figsize = (8, 6), dpi = 200)
    ax1 = fig.add_subplot(211)
    # Current axis
    ax1.tick_params(labelsize = 'large')
    ax1.set_ylabel("I [kA]")
    ax1.spines['left'].set_color('r')
    ax1.tick_params(axis = 'y', color = 'r', labelcolor = 'r')
    ax1.yaxis.label.set_color('r')
    ax1.set_xlim([-0.3, 0.7])
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top') 
    ax1.set_xlabel('t [ps]')
    ax1.set_xlim([-0.3, 0.7])
    # Crystal Signal Axis
    ax2 = ax1.twinx()
    ax2.spines['right'].set_color('b')
    ax2.spines['left'].set_color('r')
    ax2.tick_params(axis = 'y', color = "b", labelcolor = "b")
    ax2.yaxis.label.set_color("b")
    ax2.set_ylabel(r'$S_+ + S_-$ [AU]')
    ax2.set_xlim([-0.3, 0.7])
    # Signal Difference axis
    ax3 = fig.add_subplot(212)# sharex=ax1)
    xlas = (tsig - tsig[0]) * 1e-12 * c * 1e3 / np.tan(angle * np.pi / 180)
    x0 = (-0.3 - tsig[0]) * 1e-12 * c * 1e3 / (np.tan(angle * np.pi / 180))
    x1 = (0.7-tsig[0]) * 1e-12 * c * 1e3 / (np.tan(angle * np.pi / 180))
    ax3.set_xlim([x0, x1])
    ax3.set_xlabel(r'$x_{las}$ [mm]')
    ax3.set_ylabel(r'$S_{+} - S_{-}$ [AU]')
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.tick_params(labelsize = 'large')
    #ax3.spines['right'].set_color('g')
    ax3.tick_params(axis = 'y', color = 'g', labelcolor = 'g')
    ax3.yaxis.label.set_color('g')
    # Transverse offset axis
    ax4 = ax3.twinx()
    ax4.spines["left"].set_color("y")
    ax4.spines["right"].set_color("g")
    ax4.yaxis.tick_left()
    ax4.tick_params(axis = 'y', color = "y", labelcolor="y")
    ax3.yaxis.tick_right()
    ax4.yaxis.label.set_color("y")
    ax4.set_ylabel(r'$\Delta$ x [$\mu$m]')
    ax4.yaxis.tick_left()
    ax4.yaxis.set_label_position("left")
    ax4.tick_params(labelsize = "large")
    plt.subplots_adjust(hspace=0.0)
    # Plotting
    ax1.plot(ti*1e12, I, '-r', linewidth = 1)
    ax2.plot(tsig, sigA + sigB, '-b', linewidth = 1, label = r'$S_+ + S_-$')
    ax3.plot(xlas, S_bpm, '-g')
    ax4.plot(xlas, delta * r0 * 1e6, "-y")
    ax4.axhline(y = dxd*1e6, color = "k", linestyle = '--')
    ax4.axhline(y = dxw*1e6, color = "k", linestyle = '--')
    maxoff = max(abs(np.array([dxd*1e6, dxw*1e6])))
    y1     = -maxoff * 1.4
    y2     = maxoff * 1.4
    ax4.set_ylim([y1, y2])
    plt.show()

    
