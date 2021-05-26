'''
Module to compute the phase retardation experienced by a probing pulse 
propagating through an EO crystal. Allows for several different probe setups

author : keenan
'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c
from scipy.interpolate import interp1d

from crystal import crystal
from laser import laser
from plotting import makefig
import thz

def phase_retard(Ec, te, d, tau, probe, crystal, psi = 0, interp = False, \
                 n_interp = 10000, t_delay = 0, plot = False):
    """
    Function to compute the phase retardation spatially encoding into a 
    probing laser pulse. 

    Parameters:
    -----------
    Ec       : array_like
               Array of the effective THz pulse in each crystal slice as 
               created in the thz.py module
    te       : array_like
               1D array corresponding to the time axis of Ec (s)
    d        : float
               The EO crystal thickness (m)
    tau      : array_like
               Time array corresponding to the probe laser (effectively the
               probe "timing delay" due to the sweeping angle, s)
    probe    : object
               Instance of the laser class used in the EOS setup
    crystal  : object
               Instance of the crystal class used in the EOS setup
    psi      : float, optional
               Probe crossing angle (deg.), default = 0 (no spatial encoding)
    interp   : bool, optional
               Set True for interpolating gamma before returning
    n_interp : int, optional
               Length of interpolated array defauult = 10000
    t_delay  : float, optional
               Optional time delay of the laser probe (s), default = 0
    plot     : bool, optional
               Set to true to plot the computed phase retardation

    Returns:
    --------
    gamma   : array_like
              The phase retardation (rad)
    t_gamma : array_like
              The time array corresponding to gamma (s)
    """
    # Convert psi to rads
    psi = psi * np.pi / 180
    # Add shift
    tau = tau + t_delay
    # Create depth array
    nslice = np.shape(Ec)[1] + 1;
    j     = np.arange(1, nslice, 1)
    dz    = d / nslice
    d_arr = (j - 0.5) * dz
    # Compute Amplitude of gamma 
    n0 = crystal.indref(np.array([probe.y0]))[0]
    amp = 2 * np.pi * n0**3 * dz/ (probe.y0)
    # Get effective probe velocity
    Lchar, v_g_opt = probe.laser_l_char(crystal)
    vg_eff = v_g_opt * np.cos(psi)
    # Preallocate for loop
    gamma = np.zeros(len(tau))
    # Loop and sum gamma
    tau_use = tau + t_delay
    for j in range(len(d_arr)):
        fEc = interp1d(te, np.real(Ec[:, j]), bounds_error = False, \
                       fill_value = 0)
        # Ec is effectively wrapped, reset the probe timing if it is past 
        # the pulse domain
        #ind = np.argwhere(tau_use > te[-1])
        #dtau = tau_use[ind] - te[-1]
        #tau_use[ind]  = te[0] + dtau
        t_probe = d_arr[j] / vg_eff
        #if t_probe > te[-1]:
        #    diff = t_probe - te[-1]
        #    t_probe = te[0] + diff
        t_interp = t_probe + tau_use
        E_eff = fEc(t_interp)
        gamma += E_eff
    # Add amplitude
    gamma = amp * gamma
    # Artificially center gamma at t = 0
    tau = tau - tau[np.argmax(gamma)]

    if interp:
        fgamma    = interp1d(tau, gamma)
        tau_int   = np.linspace(tau[0], tau[-1], n_interp)
        gamma_int = fgamma(tau_int)
        gamma     = gamma_int
        tau       = tau_int
    if plot:
        fig, ax = makefig(xlab = "t [ps]", ylab = r'$\Gamma$ [rad]')
        ax.plot(tau * 1e15, gamma)
        plt.show() 
    return gamma, tau

def eosd(E, te, setup):
    """
    Currently in testing. Function to compute the phase retardation
    imparted on a chirped pulse (EOSD). Treats the pulse as a 'fine'
    frequency comb. 
    """
    ctype    = setup["ctype"]   # crystal type
    d        = setup["d"]       # crystal thickness
    y0       = setup["y0"]      # probe central wavelength
    tp       = setup["tp"]      # probe FWHM duration
    angle    = setup["angle"]   # probe crossing angle
    r0       = setup["r0"]      # crystal beamline distance
    method   = setup["method"]  # detection method
    nslice   = setup["nslice"]  # number of crystal slices
    tau      = setup["tau"]     # Probing time array
    t_delay  = setup["t_delay"] # Artificial probe delay
    tchirp   = setup["tc"]      # Chirped pulse FWHM
    plot     = setup["plot"]    # Plotting preference (True/False)
    
    # Initialize crystal and parameters
    cry = crystal(ctype)
    # Slice the crystal
    j     = np.arange(1, nslice, 1)
    dz    = d / nslice
    d_arr = (j - 0.5) * dz
    # Initialize probe
    dy = 27e-9;
    probe = laser({'y0' : y0, 'dy' : dy, 'tp' : tp})
    probe.chirp(tchirp)
    tlim = np.sqrt(tp*tchirp)
    # Compute Effective THz pulse
    FEr, f = thz.raw_field(E, te);
    Ec, tt = thz.cry_field(te, FEr, f, d, probe, cry, nslice = nslice)
    # Convert angle to rads
    psi = angle * np.pi / 180
    # preallocate for loop
    gamma = np.zeros(len(tau))
    upd   = int(len(tau)/10)
    for i in range(len(tau)):
        #if (i+1)%upd==0:
        #    print(i+1, "of", len(tau))
        # Get instantaneous frequency/wavelength and create a mini probe
        wi = probe.get_inst_w(tau[i])
        yi = 2*np.pi*c/wi
        mini_probe = laser({"y0" : yi, "dy" : 0, "tp" : tlim})
        n0   = cry.indref(np.array([mini_probe.y0]))[0]
        dz   = d_arr[1]-d_arr[0]
        amp1 = 2 * np.pi * n0**3 * dz/ (mini_probe.y0)
        amp2 = probe.get_inst_amp(tau[i])
        # Get effective probe velocity
        Lchar, vg_opt = mini_probe.laser_l_char(cry)
        tau_use = tau[i] + t_delay
        for j in range(len(d_arr)):
            fEc = interp1d(tt*1e-12, np.real(Ec[:, j]), bounds_error = False, \
                           fill_value = 0)
            t_probe = d_arr[j] / vg_opt
            t_interp = t_probe + tau_use
            E_eff = fEc(t_interp)
            gamma[i] += E_eff
        # Factor in amplitudes
        gamma[i] = amp1*amp2*gamma[i]
    # Artificially center gamma at t = 0
    tau = tau - tau[np.argmax(gamma)]
    return gamma, tau

def broad_signal(E, te, setup):
    """
    Currentyly in testing. Function to compute phase retardation 
    imparted on a short pulse. Takes into account the """
    # Extract variables
    ctype   = setup["ctype"]   # crystal type
    d       = setup["d"]       # crystal thickness
    y0      = setup["y0"]      # probe central wavelength
    tp      = setup["tp"]      # probe FWHM duration
    angle   = setup["angle"]   # probe crossing angle
    r0      = setup["r0"]      # crystal beamline distance
    method  = setup["method"]  # detection method
    nslice  = setup["nslice"]  # number of crystal slices
    tau     = setup["tau"]     # Probing time array
    t_delay = setup["t_delay"] # Artificial probe delay
    plot    = setup["plot"]    # Plotting preference (True/False)
    
    # Initialize crystal and parameters
    cry = crystal(ctype)
    # Initialize probe
    dy = 27e-9;
    probe = laser({'y0' : y0, 'dy' : dy, 'tp' : tp})
    # Compute Effective THz pulse
    FEr, f = thz.raw_field(E, te);
    Ec, tt = thz.cry_field(te, FEr, f, d, probe, cry, nslice = nslice)
    # Convert angle to rads
    psi = angle * np.pi / 180
    # Add shift
    tau = tau + t_delay
    # Create depth array
    nslice = np.shape(Ec)[1] + 1;
    j     = np.arange(1, nslice, 1)
    dz    = d / nslice
    d_arr = (j - 0.5) * dz
    # Compute Amplitude of gamma 
    n0 = cry.indref(np.array([probe.y0]))[0]
    amp = n0**3 * dz/ (probe.y0)
    # Get effective probe velocity
    Lchar, vg_opt = probe.laser_l_char(cry)
    Lchar2 = Lchar*Lchar
    inv_L2 = 1/Lchar2
    vg_eff = vg_opt * np.cos(psi)
    # Preallocate for loop
    gamma = np.zeros(len(tau))
    # Loop and sum gamma
    tau_use = tau + t_delay
    upd = int(len(tau)/10)
    for i in range(len(tau)):
        if (i+1)%upd==0:
            print(i+1, "of", len(tau))
        for j in range(len(d_arr)):
            sigj = probe.sigp * np.sqrt(1 + (d_arr[j] / Lchar)**2)
            f    = interp1d(tt*1e-12, np.real(Ec[:, j]),\
                            bounds_error = False, fill_value = 0)

            t_interp = tt*1e-12 + (d_arr[j] / vg_opt)
            E_eff = f(t_interp)
            t_delay = tau[i]
            integrand = E_eff * \
                        np.exp(-(tt*1e-12 - t_delay)**2 / (2 * sigj**2)) / (sigj)
            gamma[i] += np.trapz(integrand, tt*1e-12)
    # Artificially center peak gamma at t = 0
    tau = tau - tau[np.argmax(gamma)]
    return amp*gamma, tau

