# Module for studying EOS-BPM response to tilted SFQED bunches
# Standard python imports
#import json
import matplotlib.pyplot as plt 
import numpy as np 
from scipy.constants import c, epsilon_0
from scipy.interpolate import interp1d
try:
	plt.style.use("huntstone")
except:
	plt.style.use("default")
#import sys
# Get system paths from config files
#with open("conf.json") as json_conf:
#	CONF = json.load(json_conf)
#sys.path.insert(0, CONF["py_lnx"])
#sys.path.insert(0, CONF["py_wind"])
#sys.path.insert(0, CONF["py_winl"])
# Custom modules
import eo_signal as eos 

# Tilt analysis function
def get_dur(sig, tsig):
	"""
	Gets the duration (non-zero values) of a signal or current profile, 
	useful for the longer SFQED bunch. 

	Parameters:
	-----------
	sig  : array_like
		   Array of the signal/current
	tsig : array_like
		   Temporal array corresponding to sig
	Returns:
	--------
	dur : float
		  Temporal duration of the signal
	"""
	isSig = np.squeeze(np.argwhere(np.squeeze(sig) > 0))
	dur   = tsig[max(isSig)] - tsig[min(isSig)]
	return dur

def getE(I, ti, r0, tilt, cor, interp = False):
	"""
	Function to compute the Electric field of a tilted electron bunch.

	Parameters:
	-----------
	I      : array_like
			 Current profile of the electron beam
	ti     : array_like
			 Time array corresponding to
	r0     : float
			 The crystal beamline distance in m
	tilt   : float
			 The tilt of the bunch in rad 
	cor    : float
			 The "center of rotation" of the bunch in s
	interp : bool, optional
			 Whether or not to interpolate the Electric field (rec. if 
			 running in the EOS-BPM simulation)

	Returns:
	--------
	E  : array_like
		 The electric field from the tilted bunch
	te : array_like
		 Temporal array corresponding to E
	"""

	# Get Charge in each slice of the bunch, qps
	dti   = (ti[1] - ti[0])
	qps   = I * 1e3 * dti
	# Compute longitudinal, z, position of the tilted bunch
	zi    = (ti - cor) * c
	tz    = zi * np.cos(tilt) / c
	dzi   = c * (tz[1] - tz[0])
	# Compute x-offset due to the tilt
	dx = c * ti * np.tan(tilt)
	# Compute Electric field
	Qe   = np.sum(qps)
	lamz = I * 1e3 * dti / (Qe * c) # Distribution (m^-1)
	lamz = lamz / np.sum(lamz * dzi)
	E    = 2 * Qe * lamz / (4 * np.pi * epsilon_0 * (r0 + dx))
	if interp:
		te_int = np.linspace(-3, 3, 1000)*1e-12
		fE = interp1d(tz, E, bounds_error = False, fill_value = 0)
		E_int = fE(te_int)
		return E_int, te_int, dx
	else:
		return E, tz, dx
def gaussianE(Q, sigt, t, r0, tilt):
    """
    Function to compute a Gaussian electric field for use in EOS-BPM sims. 
    
    Parameters:
    -----------
    Q    : float
           Bunch charge in C
    sigt : float
           Bunch rms length in s
    t    : array_like
           Time array (s) along which the field is computed
    r0   : float
           Radial distance at which the field is computed 
           (i.e. crystal separation)
    tilt : float, optional
           Tilt, if any, of the electron beam in rad.
    Returns:
    --------
    E : array_like
        The electric field of the bunch
    """
    # Compute offset due to tilt
    dx = c * t * np.tan(tilt)
    E0 = Q / ((2*np.pi)**(1.5) * epsilon_0 * (r0 + dx) * c * sigt)
    E  = E0 * np.exp(-t**2 / (2 * sigt**2))
    return E
    
def tilt_signal(I, ti, r0, tilt, cor, setup):
	"""
	Parameters:
	-----------
	I     : array_like
			The current profile of the electron beam
	ti    : array_like
			The time profile corresponding to I
	r0    : float
			The crystal beamline distance 
	tilt  : float
			The tilt of the electron beam
	cor   : float
			The "center of rotation" of the bunch in s
	setup : dictionary object,
			The setup dictionary as described in eo_signal.
	"""
	# Compute the electric field as seen from both crystals
	Ep, tp, dm = getE(I, ti, r0, tilt, cor, interp = True)
	Em, tm, dm = getE(I, ti, r0, -tilt, cor, interp = True)
	# Compute EOS signal for each crystal
	simp = eos.E_signal(Ep, tp, setup)
	simm = eos.E_signal(Em, tm, setup)

    # return signal [0] and time [1] for both fields
	return simp[0], simp[1], simm[0], simm[1]

def comp_tilt(sigP, sigM, tsig, r0, trange, method = "cross"):
    """
    Function to compute the tilt of a bunch from the EOS-BPM signal
    Parameters:
    -----------
    sigP : array_like
           EOS-BPM signal from crystal in +x
    sigM : array_like
           EOS-BPM signal from crystal in -x
    tsig : array_like
           Time array corresponding to sigP and sigM
    r0   : float
           Crystal beamline separation in m
    trange: array_like
            Region of interest for computing the tilt [t1, t2]
         
    method : str, optional
             Detector setup, either "cross" for crossed polarizers or "bal" for
             balanced detectors
    Returns:
    --------
    dx   : array_like
           Computed transverse offset of the electron bunch
    tilt : float
           Computed tilt of the electron bunch
    """
    # Apply filter
    ind1 = np.argmin(abs(tsig - trange[0]))
    ind2 = np.argmin(abs(tsig - trange[1]))
    # Compute offset due to tilt
    if method == "cross":
        dx = (np.sqrt(sigP) - np.sqrt(sigM)) / (np.sqrt(sigP) + np.sqrt(sigM))
        dx = r0 * dx
    elif method == "bal":
        dx = r0 * ((sigP - sigM) / (sigP + sigM))
    else:
        print("Method not recognized")
        return
    # Compute tilt
    dx = dx[ind1:ind2]
    ps = np.polyfit(c*tsig[ind1:ind2], dx, 1)
    tilt = np.arctan(ps[0])
    return dx, tilt

# Some otpimization functions

def avg_g0(Q, sigt, r0s, ds, setup, plot = False):
    """
    Computes average phase retardation for a variety of beamline separations
    and crystal thicknesses and given bunch parameters
    
    Parameters:
    -----------
    Q     : float,
            Bunch charge (C)
    sigt  : float
            Bunch rms duration (s)
    r0s   : array_like
            Array of beamline separations (m)
    ds    : array_like
            Array of crystal thicknesses (m)
    setup : dictionary obj
            EOS-BPM setup dictionary as defined in eo_signal.py (omitting r0 
            and d)
    plot  : bool, optional
            Whether or not to plot results. 
    Rerturns:
    ---------
    gap_g0  : array_like
              Array of average phase retardation of a GaP crystal
    znte_g0 : array_like
              Array of average phase retardation of a ZnTe crystal
    """
    # Useful function for computing peak average phase retardation 
    def get_gamma(r0, d, ctype, Q, sigt, setup, N = 8000):
        setup["r0"]    = r0
        setup["d"]     = d
        setup["ctype"] = ctype
        dt = sigt / 10
        tb = np.linspace(-0.5*N*dt, 0.5*N*dt, N)
        setup["tau"] = tb
        E  = gaussianE(Q, sigt, tb, r0, 0)
        return max(eos.E_signal(E, tb, setup)[2])
    # Preallocate for loop
    gap_g0  = np.zeros((len(ds), len(r0s)))
    znte_g0 = np.zeros((len(ds), len(r0s)))
    for i in range(len(ds)):
        print(i+1, "of", len(ds))
        d = ds[i]
        for j in range(len(r0s)):
            r0 = r0s[j]
            gap_g0[i, j] = get_gamma(r0, d, "GaP", Q, sigt, setup)
            znte_g0[i, j] = get_gamma(r0, d, "ZnTe", Q, sigt, setup)