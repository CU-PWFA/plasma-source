"""
Module for computing EOS-BPM signal and performing analysis
"""
import matplotlib.pyplot as plt 
plt.style.use('huntstone')
import numpy as np
from scipy.constants import c
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import sys
#sys.path.insert(0,"home/keenan/eos_bpm/khuntstone/current_profiles/")
import currents as cp
from crystal import crystal
from laser import laser
import phase_retard as pr
from plotting import makefig
import thz


def get_signal(ind, ctype, d, y0, tp, angle, r0, method, th = 0, plot = False):
	"""
	Function to run a full eos simulation and compute a signal for a given 
	current profile. 

	Parameters:
	-----------
	ind : int,
		  The index of the current profile
	ctype  : str, 
			 The type of crystal to use, either 'GaP' or 'ZnTe'
	d      : float
			 The crystal thickness
	y0     : float
			 The central wavelength of the probe in m
	tp     : float
			 The FWHM of the probe in seconds
	angle  : float
			 The crossing angle of the probe in radians
	r0     : float
			 The beamline-crystal distance in m
	th     : float, optional
			 The angle of the half-wave plate for near balanced detectors (rad)
	method : str,
			 The detection method, either 'cross', 'bal', or 'near' for 
			 crossed polarizers, balanced detectors, or near for near balanced 
			 detectors. 
	
	
	
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
	# Initialize crystal and parameters
	cry = crystal(ctype)
	nslice = 100;
	j     = np.arange(1, nslice, 1)
	dz    = d / nslice
	d_arr = (j - 0.5) * dz

	# Initialize probe
	dy = 27e-9;
	probe = laser({'y0' : y0, 'dy' : dy, 'tp' : tp});

	I, ti, p2p = cp.get_current(ind)
	# Compute electric field
	E, te = cp.get_E(I, ti, r0)
	# Make symmetric
	N      = 1000
	fE     = interp1d(te, E)
	t_use  = min([abs(te[0]), te[-1]])
	te_int = np.linspace(-t_use, t_use, N)
	E_int  = np.flip(fE(te_int))

	# Compute Signal
	FEr, f = thz.raw_field(E_int, te_int);
	Ec, tt = thz.cry_field(te_int, FEr, f, d, probe, cry, nslice = nslice);
	tau    = np.linspace(0, 3, 1000) * 1e-12
	gamma, t_gamma = pr.phase_retard(Ec, tt*1e-12, d_arr, tau, probe, cry,\
									 'spatial', psi = angle, plot = False)

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
		ind1 = np.argmin(abs(t - 0.3e-12))
		ind2 = np.argmax(signal[ind1:-1])
		wind = ind1 + ind2
		drid = np.argmax(signal)
		p2p  = t[wind] - t[drid]
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


	
