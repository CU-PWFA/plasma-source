# Module for studying EOS-BPM response to tilted SFQED bunches
# Standard python imports
import json
import matplotlib.pyplot as plt 
import numpy as np 
from scipy.constants import c, epsilon_0
from scipy.interpolate import interp1d
try:
	plt.style.use("huntstone")
except:
	plt.style.use("default")
import sys
# Get system paths from config files
with open("conf.json") as json_conf:
	CONF = json.load(json_conf)
sys.path.insert(0, CONF["py_lnx"])
sys.path.insert(0, CONF["py_wind"])
sys.path.insert(0, CONF["py_winl"])
# Custom modules
import eo_signal as eos 
from plotting import makefig 

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
	rbeam = zi * np.sin(tilt)
	# Compute Electric field
	Qe   = np.sum(qps)
	lamz = I * 1e3 * dti / (Qe * c) # Distribution (m^-1)
	lamz = lamz / np.sum(lamz * dzi)
	E    = 2 * Qe * lamz / (4 * np.pi * epsilon_0 * (r0 + rbeam))
	if interp:
		te_int = np.linspace(-3, 3, 1000)*1e-12
		fE = interp1d(tz, E, bounds_error = False, fill_value = 0)
		E_int = fE(te_int)
		return E_int, te_int, rbeam
	else:
		return E, tz, rbeam

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

def tilt_delta(Sp, Sm, tsig, r0):
	"""
	Function calcuate the tilt of a bunch from the EOS-BPM signal

	Parameters:
	-----------
	Sp   : array_like
		   Array of positive (closer) crystal signal
	Sm   : array_like
		   Array of negative (further) crystal signal
	tsig : array_like
		   Temporal array corresponding to Sp and Sm
	r0   : float
		   The crystal beamline distance
	Returns:
	--------
	theta : float
			The computed beam tilt
	dx    : array_like
			The computed transverse offset of the beam
	plot  : bool, optional
	        Whether or not to plot the signal
	"""
	Rnum   = np.sqrt(Sp) - np.sqrt(Sm)
	Rden   = np.sqrt(Sp) + np.sqrt(Sm)
	R12    = Rnum / Rden
	deltas = - R12 / (1 + R12)
	dx     = deltas * r0
	return dx[~np.isnan(dx)], tsig[~np.isnan(dx)]


