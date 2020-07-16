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
#sys.path.insert(0, "/home/keenan/eos_bpm/python")
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
	
def get_peak(ind, d, th, r0, fpath):
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
	maxG : float
		   The maximum phase of the EOS-BPM setup
	res  : float
	       The response of EOS-BPM to a small offset based on maxG
	"""
	setup = {"ctype"  : "GaP", 
			 "d"      : d, 
			 "y0"     : 800e-9,
			 "tp"     : 30e-15,
			 "angle"  : th, 
			 "r0"     : r0,
			 "method" : "cross",
			 "fpath"  : fpath,
			 "tilt"   : 0,
			 "th"     : 0,
			 "nslice" : 100,
			 "plot"   : False, 
			 }
	I, ti, Idz, sig1, t_sig1, gamma, t_gamma = eos.get_signal(ind, setup)
	maxG = np.nanmax(gamma)
	res  = maxG * np.sin(maxG)
	return maxG, res

def scan1D(r0s, d, th, fpath, N = 3134, N_p = 4):
	"""
	Function to perform a 1D parameter scan of crystal-beamline distance and 
	compute the average peak signal vs. transverse offset
	
	Parameters:
	-----------
	r0s : array_like
		  Array of crystal-beamline distances to scan (m)
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
	G_peaks : array_like
			  1D array of avg peak phase vs r0
	G_std   : array_like
			  1D array of the std of phase vs r0 
	resps   : array_like
	          1D array of the average responses (increase/decrease) of the
	          signal to the offset
	res_std : array_like
	          1D array of the std of resps
	"""
	
	# Create array of indices to scan
	inds = np.arange(0, N, 1)
	# Preallocate for loo[]
	G_peaks  = np.zeros(len(r0s))
	G_std    = np.zeros(len(r0s))
	resps    = np.zeros(len(r0s))
	res_std  = np.zeros(len(r0s))
	start = time.time()
	for i in range(len(r0s)):
		if (i+1)%10 == 0:
		  print(i+1, "of", len(r0s))
		smax      = partial(get_peak, d = d, th = th, \
			                r0 = r0s[i], fpath = fpath)
		pool      = Pool(N_p)
		gmax, res = zip(*pool.map(smax, inds))
		pool.close()
		pool.join()
		G_peaks[i] = np.nanmean(gmax)
		G_std[i]   = np.nanstd(gmax)
		resps[i]   = np.nanmean(res)
		res_std[i] = np.nanstd(res)

	print("Completed in", time.time()-start, "seconds")
	return G_peaks, G_std, resps, res_std

def get_slope(ind, r0, fpath):
	"""
	Function to compute the slope of the BPM signal for drive and witness 
	bunches, all slopes are per micron of offset
	Returns:
	d_slope  : float
			   Slope of the drive beam peak bpm signal from simulation 
	w_slope  : float
			   Slope of the witness beam peak bpm signal from simulation
	cd_slope : float
			   Slope of the drive beam peak bpm signal from theory 
	cw_slope : float
			   Slope of the witness beam peak bpm signal from theory
	"""
	# Set parameters
	setup = {"ctype"  : "GaP", 
			 "d"      : 75e-6, 
			 "y0"     : 800e-9,
			 "tp"     : 30e-15,
			 "angle"  : 15, 
			 "r0"     : r0,
			 "method" : "cross",
			 "fpath"  : fpath,
			 "tilt"   : 0,
			 "th"     : 0,
			 "nslice" : 100,
			 "plot"   : False, 
			 }
	# Set offset array, low-res is fine, must be symmetric
	offsets = np.linspace(-10, 10, 21) * 1e-6
	dt      = 0.4e-12
	# Preallocate arrays
	d_peak = np.zeros(len(offsets))
	w_peak = np.zeros(len(offsets))
	# Compute crystal signal across offsets
	for i in range(len(offsets)):
		setup["r0"] = r0 + offsets[i]
		sim         = eos.get_signal(ind, setup)
		d_peak[i]   = max(sim[3])
		w_ind       = np.argmin(abs(sim[4] - dt))
		w_peak[i]   = max(sim[3][w_ind:-1])
	# Compute bpm signal and slope from simulation
	d_bpm   = d_peak - np.flip(d_peak)
	w_bpm   = w_peak - np.flip(w_peak)
	d_slope = np.nanmean(np.diff(d_bpm) / np.diff(offsets*1e6))
	w_slope = np.nanmean(np.diff(w_bpm) / np.diff(offsets*1e6)) 
	setup["r0"] = r0
	
	# Calculate slope from theory
	sim0  = eos.get_signal(ind, setup)
	w_ind = np.argmin(abs(sim0[6] - dt))
	gd0   = max(sim0[5])
	gw0   = max(sim0[5][w_ind:-1])
	cd_slope = (-np.sin(gd0) / r0) * 1e-6
	cw_slope = (-np.sin(gw0) / r0) * 1e-6

	return d_slope, w_slope, cd_slope, cw_slope

def scan_slope(r0s, fpath, N = 3135, N_p = 4):
	"""
	Function to compute the average slope of the peak bpm signal 
	(both drive and witness) for an array of r0s

	Parameters:
	
	r0s   : array_like
	        Array of crystal beamline distances to scan
	fpath : string
	        Full file path to the current profiles and dz .mat files
	N     : int, optional
	        Number of current profiles to average over, default is 3135
	N_p   : int, optional
	        Number of processers to use, default is 4

	Returns:
	--------
	"""
	# Create array of indices to scan across
	inds     = np.arange(0, N, 1)
	# Preallocate for loops
	d_slopes  = np.zeros(len(r0s))
	w_slopes  = np.zeros(len(r0s))
	cd_slopes = np.zeros(len(r0s))
	cw_slopes = np.zeros(len(r0s))
	start     = time.time()
	# Loop through r0s
	for i in range(len(r0s)):
		print(i+1, "of", len(r0s))
		pool  = Pool(N_p)
		errf  = partial(get_slope, r0 = r0s[i], fpath = fpath)
		ds, ws, cds, cws   = zip(*pool.map(errf, inds))
		pool.close()
		pool.join()
		d_slopes[i] = np.nanmean(ds)
		w_slopes[i] = np.nanmean(ws)
		cd_slopes[i] = np.nanmean(cws)
		cw_slopes[i] = np.nanmean(cds)    
	print("Completed in", time.time() - start, "seconds")
	
	return d_slopes, w_slopes, cd_slopes, cw_slopes
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

def plot_1D(r0s, g_max, g_std, res, res_std):
	"""
	Function to plot the signal vs transverse offset for each r0 in scan_1D. 
	Also computes and plots the overall drop in signal and the linearity of 
	S_peaks vs dx. 
	
	Parameters:
	-----------
	r0s     : array_like
			  Array of crystal beam distances input into scan1D
	g_max   : array_like
			  array of max phase retardation for each r0
	g_std   : array_like
			  Array of phase standard deviation per r0
	res     : array_like
			  array of max response for each r0
	res_std : array_like
			  Array of response standard deviation per r0
	"""

	fig1, ax1 = makefig(xlab = r'$r_0$ [mm]', \
		                ylab = r'Max($\Gamma_0$)')
	ax1.plot(r0s*1e3, g_max, "-k")
	ax1.fill_between(r0s*1e3, g_max - g_std, g_max+g_std, color = "r", \
		             alpha = 0.5)

	fig1, ax1 = makefig(xlab = r'$r_0$ [mm]', \
		                ylab = r'Max($\Gamma_0$ sin$\Gamma_0$)')
	ax1.plot(r0s*1e3, res, "-k")
	ax1.fill_between(r0s*1e3, res- res_std, res+res_std, color = "r", \
		             alpha = 0.5)

	plt.show()

def plot_slopes(r0s, d_slopes, w_slopes, cd_slopes, cw_slopes):
	"""
	Function to plot the simulated and calculated slopes, parameters are
	r0s the input into scan_slope and the returned slopes from scam_slope 
	"""

	fig, ax = makefig(xlab = r'$r_0$ [mm]', ylab = r'Slope [AU/ $\mu$m]')
	ax.plot(r0s*1e3, d_slopes, 'bo', label = "Simulation")
	ax.plot(r0s*1e3, cd_slopes, 'rx', label = "Theory")
	ax.legend()
	ax.set_title("Drive beam")

	fig2, ax2 = makefig(xlab = r'$r_0$ [mm]', ylab = r'Slope [AU/ $\mu$m]')
	ax2.plot(r0s*1e3, w_slopes, 'bo', label = "Simulation")
	ax2.plot(r0s*1e3, cw_slopes, 'rx', label = "Theory")
	ax2.legend()
	ax2.set_title("Witness beam")

	plt.show()