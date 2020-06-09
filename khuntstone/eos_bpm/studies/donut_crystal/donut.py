#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studing the signal from a donut-shaped EOS crystal. 

@author: keenan
"""
from cycler import cycler;
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c, epsilon_0;
eps0 = epsilon_0;
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit as cf;
import sys;
sys.path.insert(0, "../../python/");
from crystal import crystal;
from ebeam import ebeam;
from laser import laser;
import phase_retard as pr;
from plotting import makefig;
import thz;


def get_signal(drive, wit, probe, psi, ctype, d, x0, x, y0, y, td, \
				  nx = 1000, ny = 1000):
	"""
	Function to calculate the 2D EOS-BPM signal of a donut crystal.
	Unless otherwise stated inputs are in SI MKS units.
	Parameters:
	-----------
	drive : object (see python/ebeam)
			Instance of the ebeam class with drive bunch paramters
	wit   : object 
			Instance of the ebeam class with witness bunch parameters
	probe : object (see python/laser)
			Instance of the laser class with probe pulse parameters
	psi   : float
			The sweeping angle of the probe in degrees
	ctype : str
			Type of crystal being used, 'GaP' or 'ZnTe'
	d     : float
			The thickness of the crystal
	x0    : array_like
			The crystal beamline horizontal separation (drive, witness)
	x     : float
			The width of the crystal
	y0    : array_like
			The crystal beamline vertical offset (drive, witness)
	y     : float
			The height of the crystal
	td    : float
			Temporal delay of the probe pulse
	nx    : int
			The number of points to calculate signal for in x, default = 100
	ny    : int, optional
			The number of points to calcualte signal for in y, default = 100
	Returns:
	--------
	gamma_drive : array_like
				  2D array of phase retardation due to the drive bunch
	gamma_wit   : array_like
				  2D array of phase retardation due to the witness bunch
	gamma       : array_like
				  2D array of overall phase retardation.
	"""
	# Initialize crystal
	print("Initalizing for loop...");
	cry   = crystal(ctype);

	# Create spatial and temporal arrays
	y_arr = np.linspace(-y/2, y/2, ny);
	x_arr = np.linspace(-x/2, x/2, nx);
	
	tau   = np.zeros((len(x_arr), len(y_arr)));
	for i in range(len(x_arr)):
		for j in range(len(y_arr)):
			r = np.sqrt(x_arr[i]**2 + y_arr[j]**2);
			tau[i, j] = (r / c) * np.tan(psi * np.pi / 180);

	Et_drive            = np.exp(-(drive.t**2) / (2 * drive.sigt**2));
	FErt_drive, f_drive = thz.raw_field(Et_drive, drive.t);
	Ect_drive, tt_drive = thz.cry_field(drive.t, FErt_drive, f_drive, d, \
										probe, cry);

	Et_wit          = np.exp(-(wit.t + wit.del_t)**2 / (2 * wit.sigt**2));
	FErt_wit, f_wit = thz.raw_field(Et_wit, wit.t);
	Ect_wit, tt_wit = thz.cry_field(wit.t, FErt_wit, f_wit, d, probe, cry);

	# Since tau is 2D best to create an interpolation function for speed
	tau_int         = np.linspace(np.amin(tau), np.amax(tau), 1000);
	drive_int, dm   = pr.phase_retard(Ect_drive, tt_drive * 1e-12, d, \
		              tau_int + td, probe, cry, 'spatial', psi = psi);
	wit_int, dm     = pr.phase_retard(Ect_wit, tt_wit * 1e-12, d, \
		              tau_int + td, probe, cry, 'spatial', psi = psi);


	f_drive         = interp1d(tau_int, drive_int, bounds_error = True, \
		                       fill_value = 0);
	f_wit           = interp1d(tau_int, wit_int, bounds_error = True, \
		                       fill_value = 0);

	gamma_drive = np.zeros((len(x_arr), len(y_arr)));
	gamma_wit   = np.zeros((len(x_arr), len(y_arr)));

	print("Looping...");
	for i in range(len(x_arr)):
	    x_in = [x_arr[i] + x0[0], x_arr[i] + x0[1]];
	    for j in range(len(y_arr)):
	        y_in = [y_arr[j] + y0[0], y_arr[j] + y0[1]];
	        r_drive  = np.sqrt(x_in[0]**2 + y_in[0]**2);
	        r_wit    = np.sqrt(x_in[1]**2 + y_in[1]**2); 
	        E0_drive = drive.Q / (2 * np.pi * eps0 * r_drive * \
	                    np.sqrt(2 * np.pi) * c * drive.sigt);
	        E0_wit   = wit.Q / (2 * np.pi * eps0 * r_wit * \
	                  np.sqrt(2 * np.pi) * c * wit.sigt);
	        gamma_drive[i, j] = f_drive(tau[i, j]) * E0_drive;
	        gamma_wit[i, j]   = f_wit(tau[i, j]) * E0_wit;

	gamma = gamma_drive + gamma_wit;
	return gamma_drive, gamma_wit, gamma;

def radial_lineout(signal, x_arr, y_arr, dir, psi, title = '', \
	               plot = False, save = False, sname = ''):
	"""
	Function to compute a radial lineout of the signal. Since there should 
	always be a cardinal direction in which the effects of transverse offsets
	are

	Paramters:
	----------
	signal : array_like
	         an n x m array of the 2D eos signal with dimensions corresponding 
	         to x and y respectively
	x_arr  : array_like
	         The array of horizontal positions corresponding to signal
	y_arr  : array_like:
	         The array of vertical positions corresponding to signal 
	dir    : array_like (strings)
	         Array containing the information on the direction of the lineou
	         first entry is coordinate ('x' or 'y') second entry is direction 
	         ('+' or '-');
	psi    : float
	         The sweeping angle of the probe pulse
	title  : str, optional
	         The title of the plot
	plot   : bool, optional
	         Whether or not to plot the lineout
	save   : bool, otpional
	         Whether or not to save the figure if plotting
	sname  : str, optional
	         The name of the figure if saving

	Returns:
	--------
	lineout : array_like
	          The radial lineout of the signal
	t_plot  : array_like
	          The temporal array corresponding to the lineout
	"""

	# Get indices based on dir
	if dir[0] == 'x':
		if dir[1] == '+':
			x_use   = x_arr[x_arr >= 0]
			lineout = signal[np.argwhere(x_arr >= 0), np.argmin(abs(y_arr))];
			t_plot  = (x_use / c) * np.tan(psi * np.pi / 180) * 1e15;
		elif dir[1] == '-':
			x_use   = x_arr[x_arr <= 0];
			lineout = signal[np.argwhere(x_arr <= 0), np.argmin(abs(y_arr))];
			t_plot  = (x_use / c) * np.tan(psi * np.pi / 180) * 1e15;
	elif dir[0] == 'y':
		if dir[1] == '+':
			y_use   = y_arr[y_arr >= 0];
			lineout = signal[np.argmin(abs(x_arr)), np.argwhere(y_arr >=0)];
			t_plot  = (y_use / c) * np.tan(psi * np.pi / 180) * 1e15;
		elif dir[1] == '-':
			y_use   = y_arr[y_arr <= 0];
			lineout = signal[np.argmin(abs(x_arr)), np.argwhere(y_arr <=0)];
			t_plot  = (y_use / c) * np.tan(psi * np.pi / 180) * 1e15;
	# center t on 0
	t_plot = t_plot - t_plot[int(len(t_plot) / 2)];
	if plot:
		max_ind = np.argmax(lineout);
		tmin    = t_plot[max_ind] - 1000;
		tmax    = t_plot[max_ind] + 1000;
		fig, ax = makefig(x = 8, y = 6, xlab = 't [fs]', ylab = 'Signal [AU]', \
			              title = title)
		ax.plot(t_plot, lineout / max(lineout));
		ax.set_xlim([tmin, tmax])
		if save:
			fig.savefig(sname);
		plt.show();
	return np.squeeze(lineout), np.squeeze(t_plot);

def fit_lineout(lineout, theta, plot = False, title = ''):
	"""
	Function to perform a sinusoidal fit to an azimuthal lineout
	Parameters:
	-----------
	lineout : array_like
	          The azimuthal lineout
	theta   : array_like
	          The angular array corresponding to lineout
	plot    : bool, optional
	          Whether or not to plot the lineout
	title   : str, optional
	          The title of the plot
	Returns :
	---------
	fit_line : array_like
	           The sinusoidal fit line
	popt     : array_like
			   Array of curve fit parameters (amplitude, phase, offset)

	"""
	# Create function for fitting
	def sin_func(x, A, phi, B):
		return A * np.sin(x + phi) + B;

	# Perform the fit
	popt, pcov = cf(sin_func, theta, lineout);
	# Create the best fit line
	fit_line = sin_func(theta, *popt);

	# Plot if requested
	if plot:
		fig, ax = makefig(x = 8, y = 6, xlab = r'$\frac{\theta}{\pi}$', \
			              ylab = 'Signal [AU]')
		ax.plot(theta, lineout, label = 'Lineout');
		# Create label for fit line
		A   = str(np.round(popt[0], 2));
		phi = str(np.round(popt[1] / np.pi, 2)) + r'$\pi$';
		B   = str(np.round(popt[2])); 
		if np.round(popt[1], 2) == 0:
			fit_label = A + "sin(" + r'$\theta$' + ") + " + B;
		else:
			fit_label = A + "sin(" + r'$\theta$' + phi + ") + " + B;
		ax.plot(theta, fit_line, '--r', label = fit_label);
		ax.legend();

	return fit_line, popt
def azi_lineout(signal, x_arr, y_arr, off = 50, \
	             plot = False, title = '', \
	             save = False, snames = ['', '']):
	"""
	Function to perform an azimuthal lineout on the 2D DONUT-EOS signal

	Parameters:
	-----------
	signal : array_like
	         The 2D (x, y) EOS signal
	x_arr  : array_like
	         Array of horizontal positions corresponding to signal
	y_arr  : array_like
	         Array of vertical positions corresponding to signal
	off    : int, optional

	plot   : bool, optional
	         Whether or not to plot the 2D signal in polar coordinates and the
	         azimuthal lineouts
	title  : str, optional
	         titles for the 2D polar plot
	save   : bool, optional
	         Whether or not to save the figure if plotting
	sname  : str, optional
	         The name of the figures if saving

	Returns:
	--------
	sig_polar : array_like
	            The 2D array of EOS signal in polar coordinates
	azi_drive : array_like
	            The azimuthal lineout of the drive beam signal
	azi_wit   : 
	"""

	# Create polar coordinate arrays;
	r_arr = np.linspace(0, 2 * max(x_arr), len(x_arr));
	theta = np.linspace(0, 2 * np.pi, len(x_arr));

	# Interpolate signal(x,y)
	sig_int   = interp2d(x_arr, y_arr, signal);
	sig_polar = np.zeros((len(r_arr), len(theta)));
	for i in range(len(r_arr)):
		r_use = r_arr[i];
		for j in range(len(theta)):
			th_use = theta[j];
			# Compute x and y
			x_int = r_use * np.cos(th_use);
			y_int = r_use * np.sin(th_use);
			# Compute signal(r, theta)
			sig_polar[i, j] = sig_int(x_int, y_int);
	# Compute azimuthal lineouts of the drive and witness signals

	ind_drive = np.argmax(sig_polar[:, 0]);
	ind_wit   = np.argmax(sig_polar[ind_drive + off:-1, 0]) + ind_drive + off;
	azi_drive = np.squeeze(sig_polar[ind_drive, :]);
	azi_wit   = np.squeeze(sig_polar[ind_wit, :]);
	if plot:

		# Perform fit for plotting 
		drive_fit, popt_drive = fit_lineout(azi_drive, theta);
		wit_fit, popt_wit     = fit_lineout(azi_wit, theta);

		# Create labels 
		A_drive   = str(np.round(popt_drive[0], 2));
		phi_drive = str(np.round(popt_drive[1] / np.pi, 2));
		B_drive   = str(np.round(popt_drive[2], 2));
		drive_lab = A_drive + r'sin($\theta$' + phi_drive + r'$\pi$)'\
		            + " + " + B_drive; 
		A_wit = str(np.round(popt_wit[0], 2));
		phi_wit = str(np.round(popt_wit[1], 2));
		B_wit  = str(np.round(popt_wit[2], 2));
		wit_lab = A_wit + r'sin($\theta$' + phi_wit + r'$\pi$)' \
		            + " + " + B_wit; 

		fig1, ax1 = makefig(x = 8, y = 6, xlab = 'r [mm]',\
		           ylab = r'$\frac{\theta}{\pi}$', title = title);
		ext = np.array([min(r_arr)*1e3, max(r_arr)*1e3,\
		                min(theta) / np.pi, max(theta) / np.pi]);
		img = ax1.imshow(np.transpose(sig_polar), extent = ext, aspect = 'auto',\
		                cmap = 'CMRmap');
		cbar = plt.colorbar(mappable = img);

		fig2 = plt.figure(figsize = (8, 6), dpi = 200);
		ax2  = fig2.add_subplot(211);
		ax3  = fig2.add_subplot(212);


		# Pretty plot stuff
		plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933',\
		               '#DDCC77', '#CC6677', '#882255', '#AA4499'];
		cy  = cycler('color', plot_colors);
		ax2.set_prop_cycle(cy);
		ax2.set_title('Azimuthal lineouts', fontweight = 'bold');
		ax2.tick_params(labelsize = 'large');
		ax3.set_prop_cycle(cy);
		ax3.tick_params(labelsize = 'large');

		# Plotting lineouts
		ax2.set_xlabel(r'$\frac{\theta}{\pi}$');
		ax2.set_ylabel('Signal [AU]')
		ax2.plot(theta / np.pi, azi_drive, label = 'Drive');
		ax2.plot(theta / np.pi, drive_fit, '--r', label = "Fit: " + drive_lab);
		ax2.legend(handlelength=0, handletextpad=0, fancybox=True)
		ax2.set_xticks([])
		ax3.set_xlabel(r'$\frac{\theta}{\pi}$');
		ax3.set_ylabel('Signal [AU]')
		ax3.plot(theta / np.pi, azi_wit, label = 'Witness');
		ax3.plot(theta / np.pi, wit_fit, '--r', label = "Fit: " + wit_lab)
		ax3.legend(handlelength=0, handletextpad=0, fancybox=True)

		plt.subplots_adjust(hspace = 0.0)

		plt.show();
	return sig_polar, azi_drive, azi_wit;

