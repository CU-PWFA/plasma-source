#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studing the vertical resolution of spatial encoding of EOS-BPM
Created on Wed Sep 25 14:13:43 2019

@author: keenan
"""

import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c, epsilon_0;
eps0 = epsilon_0;
import sys;
sys.path.insert(0, "../../python/");
from crystal import crystal;
from ebeam import ebeam;
from laser import laser;
import phase_retard as pr;
from plotting import makefig;
import thz;

def get_2D_signal(drive, wit, probe, psi, ctype, d, x0, x, y0, y, td, \
				  nx = 1000, ny = 1000, sweep = 'horizontal'):
	"""
	Function to calculate the 2D EOS-BPM signal (i.e. imaged crystal)
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
	sweep : str, optional
			The direction in which the probe is sweeping, default horizontal
			Options are 'horizontal' (sweep in x), 'vertical' (sweep in y)
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
	cry   = crystal(ctype);

	# Create spatial and temporal arrays
	x_arr = np.linspace(0, x, nx);
	y_arr = np.linspace(0, y, ny);
	if sweep == 'horizontal':
		tau = (x_arr / c) * np.tan(psi * np.pi / 180);
	elif sweep == 'vertical':
		tau = (y_arr / c) * np.tan(psi * np.pi / 180);
	y_arr = np.linspace(-y/2, y/2, ny);
	Et_drive            = np.exp(-(drive.t**2) / (2 * drive.sigt**2));
	FErt_drive, f_drive = thz.raw_field(Et_drive, drive.t);
	Ect_drive, tt_drive = thz.cry_field(drive.t, FErt_drive, f_drive, d, \
										probe, cry);

	Et_wit          = np.exp(-(wit.t + wit.del_t)**2 / (2 * wit.sigt**2));
	FErt_wit, f_wit = thz.raw_field(Et_wit, wit.t);
	Ect_wit, tt_wit = thz.cry_field(wit.t, FErt_wit, f_wit, d, probe, cry);

	base_drive, dm  = pr.phase_retard(Ect_drive, tt_drive * 1e-12, d, tau + td,\
									  probe, cry, 'spatial', psi = psi);
	base_wit, dm    = pr.phase_retard(Ect_wit, tt_wit * 1e-12, d, tau + td, \
									  probe, cry, 'spatial', psi = psi);

	gamma_drive = np.zeros((len(x_arr), len(y_arr)));
	gamma_wit   = np.zeros((len(x_arr), len(y_arr)));


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

	        if sweep == 'horizontal':
	        	gamma_drive[i, j] = base_drive[i] * E0_drive;
	        	gamma_wit[i, j]   = base_wit[i] * E0_wit;
	        elif sweep == 'vertical':
	        	gamma_drive[i, j] = base_drive[j] * E0_drive;
	        	gamma_wit[i, j]   = base_wit[j] * E0_wit;

	gamma = gamma_drive + gamma_wit;
	return gamma_drive, gamma_wit, gamma;


def integrate_1D(sig, x, y, axis, cutoff):
	"""
	Function to get the 1D integrated signal
	Paramters:
	----------
	sig    : array_like
			 The 2D signal
	axis   : int
			 The axis of integration
	cutoff : float
			 Cutoff of integrating the signal (position in crystal)
	Returns:
	--------
	int_sig : array_like
			  The integrated signal array
	"""
	if axis == 0:
		ind = np.argwhere(x <= cutoff)[-1];
		new_sig = sig[0:int(ind), :]
		int_sig = np.sum(new_sig, axis = 0);
	elif axis == 1:
		ind = np.argwhere(y <= cutoff)[-1];
		new_sig = sig[:, 0:int(ind)];
		int_sig = np.sum(new_sig, axis = 1);
	return int_sig;  
def plot_2D_sig(sig, ext, save = False, sname = ''):
	"""
	Function to plot the 2D EOS signal
	
	Parameters:
	-----------
	sig   : array_like
			2D array of signal
	ext   : array_like
			The extent array for plt.imshow() [xmin, xmax, ymin, ymax] in mm
	save  : bool, optional
			Whether or not to save the figure
	sname : str, optional
			The filename if save == True
	"""
	
	fig, ax = makefig(x = 8, y = 6, xlab = 'x [mm]', ylab = 'y [mm]');
	
	img = ax.imshow(np.transpose(sig), extent = ext, \
					aspect = 'auto', cmap = 'CMRmap');
	plt.colorbar(mappable = img, label = 'Signal [AU]');
	if save:
		fig.savefig(sname);
	plt.show();
	
									   
	