#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:43:09 2019
Module for interpolating particle velocity and position in a plasma for betatron
radiation calculations
@author: keenan
"""

# Standard imports
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const



def velocity_interp(v, tau, tau_int, n_pad = 2):
	""" 
	Function to perform a linear 1D interpolation of a particle's velocity
		
	Parameters
	----------
	v : array_like
		Array of particle's velocity
	tau : array_like
		  Proper time array corresponding to v
	tau_int : array_like
		Array of query proper times for interpolation
	Returns
	-------
	v1n : array_like
		Array of interpolation coefficients for every time step in tau
	v_int : array_like
		Array of interpolated velocities for each time step in tau
	"""
	
	# Get indices for points in tau surrounding points in tau_int
	
	ind1 = np.squeeze(np.array([np.argwhere(tau >=i)[0] - 1 for i in tau_int]))
	ind1[np.argwhere(ind1 < 0)] = 0
	ind2 = ind1 + 1
	v0   = v[ind1]
	v1n   = (v[ind2] - v[ind1]) / (tau[ind2] - tau[ind1])

	v_int = v0 + v1n * (tau_int - tau[ind1]);
	v1n   = np.unique(v1n)
	if len(v1n) != len(tau):
		v1n = np.zeros(len(tau)) + v1n[0]
	v1n = v1n[0:-int(n_pad)]
	return v1n, v_int

def position_interp(x, tau, tau_int, n_pad = 2):
	"""
	Function to perform a quadratic 1D interpolation of a particle's position
	
	Parameters
	----------
	x : array_like
		Array of particle's position
	tau : array_like
		Array of proper times corresponding to x
	tau_int : array_like
		Array of query proper times for the interpolation
	n_pad : int
		The number of padded entries to tau

	Returns
	-------
	x0n : array_like
		Array of zeroth inteprolation coefficient (x(tau_n))
	x1n : array_like
		Array of 1st interpolation coefficient for every time step in tau
	x2n : array_like
		Array of 2nd interpolation coefficient for every time step in tau
	x_int : array_like
		Array of interpolated positions corresponding to tau_int
	"""

	# Preallocate for loop
	x_int   = np.zeros(len(tau_int))
	x0n     = np.zeros(len(tau))
	x1n     = np.zeros(len(tau))
	x2n     = np.zeros(len(tau))
	
	# Get time step (MUST BE UNIFORM)
	dtau = tau[1] - tau[0]
	cf_ind = 0
	init_ind = 0
	for i in range(len(tau_int)):
		tau_in = tau_int[i]
		inds   = np.argwhere(tau >= tau_in)[0]
		if inds > 0:
			ind1 = inds[0] - 1
		else:
			ind1 = inds[0]
		if ind1 + 1 == len(x) - 1:
			ind1 = ind1-1
		ind2 = ind1 + 1
		ind3 = ind2 + 1

		tau_l = tau[ind1]
		tau_c = tau[ind2]
		tau_u = tau[ind3]

		x_l   = x[ind1]
		x_c   = x[ind2]
		x_u   = x[ind3]

		d0    = x_c
		d2    = (x_u - x_c + x_l - x_c)
		if abs(d2) <= dtau:
			d2 = 0.0
		else:
			d2 = d2 / dtau / dtau

		d1 = (x_u - x_c) / (tau_u - tau_c) - 0.5 * d2 * (tau_u - tau_c)

		if i == 0:
			x0n[cf_ind] = d0
			x1n[cf_ind] = d1
			x2n[cf_ind] = 0.5 * d2
			cf_ind += 1
		if ind1 > init_ind:
			x0n[cf_ind] = d0
			x1n[cf_ind] = d1
			x2n[cf_ind] = 0.5 * d2
			cf_ind      += 1
			init_ind    = ind1
		x_int[i] = d0 + d1 * (tau_in - tau_c) + 0.5 * d2 * (tau_in - tau_c)**2

	x0n = x0n[0:-int(n_pad)]
	x1n = x1n[0:-int(n_pad)]
	x2n = x2n[0:-int(n_pad)]
	return x0n, x1n, x2n, x_int