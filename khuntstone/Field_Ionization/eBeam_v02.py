import sys
sys.path.insert(0, "../")
import Constants.SI as SI
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gm
from scipy.integrate import simps

global plasmaDict
# Dictionary of plasmas and their attributes 
plasmaDict = {'Ar+' : {'Vi' : 15.75962, 'Name' :  'Ar$^{+}$', 'Z' : 1},
		   'Ar2+': {'Vi' : 27.62967, 'Name' : 'Ar$^{2+}$', 'Z' : 2},
		   'Ar3+': {'Vi' : 40.74, 'Name' : 'Ar$^{3+}$', 'Z' : 3},
		   'Ar4+': {'Vi' : 59.81, 'Name': 'Ar$^{4+}$', 'Z' : 4},
		   'Ar5+': {'Vi' : 59.81, 'Name': 'Ar$^{5+}$', 'Z' : 5},
		   'He+': {'Vi' : 59.81, 'Name': 'Ar$^{4+}$', 'Z' : 1},
		   'He2+': {'Vi': 59.81, 'Name': 'Ar$^{5+}$', 'Z' : 2}
}
	
	
def peak_charge_dens(beamParams):
	''' 
	Computes the peak charge density of a Gaussian beam

	Params:
	-------
	beamParams : dictionary
		dictionary of beam parameters 
		'sigma_z' : float, longitudinal beam size, m
		'sigma_r' : array_like, transverse beam size, m
		'sigma_t' : float, sigma_z / (beta * c)
		'beta'    : relativstic beta
		'charge'  : beam charge, C
		'en'      : beam normalized emittance, m-rad
		'beta_s'  : array_like, twiss beta at waist - corresponds to sigma_r 
					array, m
		''
	Returns:
	--------
	pPK : float
		The peak charge density of the beam
	'''
	sigma_z = beamParams['sigma_z']
	sigma_r = beamParams['sigma_r']
	Q       = beamParams['charge']
	return Q / ((2*np.pi)**(3/2) * sigma_r**2 * sigma_z)

def get_sigma_r(beamParams):
	'''
	Calculates transverse beam width from beam emittance, gamma, and beta_s and 
	appends it to beamParams.
	Params:
	-------
	beamParams as described above
	'''
	
	beamParams['sigma_r'] = \
	np.sqrt(beamParams['emitt'] * beamParams['beta_s'] /beamParams['gamma'])    
def get_pos(beamParams, nr = 10, nxi = 10, npoints = 1000):
	'''
	Quick function to create position r, xi from beamParams
	'''
	sigma_r, sigma_z = beamParams['sigma_r'], beamParams['sigma_z']                         
	r_arr = np.zeros((len(sigma_r), npoints))
	for i in range(len(sigma_r)):
		r_arr[i][:] = np.linspace(-nr * sigma_r[i], nr * sigma_r[i], npoints)
	xi_arr = np.linspace(-nxi * sigma_z, nxi * sigma_z, npoints)
	return r_arr, xi_arr, npoints

def rad_E_field(pos, beamParams, eps0 = SI.permFreeSpace, \
				c = SI.lightSpeed, peak = False):
	'''
	Computes the radial electric field of an electron beam, assumes beam is 
	transversely and radially gaussian.

	Params:
	-------
	pos : dictionary
	dictionary of position and time valus 
	'r' : r - array_like, transvers position m
	'xi' : xi - array_like, longitudinal position z minus c*t, m
	beamParams : dictionary
		beamParams as described above
	peak : Boolean, optional
		Set to true if you only want the peak E-field 

	Returns:
	--------
	Er : array_like
		3D array of electric field values at r, xi, and sigma_r, GV/m
	rPeak : float
		Approximate position of peak radial electric field um
	EPeak : float
		Approximate peak radial electric field GV/m
	'''
	sigma_z = beamParams['sigma_z']
	sigma_r = beamParams['sigma_r']
	beta    = beamParams['beta']

	# peak charge density and electric field have no spatial dependence

	ppK   =  peak_charge_dens(beamParams);
	EPeak = (ppK * sigma_r / (2*eps0)) / 1e9;
	rPeak = (np.pi * sigma_r /2) * 1e6; 
	if peak:
		return np.nan, rPeak, EPeak
	else:
		r = pos['r']
		xi = pos['xi']
		Er = np.zeros((len(sigma_r), len(r[0]), len(xi)))
		for i in range(len(sigma_r)):
			rp = np.reshape(r[i], (len(r[i]), 1))
			Er[i,:,:] = (ppK[i] * sigma_r[i]**2 / (eps0 * rp)) * \
						(1 - np.exp(-rp**2/(2*sigma_r[i]**2))) * \
						np.exp(-(xi)**2 / (2 * sigma_z**2))
			Er[i, :,:] = Er[i,:,:] / 1e9;
		return Er, rPeak, EPeak
def ionization_rate(Er, beamParams, gasName):
	''' 
	Computes the ionization rate of a neutral gas due to the radial electric 
	field of a transversely and radially Gaussian beam

	Params:
	-------
	Er_z : array_like 
		3D array of the radial electric field in sigma_r, r, and z
	Er_t : array_like
		3D array of the radial electric field in sigma_r, r, and t
	beamParams : dictionary
		dictionary of beam parameters as defined above
	gasName : str
		Name of plasma to compute ionization rate, must be key in plasmaDict
	Returns:
	--------
	W : array like
		3D array of ionization rate in sigma_r, r, and xi
	'''
	Vh = 13.6; 
	Vi = plasmaDict[gasName]['Vi']
	Z  = plasmaDict[gasName]['Z']
	n = Z / np.sqrt(Vi/Vh);
	Er = abs(Er)

	W = 1.52 * ((4**n * Vi) / (n * gm(2*n))) * \
		 (20.5 * Vi**(3/2) / Er)**(2*n-1) * \
		 np.exp(-6.83 * Vi**1.5 /Er);
	return W
def ionization_frac(W, pos, beamParams, c = SI.lightSpeed):
	'''
	Calculates the plasma density of the gas ionized by the electron beam 
	by integrating W along time

	Params:
	-------
	W : array_like
		3D array of W in sigma_r, r, and xi
	pos : dictionary
		dictionary of position arrays returned from get_pos

	Returns:
	plasma_frac : array_like
		2D array of the ionization fraction along r for each sigma_r
	max_frac : array_like
		max_ionization fraction for each sigma_r
	'''

	t = pos['xi'] / (beamParams['beta']*c);
	plasma_frac = 1 - np.exp(-simps(W * 1e15, t));
	max_frac = np.max(plasma_frac, axis = 1);
	return plasma_frac, max_frac

def plot_plasma_frac(plasma_frac, pos, beamParams, gasName, ind):
	beta_s = beamParams['beta_s'][ind]
	name = plasmaDict['gasName']['Name']
	title = 'Ionization Fraction of ' + name + ' $\\beta$ = %.2f' % beta_s
	plt.plot(pos['r'][ind]*1e6, plasma_frac[ind])
	plt.xlabel('r [$\mu$m]')
	plt.ylabel('Ionization Fraction')
	plt.title(title)
	plt.show()

def neutral_ring_width(plasma_frac, pos):
	'''
	Computes the width of the neutral gas ring for each beta_s

	Params:
	-------
	plasma_frac : array_like
		array returned in ionization frac
	pos : dictionary
		dictionary of positions returned from get_pos
	Returns:
	--------
	ring_width : array_like
		width of the neutral gas ring for each beta_s, m
	'''
	r = pos['r']
	npoints = pos['npoints']
	# Preallocate
	width = np.zeros(len(plasma_frac))
	for i in range(len(plasma_frac)):
		if max(plasma_frac[i]) == 0:
			width[i] = r[i][-1] - r[i][0]
		else:
			lower_half = plasma_frac[i][0:int(npoints/2 -1)]
			upper_half = plasma_frac[i][int(npoints/2):]
			max_lower_ind = np.argwhere(lower_half == max(lower_half))[0]
			max_upper_ind = int(np.argwhere(upper_half == max(upper_half))[0])
			lower_zero = \
				   int(np.argwhere(lower_half[int(max_lower_ind):] == 0)[0])+ \
				   int(max_lower_ind)
			upper_zero = \
					 int(np.argwhere(upper_half[0:max_upper_ind] == 0)[-1]) + \
					 int(npoints/2)
			width[i] = r[i][upper_zero] - r[i][lower_zero]
	return width

def plot_width(widths, plasmaNames, beamParams, logx = False, logy = False, \
	           log = False):
	names = [plasmaDict[i]['Name'] for i in gasNames] 
	for i in range(len(widths)):
		if logx:
			plt.semilogx(beamParams['beta_s'], widths[i]*1e6,\
						 label = names[i])
		elif logy:
			plt.semilogy(beamParams['beta_s'], widths[i]*1e6,\
						 label = names[i])
		elif log:
			plt.loglog(beamParams['beta_s'], widths[i]*1e6, \
					 label = names[i])
		else:
			plt.plot(beamParams['beta_s'], widths[i]*1e6,\
						 label = names[i])
	plt.xlabel('$\\beta$ [m]')
	plt.ylabel('Neutral gas Diameter [$\\mu$m]')
	plt.legend()
	plt.show()
def plot_max_frac(max_frac, beamParams, plasmaNames, logx = False,\
	logy = False, log = False):
	names = [plasmaDict[i]['Name'] for i in plasmaNames]
	for i in range(len(max_frac)):
		if logx:
			plt.semilogx(beamParams['beta_s'], max_frac[i], \
				label = names[i])
		elif logy:
			plt.semilogy(beamParams['beta_s'], max_frac[i], \
				label = names[i])
		elif log:
			plt.loglog(beamParams['beta_s'], max_frac[i], \
				label = names[i])
		else:
			plt.plot(beamParams['beta_s'], max_frac[i], \
				label = names[i])
	plt.xlabel('$\\beta$ [m]')
	plt.ylabel('Maximum Ionization Fraction')
	plt.legend()
	plt.show()

def plot_field(field, pos, beamParams, cbar_label, beta_s, ind, \
			   lims = [], gas = False, gasName = None, c = SI.lightSpeed):
	'''
	Plots a field in the in the rz and rt planes
	'''
	if gas:
		title = 'Ionization Rate of ' + gasName +'$\\beta$ = %.2f' % beta_s[ind]
	else:
		title = 'Radial Electric Field ' + '$\\beta$ = %.2f' % beta_s[ind] ;
	# rt plane
	r = np.flipud(pos['r'][ind]) * 1e6; nr = len(r)
	t = (pos['xi'] * 1e15 / (beamParams['beta']*c)) -\
		(pos['xi'][0]*1e15 /(beamParams['beta']*c)); 
	nt = len(t)
	if not lims:
		plt.imshow(np.flipud(field[ind]), cmap = 'jet')
	else:
		plt.imshow(np.flipud(field[ind]), cmap = 'jet',\
							 vmin = lims[0], vmax = lims[1])
	x_locs = [0, nt/2, nt]
	x_labs = [0, int(t[int(len(t)/2 -1)]), int(t[-1])]
	y_locs = [0, nr/2, nr] 
	y_labs = [int(r[0]), int(r[int(nr/2 - 1)]), \
			  int(r[-1])]
	plt.xticks(x_locs, x_labs)
	plt.yticks(y_locs, y_labs)
	cbar = plt.colorbar()
	cbar.set_label(cbar_label)
	plt.xlabel('t [fs]');
	plt.ylabel('r [$\mu$m]');
	plt.title(title)
	plt.show()

	# rz
	r = r * 1e-6;
	z = pos['xi']; nz = len(z)
	sigma_r  = beamParams['sigma_r'][ind]
	sigma_z = beamParams['sigma_z']
	if not lims:
		plt.imshow(np.flipud(field[ind]), cmap = 'jet')
	else:
		plt.imshow(np.flipud(field[ind]), cmap = 'jet', \
				   vmin = lims[0], vmax = lims[1])
	x_locs = [0, nz/2, nz]
	x_labs = [int(z[0]/sigma_z), int(z[int(len(z)/2 -1)]/sigma_z),\
			  int(z[-1]/sigma_z)]
	y_locs = [0, nr/2, nr] 
	y_labs = [int(r[0]/sigma_r), int(r[int(nr/2 - 1)]/sigma_r), \
			  int(r[-1]/sigma_r)]
	plt.xticks(x_locs, x_labs)
	plt.yticks(y_locs, y_labs)
	cbar = plt.colorbar()
	cbar.set_label(cbar_label)
	plt.xlabel('z/$\sigma_z$');
	plt.ylabel('r/$\sigma_r$');
	plt.title(title)
	plt.show()
