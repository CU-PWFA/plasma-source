import sys
sys.path.insert(0, "../")
import Constants.SI as SI
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gm
from scipy.integrate import simps
import matplotlib.patches as patches
global plasmaDict
# Dictionary of plasmas and their attributes 
plasmaDict = {'Ar+' : {'Vi' : 15.75962, 'Name' :  'Ar$^{+}$', 'Z' : 1},
		   'Ar2+': {'Vi' : 27.62967, 'Name' : 'Ar$^{2+}$', 'Z' : 2},
		   'Ar3+': {'Vi' : 40.74, 'Name' : 'Ar$^{3+}$', 'Z' : 3},
		   'Ar4+': {'Vi' : 59.81, 'Name': 'Ar$^{4+}$', 'Z' : 4},
		   'Ar5+': {'Vi' : 59.81, 'Name': 'Ar$^{5+}$', 'Z' : 5},
		   'He+': {'Vi' : 59.81, 'Name': 'He$^{+}$', 'Z' : 1},
		   'He2+': {'Vi': 59.81, 'Name': 'He$^{2+}$', 'Z' : 2}
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
def rad_E_field__sigma_z(pos, beamParams, peak = False,
						 eps0 = SI.permFreeSpace, c = SI.lightSpeed):
	# Same as above but with fixed sigma_r and array of sigma_z
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
		Er = np.zeros((len(sigma_z), len(r), len(xi[0])))
		rp = np.reshape(r, (len(r), 1))
		for i in range(len(sigma_z)):
			
			Er[i,:,:] = (ppK[i] * sigma_r**2 / (eps0 * rp)) * \
						(1 - np.exp(-rp**2/(2*sigma_r**2))) * \
						np.exp(-(xi[i])**2 / (2 * sigma_z[i]**2))
			Er[i, :,:] = Er[i,:,:] / 1e9;
		return Er, rPeak, EPeak
def ionization_rate(Er, gasName):
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
def ionization_frac_sigma_z(W, pos, beamParams, c = SI.lightSpeed):
	# Same as above but with varying sigma_z, only returns max_frac
	t = pos['xi']  / (beamParams['beta']*c)
	max_frac = np.zeros(len(W))
	for i in range(len(W)):
		plasma_frac = 1 - np.exp(-simps(W[i] * 1e15, t[i]))
		max_frac[i] = np.amax(plasma_frac)
	return max_frac 

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
	names = [plasmaDict[i]['Name'] for i in plasmaNames] 
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
def plot_max_frac(max_frac, beamParams, plasmaNames, fs = 12, lw = 1,\
				  logx = False, logy = False, log = False):
	names = [plasmaDict[i]['Name'] for i in plasmaNames]
	for i in range(len(max_frac)):
		if logx:
			plt.semilogx(beamParams['beta_s'], max_frac[i], \
				label = names[i], linewidth = lw)
		elif logy:
			plt.semilogy(beamParams['beta_s'], max_frac[i], \
				label = names[i], linewidth = lw)
		elif log:
			plt.loglog(beamParams['beta_s'], max_frac[i], \
				label = names[i], linewidth = lw)
		else:
			plt.plot(beamParams['beta_s'], max_frac[i], \
				label = names[i], linewidth = lw)
	plt.xlabel('$\\beta$ [m]', fontsize = fs)
	plt.ylabel('Maximum Ionization Fraction', fontsize = fs)
	lg = plt.legend()
	lg.fontsize = fs - 2
	ax = plt.gca()
	ax.tick_params(labelsize = fs - 2)
	plt.show()

def plot_field(field, pos, beamParams, cbar_label, ind = 0, fs = 12, \
			   lw = 1,lims = [], gas = False, \
			   gasName = None, c = SI.lightSpeed):
	'''
	Plots a field in the in the rz and rt planes
	'''
	field = abs(field)
	if gas:
		title = 'Ionization Rate of ' + gasName 
	else:
		title = 'Radial Electric Field ' 
	#beta_str = 'Beta = %.2f' % beta_s[ind] + 'm' ;

	# rt plane
	r = np.flipud(pos['r'][ind]) * 1e6; nr = len(r)
	t = (pos['xi'] * 1e15 / (beamParams['beta']*c)) -\
		(pos['xi'][0]*1e15 /(beamParams['beta']*c)); 
	ext = [min(t), max(t), min(r), max(r)]
	fig1 = plt.figure(); ax1 = fig1.gca()
	cen = (int(t[-1]/2),0)
	ax1.add_artist(patches.Ellipse(cen, beamParams['sigma_t']*2e15, \
				   beamParams['sigma_r'][ind]*2e6, fc = 'none',\
				   ls = '--', lw = lw, ec = 'k'))
	img1 = ax1.imshow(np.flipud(field[ind]), cmap = 'jet',aspect = 'auto', \
		extent = ext)
	cbar1 = plt.colorbar(mappable = img1, ax = ax1)
	cbar1.set_label(cbar_label, fontsize = fs)
	cbar1.ax.tick_params(labelsize = fs - 2)
	# draw ellipse of beam
	
	
	ax1.set_xlabel('t [fs]', fontsize = fs);
	ax1.set_ylabel('r [$\mu$m]', fontsize = fs);
	ax1.set_title(title, fontsize = fs)
	ax1.tick_params(axis = 'both', labelsize = fs - 2)
	plt.show()

	# rz
	z = pos['xi'] * 1e6;
	ext = [min(z), max(z), min(r), max(r)]
	fig2 = plt.figure()
	ax2 = fig2.gca()
	ax2.add_artist(patches.Ellipse((0,0), beamParams['sigma_z']*2e6, \
				   beamParams['sigma_r'][ind]*2e6, fc = 'none',\
				   ls = '--', lw = lw, ec = 'k'))
	img2 = ax2.imshow(np.flipud(field[ind]), cmap = 'jet', aspect = 'auto', \
		extent = ext)
	cbar2 = plt.colorbar(mappable = img2, ax = ax2)
	cbar2.set_label(cbar_label, fontsize = fs)
	cbar2.ax.tick_params(labelsize = fs)
	
	ax2.set_xlabel('z [$\mu$m]', fontsize = fs);
	ax2.set_ylabel('r [$\mu$m]', fontsize = fs);
	ax2.set_title(title, fontsize = fs)
	ax2.tick_params(axis = 'both', labelsize = fs - 2)
	plt.show()

def plot_2D_plasma(W, pos, beamParams, gasName, ind = 0, \
	               lw = 1, fs = 12, c = SI.lightSpeed):
	title = gasName + ' plasma'
	r = pos['r'][ind];
	xi = pos['xi']
	t = (xi  / (beamParams['beta']*c)) -\
		(xi[0] /(beamParams['beta']*c)); 
	W_int = np.fliplr(np.cumsum(W[ind], axis = 1)) * ((t[1]-t[0]) * 1e15)
	n_rz = 1 - np.exp(-W_int)
	fig1 = plt.figure(); ax1 = fig1.gca()
	ext = [min(xi)*1e6, max(xi)*1e6, min(r)*1e6, max(r)*1e6]
	img = ax1.imshow(n_rz, cmap = 'jet', extent = ext, aspect = 'auto')
	ax1.add_artist(patches.Ellipse((0,0), beamParams['sigma_z']*2e6, \
				   beamParams['sigma_r'][ind]*2e6, fc = 'none',\
				   ls = '--', lw = lw, ec = 'k'))
	ax1.set_xlabel('z [$\mu$m]', fontsize = fs)
	ax1.set_ylabel('r [$\mu$m]', fontsize = fs)
	ax1.set_title(title, fontsize = fs)
	ax1.tick_params(axis = 'both', labelsize = fs - 2)
	cbar = plt.colorbar(mappable = img, ax = ax1)
	cbar.ax.tick_params(labelsize = fs)
	plt.show()