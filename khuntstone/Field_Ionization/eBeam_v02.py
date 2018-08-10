import sys
sys.path.insert(0, "../")
import Constants.SI as SI
import numpy as np
from scipy.special import gamma as gm
from scipy.integrate import simps
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#plt.style.use('presentation')
global plasmaDict, c
c = SI.lightSpeed
# Dictionary of plasmas and their attributes 
plasmaDict = {'Ar+' : {'Vi' : 15.75962, 'Name' :  'Ar$^{+}$', 'Z' : 1},
		   'Ar2+': {'Vi' : 27.62967, 'Name' : 'Ar$^{2+}$', 'Z' : 2},
		   'Ar3+': {'Vi' : 40.74, 'Name' : 'Ar$^{3+}$', 'Z' : 3},
		   'Ar4+': {'Vi' : 59.81, 'Name': 'Ar$^{4+}$', 'Z' : 4},
		   'Ar5+': {'Vi' : 59.81, 'Name': 'Ar$^{5+}$', 'Z' : 5},
		   'He+': {'Vi' : 59.81, 'Name': 'He$^{+}$', 'Z' : 1},
		   'He2+': {'Vi': 59.81, 'Name': 'He$^{2+}$', 'Z' : 2}
}
def get_sigma_r(beamParams):
	'''
	Calculates transverse beam width from beam emittance, gamma, and beta_s and 
	appends it to beamParams.
	Params:
	-------
	beamParams : dictionary
		Dictionary of beam parameters (more detail below)
	'''
	
	beamParams['sigma_r'] = \
	np.sqrt(beamParams['emitt'] * beamParams['beta_s'] /beamParams['gamma'])

def get_beam(gamma = 20000.,en = 5.3e-6, beta_s = np.array([.1]),\
            sigma_z = 5.2e-6, Q = 1.5e-9, ):
	'''
	Creates a dictionary of beam params based on input values. Defualt is a 
	Facet II like beam

	Parameters:
	-----------
	gamma : float, optional
		Relativistic factor for the centroid beam energy, default 20000.
	en : float, optional
		Normalized emittance of the beam [mm-rad], default 5.3e-6
	beta_s : array_like, optional
		Array of beam waist betas [m], default [.1]
	sigma_z : float, optional
		Longitudinal beam size [m], default 5.2e-6
	Q : float, optional
		Charge of the beam [C], default 1.5e-9

	Returns:
	--------
	beamParams : dictionary
		Dictionary of beam parameters
	'''

	beamParams = {
              'gamma'   : gamma,  
              'sigma_z' : sigma_z,  
              'charge'  : Q, 
              'emitt'   : en, 
              'beta_s'  : beta_s
             }
    # get sigma_r from current gamma, beta_s, and en
	get_sigma_r(beamParams)
	
	# get beta from gamma
	beamParams['beta'] = np.sqrt(1 - 1/gamma**2); 
	# get sigma t from beta and sigma_z
	beamParams['sigma_t'] = sigma_z / (beamParams['beta'] * c)
	return beamParams

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

    
def get_pos(beamParams, nr = 10, nxi = 10, npoints = 1000):
	'''
	Quick function to create position r, xi from beamParams
	'''
	sigma_r, sigma_z = beamParams['sigma_r'], beamParams['sigma_z']                         
	r_arr = np.zeros((len(sigma_r), npoints))
	for i in range(len(sigma_r)):
		r_arr[i][:] = np.linspace(-nr * sigma_r[i], nr * sigma_r[i], npoints)
	xi_arr = np.linspace(-nxi * sigma_z, nxi * sigma_z, npoints)
	pos = {'r' : r_arr, 'xi' : xi_arr, 'npoints' : npoints}
	return pos

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

