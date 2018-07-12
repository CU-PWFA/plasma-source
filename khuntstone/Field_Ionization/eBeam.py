# Suite of functions useful for plasma production via field Ionization
# Assume Gaussian electron beam for beamParams
import sys
sys.path.insert(0, "../")
import Constants.SI as SI
import numpy as np
def peak_charge_dens(beamParams):
	''' 
	Computes the peak charge density of a Gaussian beam

	Params:
	-------
	beamParams : array_like
		Array of beam parameters 
		[0] - sigma_z, longitudinal beam size m
		[1] - sigma_r, transverse beam size m
		[2] - beta, relativistic factor of beam (no energy spread)
		[3] - Q, beam charge

	Returns:
	--------
	pPK : float
		The peak charge density of the beam
	'''
	sigma_z = beamParams[0]
	sigma_r = beamParams[1]
	Q       = beamParams[3]
	pPK = Q / ((2*np.pi)**(3/2) * sigma_r**2 * sigma_z)
	return pPK
def rad_E_field(pos, beamParams, t = 0.0, eps0 = SI.permFreeSpace,\
                c = SI.lightSpeed, peak = False):
	'''
	Computes the radial electric field of a Gaussian eBeam
	
	Params:
	-------
	pos : array_like
	array of position and time valus 
	[0] r - array_like, transvers position m
	[1] z - array_like, longitudinal position m
	
	t : float, optional
		time in seconds, default is 0
		
	beamParams : array_like
		Array of beam parameters
	peak : Boolean, optional
	    Set to true if you only want the peak E-field 

	Returns:
	--------
	Er : array_like
		2D array of electric field values at r and z
	rPeak : float
		Approximate position of peak radial electric field
	EPeak : float
		Approximate peak radial electric field
	'''
	sigma_z = beamParams[0]
	sigma_r = beamParams[1]
	beta    = beamParams[2]
	r       = pos[0]
	z       = pos[1]
	# peak charge density and E-field has no positional dependence
	pPK = peak_charge_dens(beamParams);
	EPeak = pPK * sigma_r / (2*eps0);
	rPeak = np.pi * sigma_r / 2;
	# Preallocate for loop
	if peak:
	    return np.nan, rPeak, EPeak
	else:
	    Er = np.zeros((len(r), len(z)));
	
	    for i in range(len(r)):
		    for j in range(len(z)):
			    Er[i,j] = (pPK * sigma_r**2 / (eps0 * r[i])) * \
			              (1 - np.exp(-r[i]**2/(2*sigma_r**2))) * \
			              np.exp(-(z[j] - beta*c*t)**2 / (2*sigma_z**2))
			    # Put Er in GV/m
			    Er[i,j] = Er[i,j] / 1e9;
	
	    return Er, rPeak, EPeak;
