import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit as cf
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "../../python/")
sys.path.insert(0,"../")
from Constants import SI 
from vsim import load 
from vsim import analyze
import h5py
global kb, me, c
kb = SI.boltzmann; # Boltzmann constant in J/K
me = SI.elecMass; # electron mass in kg
c  = SI.lightSpeed; # speed of light in m/s
def get_ArIonize_data(ionize_start, double_start,\
                      nfiles, simPath, simName):
	''' Function that gets data for edgeE, NeutralGas, ionized_electrons, Ar1, 
		and Ar2 from every dump of VSim Argon ionization code.

		Parameters:
		-----------
		ionization_start : int
			dump number when ionization began
		double_start : int
			dump number when double ionization began
		nfiles : int
			number of dumps 
		simPath : string
			file path to the simulation
		simName : string
			name of the simulation
		Returns:
		--------
		edgeE_data : array_like
			Array of data for the edge E-field
		gas_data : array_like
			Array of data for the neutral gas source
		elec_data : array_like
			Array of data for the ionized electrons
		Ar1_data : array_like
			Array of data for Ar1 ions
		Ar2_data : array_like
			Array of data for Ar2 ions
		t : array_like
			Array of times of each dump
	'''
	edgeE_data  = []
	gas_data    = []
	elec_data   = []
	Ar1_data    = []
	Ar2_data    = []
	attrs       = []
	for i in range(nfiles):
		edgeE_data.append(load.get_field_data(simPath + simName + \
						  '_edgeE_' + str(i) + '.h5', 'edgeE'))
		attrs.append(load.get_field_attrs(simPath + simName + \
						  '_edgeE_' + str(i) + '.h5', 'edgeE'))
		gas_data.append(load.get_field_data(simPath + simName + \
						'_NNeutralGas_' + str(i) + '.h5', 'NNeutralGas'))
		if i >= ionize_start:
			elec_data.append(load.get_species_data(simPath + simName + \
				   '_ionizedElectrons_' + str(i) + '.h5', 'ionizedElectrons'))
			Ar1_data.append(load.get_species_data(simPath + simName + \
										   '_Ar1_' + str(i) + '.h5', 'Ar1'))
		if i >= double_start:
			Ar2_data.append(load.get_species_data(simPath + simName + \
										   '_Ar2_' + str(i) + '.h5', 'Ar2'))
	t = np.zeros(len(attrs))
	for i in range(len(attrs)):
		t[i] = attrs[i]['time']
	return np.array(edgeE_data), np.array(gas_data), np.array(elec_data),\
	       np.array(Ar1_data), np.array(Ar2_data), t

def get_ArRecomb_data(recomb_start, nfiles, simPath, simName):
	''' Function that gets data for, Ar0, electrons, Ar1, 
		and Ar2 from every dump of VSim Argon ionization code.

		Parameters:
		-----------
		recomb_start : int
			dump number when recombination began
		nfiles : int
			number of dumps 
		simPath : string
			file path to the simulation
		simName : string
			name of the simulation

		Returns:
		--------
		elec_data : array_like
			Array of data for the ionized electrons
		Ar0_data : array_like
			Array of data for neutral Ar0
		Ar1_data : array_like
			Array of data for Ar1 ions
		Ar2_data : array_like
			Array of data for Ar2 ions
		t : array_like
			Array of times at each dump
	'''
	elec_data   = []
	Ar0_data    = []
	Ar1_data    = []
	Ar2_data    = []
	attrs       = []
	for i in range(nfiles):
		elec_data.append(load.get_species_data(simPath + simName + \
				   '_electrons_' + str(i) + '.h5', 'electrons'))
		attrs.append(load.get_field_attrs(simPath + simName + \
						  '_edgeE_' + str(i) + '.h5', 'edgeE'))
		Ar1_data.append(load.get_species_data(simPath + simName + \
										   '_Ar1_' + str(i) + '.h5', 'Ar1'))
		Ar2_data.append(load.get_species_data(simPath + simName + \
										   '_Ar2_' + str(i) + '.h5', 'Ar2'))
		if i >= recomb_start:
			Ar0_data.append(load.get_species_data(simPath + simName + \
				'_Ar0_' + str(i) + '.h5', 'Ar0'))	
			
	t = np.zeros(len(attrs))
	for i in range(len(attrs)):
		t[i] = attrs[i]['time']
	# Cast to numpy arrays for faster analysis.
	return np.array(elec_data), np.array(Ar0_data), np.array(Ar1_data), \
	       np.array(Ar2_data), t

def double_peak(v, A, sig, B):
	'''
	Function to calculate two Gaussians of the same standard deviation offset 
	from 0 by the same amount in opposite direction. Used to fit species VDFs.

	Paramters:
	----------
	v : array_like
		Array of velocities 
	A : float
		The peak of Gaussian, important for getting a good fit not for analysis
	sig : float, 
		the standard deviation of the Gaussian. 
	B : float
		The offset of the Gassians from 0. Usually the absolute value of the 
		mean of v. 
	'''
	return ((A) * (2*np.pi)**(-.5)) * (np.exp(-(v-B)**2 / (2*sig**2)) + \
		    np.exp(-(v+B)**2 / (2*sig**2)))
def correction(uy):
	'''
	corrects the y VDF to be centered at 0 when strong a0 causes a double peak 
	distribution

	Parameters:
	-----------
	uy : array_like
		the uncorrected y velocities

	Returns:
	--------
	uy_correct : array_like
		the corrected y-velocities
	'''
	uy_correct = np.zeros(len(uy))
	mean = np.mean(abs(uy))
	for i in range(len(uy)):
		if uy[i] < 0:
			uy_correct[i] = uy[i] + mean
		else:
			uy_correct[i] = uy[i] - mean
	return uy_correct
def get_T(species, ndim = 2, correct_y = False):
	'''
	Function to compute the temperature of a species from its VDF. 

	Parameters:
	-----------
	species : array_like
		array of data returned from vsim.load for the species (electrons, ions
		etc...)
	ndim : boolean, otpional
		Number of dimensions of the ionizing simulation
	correct_y : boolean, optional
		Whether or not to adjust the y -VDF for performing 2D or 3D temperature
		determination
	Returns:
	--------
	T_1D : float
		The temperature of the species using the polarization direction VDF
	T_2D : float
		The temperature from the 2D VDF
	T_3D : float
		The temperature from the 3D VDF
	'''
	# species velocity data
	if ndim == 2:
		ux = analyze.get_ux(species)
		uy = analyze.get_uy(species)
		gammas = analyze.get_ptc_gamma(species)
		ux = ux/gammas;
		uy = uy/gammas;
		T_3D = np.nan
	if ndim == 3:
		ux = species[:,3]
		uy = species[:,4]
		uz = species[:,5]
		u2 = ux**2 + uy**2 + uz**2
		gammas = 1/np.sqrt(1 - u2/(c**2 + u2))
		ux = ux/gammas;
		uy = uy/gammas;
		uz = uz/gammas;

	# particle weights
	weights = analyze.get_weights(species)

	# 1D temperature fit
	# nominal guesses to help the fit 
	p0 = [7e-8, .5e7, .5e7]
	# histogram of uy
	n, bins, patches = plt.hist(uy, 200, weights = weights, normed = True)
	popt, pcov = cf(double_peak, bins[1:], n, p0 = p0)
	sigma_1D = popt[1]
	T_1D = (me * sigma_1D**2/kb)/11600;

	# 2D and 3D temperature fit
	if correct_y:
		uy_fit = correction(uy)
	else:
		uy_fit = uy
	u2d = np.sqrt(ux**2 + uy_fit**2)
	(mu, sigma_2D) = norm.fit(u2d)
	T_2D = (me * sigma_2D**2/kb)/11600;
	if ndim == 3:
		u3d = np.sqrt(ux**2 + uy_fit**2 + uz**2)
		# make isotropic
		um = -u3d
		u3d = um.tolist() + u3d.tolist
		(mu, sigma_3D) = norm.fit(u3d)
		T_3D = (me * sigma_3D**2/kb)/11600;

	return T_1D, T_2D, T_3D

	
