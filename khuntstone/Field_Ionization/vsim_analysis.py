import numpy as np
from scipy.stats import norm, chisquare
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
sys.path.insert(0, "../../python/")
sys.path.insert(0, "../")
from Constants import SI
import vsim.load as load
import vsim.analyze as analyze
import h5py
import matplotlib.mlab as mlab
global kb, me
kb = SI.boltzmann; # J/K
me  = SI.elecMass; # kg
plt.style.use('notes')

def get_data(ionize_start, nfiles, simPath, simName):
	''' Function that gets data for edgeE, NeutralGas, ionized_electrons, Ar1, 
	    and Ar2 from every dump of VSim Argon ionization code.

	    Parameters:
	    -----------
	    ionization_start : int
	    	dump number when ionization began
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
	'''
	edgeE_data = []
    gas_data = []
    elec_data =  []
    Ar1_data =   []
    Ar2_data =   []
    for i in range(nfiles):
        edgeE_data.append(load.get_field_data(simPath + simName + \
        	              '_edgeE_' + str(i) + '.h5', 'edgeE'))
        gas_data.append(load.get_field_data(simPath + simName + \
        	            '_NNeutralGas_' + str(i) + '.h5', 'NNeutralGas'))
    	if i >= ionize_start:
    		elec_data.append(load.get_species_data(simPath + simName + \
                   '_ionizedElectrons_' + str(i) + '.h5', 'ionizedElectrons'))
        Ar1_data.append(load.get_species_data(simPath + simName + \
                                           '_Ar1_' + str(i) + '.h5', 'Ar1'))
        if i >= ionize_start + 1:
        	Ar2_data.append(load.get_species_data(simPath + simName + \
                                           '_Ar2_' + str(i) + '.h5', 'Ar2'))
    return edgeE_data, gas_data, elec_data, Ar1_data, Ar2_data2
