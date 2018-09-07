import numpy as np
from scipy.stats import norm, chisquare
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
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
		if i >= ionize_start + 1:
			Ar2_data.append(load.get_species_data(simPath + simName + \
										   '_Ar2_' + str(i) + '.h5', 'Ar2'))
	return np.array(edgeE_data), np.array(gas_data),\
		   np.array(elec_data), np.array(Ar1_data), np.array(Ar2_data), attrs

def alpha_colormap(cmap, cutoff, flip=True):
	N = cmap.N
	cmapt = cmap(np.arange(N))
	alpha = np.ones(N)
	if flip:
		temp = alpha[:int(cutoff*N)]
		M = len(temp)
		alpha[:int(cutoff*N)] = np.linspace(0, 1, M)
	else:
		alpha[int((1-cutoff)*N):] = 0.0
	cmapt[:, -1] = alpha
	cmapt = colors.ListedColormap(cmapt)
	return cmapt

def plot_gas_evolution(gas, pulse, t, LX, LY, savePath, aniName):
	'''
	Plots the gas density and laser pulse throughout the simulation

	Parameters:
	-----------
	gas : array_like
		data of the gas to be plotted
	pulse : array_like
		data of the laser pulse to be plotted
	t : array_like
		time at each dump used in the plot in fs
	LX : int
		longitudinal length of the simulation (um)
	LY : int 
		transverse length of the simulation (um)
	savePath : string
		path to save the animation to
	aniName : string
		name of the animation
	'''
	gas = np.squeeze(gas)
	fig = plt.figure(figsize = (5,3.5), dpi = 150)
	ax  = fig.gca()
	ext = np.array([-LX/2, LX/2, -LY/2, LY/2])

	
	# y-intensity (for pulse)
	pulse = pulse[:,:,:,1]**2
	pulse = pulse.transpose()
	cmap_p = alpha_colormap(plt.cm.get_cmap('viridis'), .01 * np.amax(pulse))
	im = plt.imshow(gas[0].transpose(), extent = ext, cmap = 'jet', \
					aspect = 'auto', animated = True)
	im2 = plt.imshow(pulse[:,:,0], extent = ext,\
					 aspect = 'auto', animated = True, alpha = .5, \
					 cmap = cmap_p, \
					 vmin = np.amin(pulse), vmax = np.amax(pulse))
					 
	plt.xlabel(r'x [$\mu$m]')
	plt.ylabel(r'y [$\mu$m]')
	cbar = plt.colorbar(mappable = im)
	cbar.set_label(r'm$^{-3}$')
	def updatefig(i, gas, pulse):
		im.set_data(gas[i].transpose())
		im2.set_data(pulse[:,:,i])
		plt.title('t = %.2f' %t[i] + ' fs \n Neutral gas density')	
	ani = animation.FuncAnimation(fig, updatefig, frames = len(gas), \
					fargs = (gas, pulse), repeat = True, interval = 500)
	ani.save(savePath + aniName)
def plot_pulse(pulse, t, LX, LY, savePath, aniName):
	fig = plt.figure(figsize = (5,3.5), dpi = 150)
	ax  = fig.gca()
	ext = np.array([-LX/2, LX/2, -LY/2, LY/2])
	pulse = pulse[:,:,:,1]#**2
	pulse = pulse.transpose()
	im = plt.imshow(pulse[:,:,0], extent = ext,\
					 aspect = 'auto', animated = True, cmap = 'jet', \
					 vmin = np.amin(pulse), vmax = np.amax(pulse) )
	plt.xlabel(r'x [$\mu$m]')
	plt.ylabel(r'y [$\mu$m]')
	#cbar = plt.colorbar(mappable = im)
	#cbar.set_label(r'm$^{-3}$')
	def updatefig(i, pulse = pulse):
		im.set_array(pulse[:,:,i])
		plt.title('t = %.2f' %t[i] + ' fs \n Laser pulse')	
	ani = animation.FuncAnimation(fig, updatefig, frames = len(pulse[0,0,:]), \
					repeat = True, interval = 500)
	ani.save(savePath + aniName)

def vdf_fit(species):
	''' Plots the VDF (total and in individual directions) at the end of the
	    simulation
	'''
	ux = species[:,3]
	uy = species[:,4]
	uz = species[:,5]
	weights = species[:,-1]
	fig1 = plt.figure(); fig2 = plt.figure(); fig3 = plt.figure;
	fig4 = plt.figure(); ax1 = fig1.gca(); ax2 = fig2.gca(); 
	ax3  = 