import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit as cf
import ionize_analysis as ia
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib.ticker as tkr
import sys
sys.path.insert(0, "../../python/")
sys.path.insert(0, "../")
import vsim.load as load
import vsim.analyze as analyze
from Constants import SI
import h5py
import matplotlib.mlab as mlab
import ionize_analysis as ia
global kb, me
plt.style.use('Presentation')
kb = SI.boltzmann;
me = SI.elecMass;

def plot_EDF_fit(v, weights, dim, name, pol = False, path = '', save = False):
	'''
	Function to plot a 1D fit of the energy distribution function

	Parameters:
	v : array_like
		velocities of the particles in m/s
	weight : array_like
		weights of the particles
	dim : str
		the dimension being plotted 
	name : str
		name of the figure
	pol : bool, optional
		whether or not the velocity is in the polarization direction, 
		in which case a custom fit is used. 
	path : str, optional
		path to save the figure to if wanted
	save : bool, optional
		whether or not to save the figure
	'''
	fig = plt.figure(figsize = (6,4), dpi = 150)
	ax  = fig.gca()
	n, bins, patches = ax.hist(v, 200, weights = weights, normed = True)
	if pol:
		p0 = [7e-8, .5e7, .5e7]
		popt, pcov = cf(ia.double_peak, bins[1:], n, p0 = p0)
		sigma = popt[1]
		T = (me * sigma**2 / kb)/11600;
		y = ia.double_peak(bins[1:], *popt)
		l = ax.plot(bins[1:], y, '--r', label = 'T = %.3f' % T + ' eV')
		print(popt[2])
	else:
		(mu, sigma) = norm.fit(v)
		T = (me * sigma**2 / kb)/11600;
		y = mlab.normpdf(bins, mu, sigma)
		l = ax.plot(bins, y, '--r', label = 'T = %.3f' % T + ' eV')
	ax.set_title(dim + '-EDF')
		# Change ticks to energy, normalize y ticks to a max of 1
	xlocs = ax.get_xticks()
	ylocs = ax.get_yticks()
	ylocs = np.round(np.linspace(0,1, len(ylocs)),2)
	ax.set_yticklabels(ylocs)
	xlocs = .5 * me * xlocs**2 * 6.242e18 
	ax.set_xticklabels(np.round(xlocs,2))
	ax.set_xlabel('Energy (eV)')
	ax.set_ylabel('Counts (AU)')
	ax.legend()
	if save:
		plt.savefig(path + name + '.png')
	plt.show()

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



def make_color(weights):
	# Make a color scale by particle weight
	minima = min(weights)
	maxima = max(weights)
	norm = colors.Normalize(vmin = minima, vmax = maxima, clip = True)
	mapper = cm.ScalarMappable(norm = norm, cmap = cm.jet)
	return mapper.to_rgba(weights)

def plot_gas_evolution(gas, pulse, t, LX, LY, savePath, aniName,save = False):
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
	save : boolean, optional
		set True if you want to save the animation
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
					 aspect = 'auto', animated = True, alpha = .75, \
					 cmap = cmap_p, \
					 vmin = np.amin(pulse), vmax = np.amax(pulse))

	plt.xlabel(r'x ($\mu$m)')
	plt.ylabel(r'y ($\mu$m)')
	cbar = plt.colorbar(mappable = im)
	cbar.set_label(r'm$^{-3}$')
	def updatefig(i, gas, pulse):
		im.set_data(gas[i].transpose())
		im2.set_data(pulse[:,:,i])
		plt.title('t = %.2f' %t[i] + ' (fs) \n Neutral gas density')
	ani = animation.FuncAnimation(fig, updatefig, frames = len(gas), \
					fargs = (gas, pulse), interval = 1000)
	if save:
		ani.save(savePath + aniName)

def plot_pulse(pulse, t, LX, LY, savePath, aniName, save = False):
	fig = plt.figure(figsize = (5,3.5), dpi = 150)
	ax  = fig.gca()
	ext = np.array([-LX/2, LX/2, -LY/2, LY/2])
	pulse = pulse[:,:,:,1]#**2
	pulse = pulse.transpose()
	im = plt.imshow(pulse[:,:,0], extent = ext,\
					 aspect = 'auto', animated = True, \
					 vmin = np.amin(pulse), vmax = np.amax(pulse) )
	plt.xlabel(r'x ($\mu$m)')
	plt.ylabel(r'y ($\mu$m)')
	#cbar = plt.colorbar(mappable = im)
	#cbar.set_label(r'm$^{-3}$')
	def updatefig(i, pulse = pulse):
		im.set_array(pulse[:,:,i])
		plt.title('t = %.2f' %t[i] + ' (fs) \n Laser pulse')	
	ani = animation.FuncAnimation(fig, updatefig, frames = len(pulse[0,0,:]),\
					interval = 1000)
	if save:
		ani.save(savePath + aniName)

def plot_phase_space(species, ts, t, i_start,savePath, aniName, save = False):
	# Set dimensions for plot
	fig = plt.figure(figsize = (5,3.5), dpi = 150)
	ax  = fig.gca()
	ax.set_xlabel(r'$\gamma v_x$ (1e6 m/s)') 
	ax.set_ylabel(r'$\gamma v_y$ (1e6 m/s)')
	scatt = ax.scatter([], [], s = .25, marker = 'o')
	ax.set_ylim([-8, 8])
	ax.set_xlim([-8, 8])
	xdata, ydata = [], []
	def update(i, species = species):
		#print(i)
		if i >= i_start:
			xdata = analyze.get_ux(species[i-i_start])
			ydata = analyze.get_uy(species[i-i_start])
			weights = analyze.get_weights(species[i-i_start])
			scatt.set_color(make_color(weights))
			scatt.set_offsets(np.c_[xdata*1e-6,ydata*1e-6])
		ax.set_title(' t = %.2f' % t[i] + ' (fs) \n' + ts)
	ani = animation.FuncAnimation(fig, update, frames = len(t), \
								 interval = 1000)
	if save:
		ani.save(savePath + aniName)

def plot_v_dist(species, ts, t, i_start, savePath, aniName, save = False):
	fig = plt.figure(figsize = (5, 7), dpi = 150)
	ax1  = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	ax1.set_xlabel(r'$\gamma v_x (1e6 m/s)$')
	ax1.set_ylabel('Counts (A.U.)')
	ax2.set_xlabel(r'$\gamma v_y (1e6 m/s)$')
	ax2.set_ylabel('Counts (A.U.)')
	plt.subplots_adjust(hspace = 0.2)
	def update(i, species = species):
		if i >= i_start:
			ax1.clear()
			ax2.clear()
			ux = analyze.get_ux(species[i-i_start])
			uy = analyze.get_uy(species[i-i_start])
			weights = analyze.get_weights(species[i-i_start])

			
			nx, binsx, px = ax1.hist(ux, 200, weights = weights)#, normed = True)
			ny, binsy, py = ax2.hist(uy, 200, weights = weights)#, normed = True)
			
			xlocsx = ax1.get_xticks()
			ax1.set_xticklabels(xlocsx/1e6)
			
			xlocsy = ax2.get_xticks()
			ax2.set_xticklabels(xlocsy/1e6)
			
		ax1.set_title('t = %.2f' % t[i] + '(fs) \n' + ts)
		ax1.set_xlabel(r'$\gamma v_x (1e6 m/s)$')
		ax1.set_ylabel('Counts (A.U.)')
		ax2.set_xlabel(r'$\gamma v_y (1e6 m/s)$')
		ax2.set_ylabel('Counts (A.U.)')
	ani = animation.FuncAnimation(fig, update, frames = len(t), \
								 interval = 1000)
	if save:
		ani.save(savePath + aniName)
def ionization_vis(gas, pulse, electrons, i_start,\
				   t, LX, LY, savePath, \
				   aniName, save = False):
	fig = plt.figure(figsize = (5, 10))
	xl = [-LX/2, LX]
	yl = [-LY/2, LY/2]
	ax1 = fig.add_subplot(211)
	ax1.xaxis.set_ticklabels([])
	ax1.set_xlim(xl)
	ax1.set_ylim(yl)
	ax1.set_ylabel(r'y ($\mu$m)')
	
	ax2 = fig.add_subplot(212)
	ax2.set_xlim(xl)
	ax2.set_ylim(yl)
	ax2.set_ylabel(r'y ($\mu$m)')
	ax2.set_xlabel(r'x ($\mu$m)')

	
	
	gas = np.squeeze(gas)
	ext = np.array([-LX/2, LX/2, -LY/2, LY/2])

	plt.subplots_adjust(hspace = 0.2)
	# y-intensity (for pulse)
	pulse = pulse[:,:,:,1]**2
	pulse = pulse.transpose()
	cmap_p = alpha_colormap(plt.cm.get_cmap('viridis'), .01 * np.amax(pulse))
	im = ax1.imshow(gas[0].transpose(), extent = ext, cmap = 'jet', \
					 aspect = 'auto', animated = True)
	im2 = ax1.imshow(pulse[:,:,0], extent = ext,\
					 aspect = 'auto', animated = True, alpha = .5, \
					 cmap = cmap_p, \
					 vmin = np.amin(pulse), vmax = np.amax(pulse))

	#cbar = plt.colorbar(mappable = im)
	#cbar.set_label(r'm$^{-3}$')
	scatt1 = ax2.scatter([], [], s = 0.25, marker = 'o')
	
	def update(i, gas = gas, pulse = pulse, i_start = i_start, \
				  electrons = electrons ):
		im.set_data(gas[i].transpose())
		im2.set_data(pulse[:,:,i])
		if i >= i_start:
			xdata_e = analyze.get_x(electrons[i-i_start])
			ydata_e = analyze.get_y(electrons[i-i_start])
			weights_e = analyze.get_weights(electrons[i-i_start])
			uy = analyze.get_uy(electrons[i-i_start])
			ux = analyze.get_ux(electrons[i-i_start])
			scatt1.set_color(make_color(weights_e))
			scatt1.set_offsets(np.c_[xdata_e*1e6 - 60, ydata_e*1e6])


		
		
	ani = animation.FuncAnimation(fig, update, frames = len(t), \
								 interval = 1000)
	if save:
		ani.save(savePath + aniName)

def plot_recomb(elec, neutrals, i_start, ts, t, savePath, aniName, \
				save = False):
	fig = plt.figure(figsize = (5,3.5), dpi = 150)
	ax  = fig.gca()
	ax.set_xlabel(r'x ($\mu$m)') 
	ax.set_ylabel(r'y ($\mu$m)')
	scatt = ax.scatter([], [], s = .25, marker = 'o', label = 'Electrons')
	scatt1 = ax.scatter([], [], s = .25, marker = 'o', c = 'r', label = 'Ar')
	ax.set_ylim([-300, 300])
	ax.set_xlim([-300, 300])
	ax.axhline(y= LY/2, color='r', linestyle='-')
	ax.axhline(y= -LY/2, color='r', linestyle='-')
	ax.axvline(x= LX/2, color='r', linestyle='-')
	ax.axvline(x= -LX/2, color='r', linestyle='-')
	xdata, ydata = [], []
	def update(i, elec = elec, neutrals = neutrals):
		xdata = analyze.get_x(elec[i])* 1e6
		ydata = analyze.get_y(elec[i]) * 1e6
		if i >= i_start:
			xdata1 = analyze.get_x(neutrals[i-i_start])*1e6
			ydata1 = analyze.get_y(neutrals[i-i_start])*1e6
			scatt1.set_offsets(np.c_[xdata1, ydata1])
		weights = analyze.get_weights(elec[i])
		scatt.set_color(make_color(weights))
		scatt.set_offsets(np.c_[xdata,ydata])
		ax.set_title(' t = %.2f' % t[i] + ' (ps) \n' + ts)
	ani = animation.FuncAnimation(fig, update, frames = len(elec), \
								 interval = 1000)
	if save:
		ani.save(savePath + aniName)
def plot_weights(elec, neutrals, recomb_start, t):
	w_elec = np.zeros(len(elec))
	w_neutral = np.zeros(len(elec))
	diff = np.zeros(len(elec))
	for i in range(len(elec)):
		we = analyze.get_weights(elec[i])
		w_elec[i] = sum(we)

		if i >= 1:
			diff[i] = w_elec[0] - w_elec[i]
		if i >= recomb_start:
			wn = analyze.get_weights(neutrals[i-recomb_start])
			w_neutral[i] = sum(wn)
	# Normalize these
	w_elec = w_elec/ np.amax(w_elec)
	w_neutral = w_neutral / np.amax(w_neutral)
	diff = diff/np.amax(diff);
	w_diffuse = diff - w_neutral
	fig = plt.figure(figsize = (5,3), dpi = 150);
	ax = fig.gca()
	ax.plot(t*1e9, w_elec, '-b', label = 'Electrons')
	ax.plot(t*1e9, w_neutral, '-r', label = 'Ar')
	ax.plot(t*1e9, w_diffuse, '-g', label = 'Diffused electrons')
	ax.set_xlabel('t (ns)')
	ax.set_ylabel('Sum of weights')
	ax.set_ylim([0, 1])
	plt.legend()
	plt.show()
