import Constants.SI as SI
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
global plasmaDict
#plt.style.use('Presentation')
#plt.style.use('fivethirtyeight')
# Dictionary of plasmas and their attributes 
plasmaDict = {'Ar+' : {'Vi' : 15.75962, 'Name' :  'Ar$^{+}$', 'Z' : 1},
		   'Ar2+': {'Vi' : 27.62967, 'Name' : 'Ar$^{2+}$', 'Z' : 2},
		   'Ar3+': {'Vi' : 40.74, 'Name' : 'Ar$^{3+}$', 'Z' : 3},
		   'Ar4+': {'Vi' : 59.81, 'Name': 'Ar$^{4+}$', 'Z' : 4},
		   'Ar5+': {'Vi' : 59.81, 'Name': 'Ar$^{5+}$', 'Z' : 5},
		   'He+': {'Vi' : 59.81, 'Name': 'He$^{+}$', 'Z' : 1},
		   'He2+': {'Vi': 59.81, 'Name': 'He$^{2+}$', 'Z' : 2}
}

def plot_plasma_frac(plasma_frac, pos, beamParams, gasName, ind):
	beta_s = beamParams['beta_s'][ind]
	name = plasmaDict[gasName]['Name']
	title = 'Ionization Fraction of ' + name + ' $\\beta$ = %.2f' % beta_s
	plt.figure(figsize = (5,3), dpi = 150)
	plt.plot(pos['r'][ind]*1e6, plasma_frac[ind])
	plt.xlabel('r [$\mu$m]')
	plt.ylabel('Ionization Fraction')
	plt.title(title)
	plt.show()

def plot_width(widths, plasmaNames, beamParams, logx = False, logy = False, \
			   log = False):
	names = [plasmaDict[i]['Name'] for i in plasmaNames] 
	plt.figure(figsize = (5,3), dpi = 150)
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
def plot_max_frac(max_frac, beamParams, plasmaNames,\
				  logx = False, logy = False, log = False):
	names = [plasmaDict[i]['Name'] for i in plasmaNames]
	plt.figure(figsize = (5,3), dpi = 150)
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
	lg = plt.legend()
	ax = plt.gca()
	plt.show()

def plot_field(field, pos, beamParams, cbar_label, ind = 0, \
			   gas = False, gasName = None, c = SI.lightSpeed):
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
	fig1 = plt.figure(figsize = (5,3), dpi = 150); ax1 = fig1.gca()
	cen = (int(t[-1]/2),0)
	ax1.add_artist(patches.Ellipse(cen, beamParams['sigma_t']*2e15, \
				   beamParams['sigma_r'][ind]*2e6, fc = 'none',\
				   ec = 'w', ls = '--'))
	img1 = ax1.imshow(np.flipud(field[ind]), cmap = 'jet',aspect = 'auto', \
		extent = ext)
	cbar1 = plt.colorbar(mappable = img1, ax = ax1)
	cbar1.set_label(cbar_label)
	
	
	
	ax1.set_xlabel('t [fs]');
	ax1.set_ylabel('r [$\mu$m]');
	ax1.set_title(title)
	
	plt.show()

	# rz
	z = pos['xi'] * 1e6;
	ext = [min(z), max(z), min(r), max(r)]
	fig2 = plt.figure(figsize = (5,3), dpi  = 150)
	ax2 = fig2.gca()
	ax2.add_artist(patches.Ellipse((0,0), beamParams['sigma_z']*2e6, \
				   beamParams['sigma_r'][ind]*2e6, fc = 'none',\
				   ec = 'w', ls = '--'))
	img2 = ax2.imshow(np.flipud(field[ind]), cmap = 'jet', aspect = 'auto', \
		extent = ext)
	cbar2 = plt.colorbar(mappable = img2, ax = ax2)
	cbar2.set_label(cbar_label)
	
	
	ax2.set_xlabel('z [$\mu$m]');
	ax2.set_ylabel('r [$\mu$m]');
	ax2.set_title(title)
	
	plt.show()

def plot_2D_plasma(W, pos, beamParams, gasName, ind = 0, \
	               c = SI.lightSpeed):
	title = gasName + ' plasma'
	r = pos['r'][ind];
	xi = pos['xi']
	t = (xi  / (beamParams['beta']*c)) -\
		(xi[0] /(beamParams['beta']*c)); 
	W_int = np.fliplr(np.cumsum(W[ind], axis = 1)) * ((t[1]-t[0]) * 1e15)
	n_rz = 1 - np.exp(-W_int)
	fig1 = plt.figure(figsize = (5,3), dpi = 150); ax1 = fig1.gca()
	ext = [min(xi)*1e6, max(xi)*1e6, min(r)*1e6, max(r)*1e6]
	img = ax1.imshow(n_rz, cmap = 'jet', extent = ext, aspect = 'auto')
	ax1.add_artist(patches.Ellipse((0,0), beamParams['sigma_z']*2e6, \
				   beamParams['sigma_r'][ind]*2e6, fc = 'none',\
				   ec = 'w', ls = '--'))
	ax1.set_xlabel('z [$\mu$m]')
	ax1.set_ylabel('r [$\mu$m]')
	ax1.set_title(title)
	cbar = plt.colorbar(mappable = img, ax = ax1)
	cbar.ax.tick_params()
	plt.show()
	return n_rz