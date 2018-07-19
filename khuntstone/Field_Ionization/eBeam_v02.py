import sys
sys.path.insert(0, "../")
import Constants.SI as SI
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gm
from scipy.integrate import simps


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
def get_pos(beamParams, nr = 30, nt = 10, nz = 2, npoints = 1000):
    '''
    Quick function to create position r, z, t from beamParams
    '''
    sigma_r, sigma_z, sigma_t = beamParams['sigma_r'], beamParams['sigma_z'],\
                                beamParams['sigma_t']
    r_arr = np.zeros((len(sigma_r), npoints))
    for i in range(len(sigma_r)):
        r_arr[i][:] = np.linspace(-nr * sigma_r[i], nr * sigma_r[i], npoints)
    z_arr = np.linspace(-nz * sigma_z, nz * sigma_z, npoints)
    t_arr = np.linspace(-nt * sigma_t/2, nt*sigma_t/2, npoints)
    return r_arr, z_arr, t_arr
def rad_E_field(pos, beamParams, eps0 = SI.permFreeSpace, \
                c = SI.lightSpeed, peak = False, rz = False):
    '''
    Computes the radial electric field of an electron beam, assumes beam is 
    transversely and radially gaussian.

    Params:
    -------
    pos : dictionary
    dictionary of position and time valus 
    'r' : r - array_like, transvers position m
    'z' : z - array_like, longitudinal position, m
    't' : t - array_like, time, s    
    beamParams : dictionary
        beamParams as described above
    peak : Boolean, optional
        Set to true if you only want the peak E-field 
    rz : Boolean, optional
        Set to true to compute the E(r,z) with t = 0.0, default false

    Returns:
    --------
    Er : array_like
        3D array of electric field values at r, t/z, and sigma_r/beta_s, GV/m
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
        z = pos['z']
        t = pos['t']
        Erz = np.zeros((len(sigma_r), len(r[0]), len(z)))
        Ert = np.zeros((len(sigma_r), len(r[0]), len(t)))
        for i in range(len(sigma_r)):
            rp = np.reshape(r[i], (len(r[i]), 1))
            Erz[i,:,:] = (ppK[i] * sigma_r[i]**2 / (eps0 * rp)) * \
                        (1 - np.exp(-rp**2/(2*sigma_r[i]**2))) * \
                        np.exp(-(z)**2 / (2 * sigma_z**2))
            Erz[i, :,:] = Erz[i,:,:] / 1e9;

            Ert[i,:,:] = (ppK[i] * sigma_r[i]**2 / (eps0 * rp)) * \
                        (1 - np.exp(-rp**2/(2*sigma_r[i]**2))) * \
                        np.exp(-(- beta * c* t)**2 / (2 * sigma_z**2))
            Ert[i, :, :] = Ert[i,:,:] / 1e9;
        return Erz, Ert, rPeak, EPeak
def ionization_rate(Er_z, Er_t, beamParams, Vi, Z = 1):
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
	Vi : float
		the ionization potential of the neutral gas
	Z : int, optional
		atomic residue charge after ionization, default 1

	Returns:
	--------
	W_t : array like
		3D array of ionization rate in sigma_r, r, and t
	W_z : array_like
		3D array of ionization rate in sigma_r, r, and z
	'''
	Vh = 13.6; 
	n = Z / np.sqrt(Vi/Vh);
	Er_z = abs(Er_z)
	Er_t = abs(Er_t)

	W_z = 1.52 * ((4**n * Vi) / (n * gm(2*n))) * \
         (20.5 * Vi**(3/2) / Er_z)**(2*n-1) * \
         np.exp(-6.83 * Vi**1.5 /Er_z);

	W_t = 1.52 * ((4**n * Vi) / (n * gm(2*n))) * \
         (20.5 * Vi**(3/2) / Er_t)**(2*n-1) * \
         np.exp(-6.83 * Vi**1.5 /Er_t);

	return W_z, W_t
def plot_field(field_t, field_z, pos, beamParams, cbar_label, beta_s, ind, \
               lims = [], gas = False, gasName = None):
    '''
    Plots a field in the in the rz and rt planes
    '''
    if gas:
        title = 'Ionization Rate of ' + gasName +'$\\beta$ = %.2f' % beta_s[ind]
    else:
        title = 'Radial Electric Field ' + '$\\beta$ = %.2f' % beta_s[ind] ;
    # rt plane
    r = np.flipud(pos['r'][ind]) * 1e6; nr = len(r)
    t = pos['t'] * 1e15 - pos['t'][0]*1e15; nt = len(t)
    if not lims:
        plt.imshow(np.flipud(field_t[ind]), cmap = 'jet')
    else:
        plt.imshow(np.flipud(field_t[ind]), cmap = 'jet',\
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
    z = pos['z']; nz = len(z)
    sigma_r  = beamParams['sigma_r'][ind]
    sigma_z = beamParams['sigma_z']
    if not lims:
        plt.imshow(np.flipud(field_z[ind]), cmap = 'jet')
    else:
        plt.imshow(np.flipud(field_z[ind]), cmap = 'jet', \
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
