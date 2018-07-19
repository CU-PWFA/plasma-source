# Suite of functions useful for plasma production via field Ionization
# Assume Gaussian electron beam for beamParams
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
    beamParams : array_like, array of beam parameters
        dictionary of beam parameters 
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
def get_sigma_r(en, beta_s, gamma):
    '''
    Computes the width of the electron beam.
    
    Params:
    -------
    en : float
        Beam normalized emittance in m-rad
    beta_s : float
        beam beta function at waist
    gamma : float
        Relativistic factor of the beam
    '''
    return np.sqrt(en * beta_s /gamma)
def rad_E_field(pos, beamParams, eps0 = SI.permFreeSpace,\
                c = SI.lightSpeed, peak = False, rz = False):
    '''
    Computes the radial electric field of a Gaussian eBeam
    
    Params:
    -------
    pos : array_like
    array of position and time valus 
    [0] r - array_like, transvers position m
    [1] t/z - array_like, t in s, z in m
        
    beamParams : array_like
        Array of beam parameters
    peak : Boolean, optional
        Set to true if you only want the peak E-field 

    Returns:
    --------
    Er : array_like
        2D array of electric field values at r and t/z GV/m
    rPeak : float
        Approximate position of peak radial electric field um
    EPeak : float
        Approximate peak radial electric field GV/m
    '''
    sigma_z = beamParams[0]
    sigma_r = beamParams[1]
    beta    = beamParams[2]
    # peak charge density and E-field has no positional dependence
    pPK = peak_charge_dens(beamParams);
    EPeak = (pPK * sigma_r / (2*eps0)) / 1e9;
    rPeak = (np.pi * sigma_r / 2) * 1e6;
    # Preallocate for loop
    if peak:
        return np.nan, rPeak, EPeak
    else:
        if rz:
            r = pos[0]
            z = pos[1]
            t = 0.0
            Er = np.zeros((len(r), len(z)));
            
            for i in range(len(r)):
                for j in range(len(z)):
                    Er[i,j] = (pPK * sigma_r**2 / (eps0 * r[i])) * \
                              (1 - np.exp(-r[i]**2/(2*sigma_r**2))) * \
                              np.exp(-(z[j] - beta*c*t)**2 / (2*sigma_z**2))
                # Put Er in GV/m
                    Er[i,j] = Er[i,j] / 1e9;
        else:
            r = pos[0]
            t = pos[1]
            z = 0.0
            Er = np.zeros((len(r), len(t)));
    
            for i in range(len(r)):
                for j in range(len(t)):
                    Er[i,j] = (pPK * sigma_r**2 / (eps0 * r[i])) * \
                              (1 - np.exp(-r[i]**2/(2*sigma_r**2))) * \
                              np.exp(-(z - beta*c*t[j])**2 / (2*sigma_z**2))
                # Put Er in GV/m
                    Er[i,j] = Er[i,j] / 1e9;
    
        return Er, rPeak, EPeak;
def ionization_rate(Er, beamParams, Vi, Z = 1):
    '''
    Computes the ionization rate of a neutral gas due to the radial electric
    field of a Gaussian electron beam
    
    Params:
    -------
    Er : array_like
        The radial E-field in either the r-t or r-z plane
    Vi : float
        The ionization energy of the neutral gas in eV
    Z : int, optional
        The atomic residue charge (1,2,...) default 1
    Returns:
    --------
    W : array_like
        The ionization rate in the r-t plane in inverse femtoseconds
    '''
    # Constants
    Vh = 13.6;    # eV
    n  = Z / np.sqrt(Vi/Vh); # Effective quantum number
    Er = abs(Er)
    W  = 1.52 * ((4**n * Vi) / (n * gm(2*n))) * \
         (20.5 * Vi**(3/2) / Er)**(2*n-1) * \
         np.exp(-6.83 * Vi**1.5 /Er);
    return W
def maxIonization(beta_s, sigma_z, gamma, en, Q, Vi, c = SI.lightSpeed):
    max_ion = np.zeros(len(beta_s))
    beta = np.sqrt(1 - 1 / gamma**2)
    npoints = 1000
    n_sigma_r = 10;
    n_sigma_t = 10;
    sigma_t = sigma_z / (beta * c)
    t = np.linspace(-(sigma_t/2)*n_sigma_t, (sigma_t/2) * n_sigma_t, npoints)
    for i in range(len(beta_s)):
        sigma_r = get_sigma_r(en, beta_s[i], gamma)
        r = np.linspace(-sigma_r * n_sigma_r, sigma_r * n_sigma_r, npoints)
        pos = [r, t]
        beamParams = [sigma_z, sigma_r, beta, Q]
        dum, dum, Er = rad_E_field(pos, beamParams, peak = True)
        W = ionization_rate(Er, beamParams, Vi)
        n = plasmaDens([W], t, sigma_t)
        max_ion[i] = np.amax(n)
    return max_ion
def plasmaDens(W, t, sigma_t):
    '''
    Calculates the plasma density of the gas ionized by the electron beam in 
    the r-t plane
    
    Parameters:
    ---------- 
    W : array_like
        The ionization rate in the r-t plane
    t : array_like 
        times corresponding to W
    Returns:
    n_plasma : array_like
        The ratio of the ionized plasma to the intial gas density
    '''
    
    # Preallocate for loop
    if len(W) == 1:
        n_plasma = 1 - np.exp(-W[0] * 1e15 * sigma_t)
    else:
        n_plasma = np.zeros((len(W)))
        for i in range(len(W)):
            n_plasma[i] = 1 - np.exp(-simps(W[i,:]*1e15, t))
    return n_plasma
def plot_plasma(r, n_plasma, limits, gas, betas):
    '''
    Plots the ionization fraction of a neutral gas along r. 
    '''
    label1 = r'$\beta$ = ' + str(betas[0])
    plt.plot(r[0] * 1e6, n_plasma[0], '-b', label = label1)
    if len(betas) > 1:
        label2 = r'$\beta$ = ' + str(betas[1])
        plt.plot(r[1] * 1e6, n_plasma[1], '--r', label = label2)
    plt.ylabel('$n/n_0$')
    plt.xlabel('r [$\mu$m]')
    plt.title('Ionization Fraction of ' + gas)
    plt.xlim(limits)
    plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0)
    plt.show()
def plot_field_rt(field, pos, cbar_label, lims = [], gas = False, gasName = None):
    '''
    Plots a field in the in the rz plane
    '''
    if gas:
        title = 'Ionization Rate of ' + gasName
    else:
        title = '';
    r = pos[0] * 1e6; nr = len(r)
    t = pos[1] * 1e15 + pos[1][-1]*1e15; nt = len(t)
    if not lims:
        plt.imshow(field, cmap = 'jet')
    else:
        plt.imshow(field, cmap = 'jet', vmin = lims[0], vmax = lims[1])
    x_locs = [0, nt/2, nt]
    x_labs = [0, 2*int(t[int(len(t)/2 -1)]), 2*int(t[-1])]
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
def plot_field_rz_norm(field, pos, beamParams, cbar_label):
    '''
    Plots a field in the in the rz plane
    '''
    r = pos[0]; nr = len(r)
    z = pos[1]; nz = len(z)
    sigma_r  = beamParams[1]
    sigma_z = beamParams[0]
    plt.imshow(abs(field), cmap = 'jet')
    x_locs = [0, nz/2, nz]
    x_labs = [int(z[0]/sigma_z), int(z[int(len(z)/2 -1)]/sigma_z),\
              int(z[-1]/sigma_z)]
    y_locs = [0, nr/2, nr] 
    y_labs = [int(r[-1]/sigma_r), int(r[int(nr/2 - 1)]/sigma_r), \
              int(r[-0]/sigma_r)]
    plt.xticks(x_locs, x_labs)
    plt.yticks(y_locs, y_labs)
    cbar = plt.colorbar()
    cbar.set_label(cbar_label)
    plt.xlabel('z/$\sigma_z$');
    plt.ylabel('r/$\sigma_r$');
    plt.show()
def plot_2D_plasma(W, r, t, title, lims = [], c = SI.lightSpeed):
    W_int = np.fliplr(np.cumsum(W, axis = 1)) * ((t[1]-t[0]) * 1e15)
    n_rz = 1 - np.exp(-W_int)
    if not lims:
        plt.imshow(n_rz, cmap = 'jet')
    else:
        plt.imshow(n_rz, cmap = 'jet', vmin = lims[0], vmax = lims[1])
    nt = len(t)
    nr = len(r)
    t_arr = t - t[0]
    y_locs = np.array([0, nr/4, nr/2, 3*nr/4, nr-1]); y_locs = [int(i) for i in y_locs]
    x_locs = [0, nt/4, nt/2, 3*nt/4, nt-1]; x_locs = [int(i) for i in x_locs]
    y_labs = [int(i*1e6) for i in r[y_locs]]
    x_labs = [int(i*1e6*c) for i in t_arr[x_locs]]
    plt.yticks(y_locs, y_labs)
    plt.xticks(x_locs, x_labs)
    plt.xlabel('z [$\mu$m]')
    plt.ylabel('r [$\mu$m]')
    plt.title(title)
    plt.colorbar()
    plt.show()
    return n_rz;