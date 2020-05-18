#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing functions for computing the properties of an EO crystal
transcoeff

@author: keenan
"""
from cycler import cycler;
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
eps0 = 8.854187817e-12;


# Colors for plotting.
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77','#CC6677', '#882255',  '#AA4499'];
cy = cycler('color', plot_colors);

def dielec(params):
    '''
    Calculates the complex dielectric function for a range of THz frequencies 
    
    Parameters:
    -----------
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            f     : array_like
                    The range of frequencies to compute the dielectric function
                    in Hz
            plot  : bool, optional
                    Whether or not to plot results
            log   : bool, optional
                    Whether or not to have the xscale of plots be log
    Returns:
    --------
    eps : array_like
            complex dielectric
    n : array_like
        index of refraction (Re[sqrt(eps)])
    kappa : array_like
        attenuation (Im[sqrt(eps)])
    '''
    ctype = params['ctype']; f = params['f'];
    if not bool(params.get('log')):
        log = False;
    else:
        log = params['log'];
    if ctype.lower() == 'gap':
        epsel  = 8.7;
        f0     = 10.98e12; # Hz
        s0     = 1.8;
        gamma0 = 0.02e12; # Hz
    elif ctype.lower() == 'znte':
        epsel  = 7.4;
        f0     = 5.3e12;
        s0     = 2.7;
        gamma0 = 0.09e12; 
    else:
        print("Crystal type not specified")
        epsel  = 1;
        f0     = 0;
        s0     = 0;
        gamma0 = 0;
    
    eps = epsel + s0 * f0**2 / (f0**2 - f**2 - 1j * gamma0 * f);
    n   = np.real(np.sqrt(eps));
    kappa = np.imag(np.sqrt(eps));
    if bool(params.get('plot')):
        if params['plot']:
            fig1 = plt.figure(figsize = (4,4), dpi = 200); ax1 = fig1.gca();
            fig2 = plt.figure(figsize = (4,4), dpi = 200); ax2 = fig2.gca();
            
            ax1.plot(f*1e-12, n);
            ax1.set_xlabel('f [THz]');
            ax1.set_ylabel('n(f)');
            ax1.set_ylim([0, 6]);
            ax1.set_title(ctype + " index of refraction");
            
            ax2.plot(f*1e-12, kappa);
            ax2.set_xlabel('f [THz]');
            ax2.set_ylabel(r'$\kappa$(f)');
            ax2.set_ylim([0, 0.014]);
            ax2.set_title(ctype + " attenuation coefficient");
            if log:
                ax1.set_xscale('log');
                ax2.set_yscale('log');
    return eps, n, kappa

def transcoeff(params):
    '''
    Function to calculate the transmission coefficient as a function of 
    frequency. 
    
    Parameters:
    -----------
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            f     : array_like
                    The range of frequencies to compute the dielectric function
                    in Hz
            plot  : bool, optional
                    Whether or not to plot results
   
    Returns:
    --------
    A : array_like
        The transmission coefficient 
    '''
    
    eps, n, kappa = dielec(params);
    A = 2 / (1 + n + 1j * kappa);
    
    if bool(params.get('plot')):
        if params['plot']:
            fig = plt.figure(figsize = (4,4), dpi = 200)
            ax  = fig.gca();
            ax.plot(params['f']*1e-12, np.real(A));
            ax.set_ylabel('Re[A]')
            ax.set_xlabel('f [THz]')
            ax.set_title(params['ctype'] + " transmission coefficient");
    return A
def eocoeff(params):
    '''
    Parameters:
    -----------
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            f     : array_like
                    The range of frequencies to compute the dielectric function
                    in Hz
            plot  : bool, optional
                    Whether or not to plot results

    Returns:
    ---------
    r41 : array_like
        The electro-optic coefficient
    '''
    if params['ctype'].lower() == 'gap':
        dE     = 1e-12; 
        C      = -0.53;
        f0     = 10.98e12;
        gamma0 = 0.02e12;
    elif params['ctype'].lower() == 'znte':
        dE     = 4.25e-12;
        C      = -0.07;
        f0     = 5.3e12;
        gamma0 = 0.09e12;
    else:
        print('Crystal type not specified');
        dE     = 1;
        C      = 0;
        f0     = 0;
        gamma0 = 0;
        
    r41 = dE * (1 + C * f0**2 \
        / (f0**2 - params['f']**2 - 1j * gamma0 * params['f']));

    if bool(params.get('plot')):
        if params['plot']:
            fig = plt.figure(figsize = (4,4), dpi = 200);
            ax  = fig.gca();
            ax.loglog(params['f']*1e-12, abs(np.real(r41)));
            ax.set_ylim([.01e-12, 40e-12]);
            ax.set_xlabel('f (THz)')
            ax.set_ylabel(r'$r_{41}$')
            ax.set_title(params['ctype'] + " EO coefficient")
            plt.show();
    return r41;
    
def indref(params):
    '''
    Function to calculate the optical index of refrection of a crystal as a 
    function of wavelength. 
    
    Parameters:
    -----------
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            lam   : array_like
                    Wavelength in m
            plot  : bool, optional
                    Whether or not to plot results
    Returns:
    --------
    n : array_like
        index of refraction
    '''
    # If lam is not an array convert it
    lam = params['lam']; ctype = params['ctype'];
    if not (isinstance(lam, list) or type(lam) is np.ndarray):
        lam = np.array([lam]);
    # convert to microns
    lam = lam * 1e6;
    if ctype.lower() == 'gap':
        n = np.sqrt(2.68 + 6.40 * lam**2 / (lam**2 - 0.0903279));
    elif ctype.lower() == 'znte':
        und = np.argwhere(lam < 30);
        ove = np.argwhere(lam >= 30);
        n = np.zeros(len(lam));
        n[und] = np.sqrt(9.92 + 0.42530 / (lam[und]**2 - 0.37766**2) \
                     + 2.63580 / ((lam[und]**2) / (56.5**2) - 1));
        n[ove] = np.sqrt(4.270 + 3.01 * lam[ove]**2 / (lam[ove]**2 - 0.142));
        
    else:
        n = np.zeros(len(lam)) + 1.0;
    if bool(params.get('plot')):
        if params['plot']:
            fig = plt.figure(figsize = (4,4), dpi = 200);
            ax  = fig.gca();
            ax.plot(lam*1e3, n);
            ax.set_xlabel(r'$\lambda$ [nm]');
            ax.set_ylabel(r'n($\lambda$)');
            ax.set_title(ctype + " optical index of refraction");
            plt.show();
    return n
    
def velocities(params):
    '''
    Function to calculate the phase and group velocity, and group velocity 
    dispersion as a function of frequency.
    
    Parameters:
    -----------
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            f   : array_like
                    Frequency in Hz
            plot  : bool, optional
                    Whether or not to plot results
    
    Returns:
    --------
    v_ph : array_like
        Phase velocity
    v_g : array_like
        Group velocity
    g_vd : array_like
        Group velocity dispersion
    '''
    f = params['f'];
    just_one = False;
    if len(f) == 1:
        just_one = True;
        f = np.array([0.999 * f[0], f[0], 1.001 * f[0]])
    # index of refraction
    if all(i < 1e14 for i in f): # THz range
        eps, n, kappa = dielec(params);
    else:
        # optical range
        lam = (c / (f)); # m
        lam[f == 0] = np.nan;
        params['lam'] = lam;
        n   = indref({'ctype' : params['ctype'], 'lam' : lam});
    # get dn/df
    dndf_calc = np.diff(n) / np.diff(f);
    dndf      = np.append(dndf_calc, dndf_calc[-1]);
    
    if just_one:
        n    = n[1];
        f    = f[1];
        dndf = dndf[1];
    v_ph = c / n;
    v_g  = c / (n + f * dndf);
    g_vd = (2 / (c * np.pi)) * dndf * 1e15; # fs^2/mm
    
    if bool(params.get('plot')):
        if params['plot']:
        # Add 800 nm line for comparison 
            lam_opt   = np.array([0.999 * 800, 800, 1.001 * 800]) * 1e-9; # m
            f_opt     = c / lam_opt;
            n_opt     = indref({'ctype': params['ctype'],'lam' : lam_opt});
            dndf_calc = np.diff(n_opt) / np.diff(f_opt);
            dndf      = np.append(dndf_calc, dndf_calc[-1]);
            n_opt     = n_opt[1];
            f_opt     = f_opt[1];
            v_g_opt   = c / (n_opt + f_opt * dndf[1]);
            
            fig = plt.figure(figsize = (4,4), dpi = 200);
            ax  = fig.gca();
            
            ax.plot(f*1e-12, v_g / c, '-b', label = r'$v_g$');
            ax.plot(f*1e-12, v_ph / c, '--r', label = r'$v_{ph}$');
            ax.axhline(v_g_opt / c, linestyle = '-.', color = 'k', \
                       label = r'$v_g$ (800 nm)');
            ax.set_xlabel('f [THz]');
            ax.set_ylabel('v/c');
            ax.set_ylim([0, v_ph[0]/c]);
            ax.set_title(params['ctype'] + " velocities");
            plt.legend();
            plt.show();
    return v_ph, v_g, g_vd  

def georesponse(params):
    '''
    Function to calculate the geometric response function as a functon of
    frequency and crystal thickness
     
    Parameters:
    -----------
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            f   : array_like
                    Frequency in Hz
            d   : array_like
                  Crystal thicknesses in m
            y0  : float
                  Central wavelength of the probe pulse in m
            plot  : bool, optional
                    Whether or not to plot results
    
    Returns:
    --------
    G : array_like
        The crystal geometric response function
    '''
    # Supress extra plots
    if bool(params.get('plot')):
        plot = params['plot'];
        params['plot'] = False;
    else:
        plot = False;

    # THz range
    v_ph, v_g, g_vd = velocities(params);
    
    # Optical
    lam_opt = np.linspace(params['y0'] - 25, params['y0'] + 25, 51) * 1e-9; # m
    f_opt   = (c / lam_opt); # Hz;
    v_ph_opt, v_g_opt, d_vg_opt = velocities({'ctype' : params['ctype'], \
                                              'f' : f_opt});
    v_g_opt = np.mean(v_g_opt);
    
    # Attenuation coeff
    eps, n, kappa = dielec(params);
    a = 2 * np.pi * params['f'] * kappa / c;
    # Set longitudinal distance array
    nz = 1000;
    G  = np.zeros((len(params['f']), len(params['d'])), dtype = 'complex128');
    
    inv_del_v = (1 / v_ph - 1/v_g_opt);
    
    for i in range(len(params['d'])):
        z  = np.linspace(0, params['d'][i], nz);
        dz = abs(z[1]-z[0]);
        for j in range(len(z)):
            G[:, i] += (1 / params['d'][i]) \
                       * np.exp(1j * 2 * np.pi * params['f'] * z[j] \
                                * inv_del_v) * np.exp(-a * z[j]) * dz;
    params['plot'] = plot;
    if bool(params.get('plot')):
        if params['plot']:
            fig = plt.figure(figsize = (6,6), dpi = 200);
            ax  = fig.gca();
            for i in range(len(params['d'])):
                ax.plot(params['f'] * 1e-12, abs(G[:,i]), label = 'd = ' \
                        + str(params['d'][i] * 1e6) \
                        + r'$\mu m$');
            ax.set_xlabel('f (THz)');
            ax.set_ylabel('|G(f)|');
            plt.legend()
            plt.show();
    return G;  
    
def eoresponse(params):
    '''
    Function to compute the electro optic response function as a function of 
    frequency and crystal thickness. 
    
    params : dictionary
        Params should have the following items:
            ctype : str
                    Specifies the EO crystal, either GaP or ZnTe
            f   : array_like
                    Frequency in Hz
            d   : array_like
                  Crystal thicknesses in m
            y0  : float
                  Central wavelength of the probe pulse in m
            plot  : bool, optional
                    Whether or not to plot results
    
    Returns:
    --------
    G_eo : array_like
        The electro-optic response function;
    '''
    # Supress extra plots
    if bool(params.get('plot')):
        print("ok")
        plot = params['plot'];
        params['plot'] = False;
    else:
        plot = False;
    G   = georesponse(params);
    A   = transcoeff(params);
    r41 = eocoeff(params);
    params['plot'] = plot;
    if bool(params.get('plot')) and params['plot']:
        fig = plt.figure(figsize = (6,6), dpi = 200);
        ax  = fig.gca();
        ax.set_prop_cycle(cy);
        ax.set_xlabel('f (THz)')
        ax.set_ylabel(r'$G_{EO}$ [pV/m]');
        ax.set_title(params['ctype'] + " electro-optic coefficient");
    G_eo = np.zeros((len(params['f']), len(params['d'])), dtype = 'complex128');
    for i in range(len(params['d'])):
        G_eo[:, i] = G[:, i] * A * r41;
        if bool(params.get('plot')) and params['plot']:
            ax.plot(params['f'] * 1e-12,abs(G_eo[:, i]) * 1e12, label = 'd = ' \
                    + str(np.round(params['d'][i] * 1e6)) + r'$\mu$m');
    if bool(params.get('plot')) and params['plot']:
        ax.set_ylim([0, G_eo[0, 0] * 1e12]);
        plt.legend();
        plt.show();
    return G_eo

