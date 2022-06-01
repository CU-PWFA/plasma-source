# -*- coding: utf-8 -*-
"""
Created on Thu May  5 12:51:10 2022

@author: Robert
"""

import numpy as np
from scipy.special import gamma
from scipy.constants import physical_constants
c = physical_constants['speed of light in vacuum'][0]
eps_0 = physical_constants['vacuum electric permittivity'][0]

def E_from_energy(energy, tau, w_0, m, n=1.0):
    """ Calculate the peak electric field strength of a super-Gaussian pulse
    
    Parameters
    ----------
    energy : float
        Energy of the pulse, J.
    tau : float
        FWHM pulse length of the pulse, s.
    w_0 : float
        Spot size of the pulse, m.
    m : int
        Super-Gaussian order of the transverse pulse shape, must be even.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.
    
    Returns
    -------
    E_0 : float
        The peak electric field strength of the pulse.
    """
    return np.sqrt(energy*2*np.sqrt(np.log(2))*m/(c*n*eps_0*np.pi**1.5*4**(-1/m)*gamma(2/m)*tau*w_0**2))

def fluence_from_intensity(I_0, tau):
    """ Calculate the laser fluence of a Gaussian pulse from the peak intensity.
    
    Parameters
    ----------
    I_0 : float
        Peak intensity of the laser pulse, W/m^2.
    tau : float
        FWHM pulse length of the pulse, s.
    
    Returns
    -------
    fluence : float
        Laser fluence of the pulse, J/m^2.
    """
    return I_0*tau*np.sqrt(np.pi/(4*np.log(2)))