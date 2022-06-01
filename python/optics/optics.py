# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:37:01 2022

@author: Robert
"""

import numpy as np
from scipy.constants import physical_constants
c = physical_constants['speed of light in vacuum'][0]
eps_0 = physical_constants['vacuum electric permittivity'][0]

def lam_from_omega(omega, n=1.0):
    """ Calculate the wavelength of a light wave from the angular frequency. 
    
    Parameters
    ----------
    omega : float
        Angular frequency of the wave, rad/s.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.
    
    Returns
    -------
    lam : float
        Wavelength of the wave, m.
    """
    return 2*np.pi*c/(omega*n)

def lam_from_f(f, n=1.0):
    """ Calculate the wavelength of a light wave from the frequency.

    Parameters
    ----------
    f : float
        Frequency of the wave, Hz.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.

    Returns
    -------
    lam : float
        Wavelength of the wave, m.
    """
    return c/(f*n)

def omega_from_lam(lam, n=1.0):
    """ Calculate angular frequency of a light wave from the wavelength. 
    
    Parameters
    ----------
    lam : float
        Wavelength of the wave, m.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.
    
    Returns
    -------
    omega : float
        Angular frequency of the wave, rad/s.
    """
    return 2*np.pi*c/(lam*n)

def f_from_lam(f, n=1.0):
    """ Calculate the frequency of a light wave from the wavelength.

    Parameters
    ----------
    lam : float
        Wavelength of the wave, m.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.

    Returns
    -------
    f : float
        Frequency of the wave, Hz.
    """
    return c/(f*n)

def I_from_E(E, n=1.0):
    """ Calculate the laser intensity from the scalar electric field.

    Parameters
    ----------
    E : float or complex
        Electric field of the wave, either amplitude or complex, V/m.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.

    Returns
    -------
    I : float
        Intensity of the optical wave, W/m^2.
    """
    return 0.5*c*n*eps_0*abs(E)**2

def E_from_I(I, n=1.0):
    """ Calculate the amplitude of the electric field from the intensity.

    Parameters
    ----------
    I : float
        Intensity of the optical wave, W/m^2.
    n : float, optional
        Index of refraction of the medium the wave is in, defaults to 1.0.

    Returns
    -------
    E : float
        Electric field amplitude of the wave, V/m.
    """
    return np.sqrt(2*I/(c*n*eps_0))