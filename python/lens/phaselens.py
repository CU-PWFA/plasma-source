#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:43:59 2020

@author: valentina_lee
"""

import numpy as np
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
import pyximport; pyximport.install()
from beam.beams import laserbeam
from ionization import ionization
import matplotlib.pyplot as plt

def phase_lens_1d(X, Ein, F, Zin, Zout, lam):
    """ Generates the intensity profile that is an initial E field goes through a phase only lens at distance Z after the lens

    Parameters
    ----------
    X : double
        Physical width of the grid in the x direction, the grid goes from [-X/2, X/2].
    Ein : array-like, 1d
        E field at the lens
    F : double
        focal lenth of the length
    Zin : double
        The distance from the input E field to the phase mask
    Zout : double
        Intensity at Zout where the intensity output will be. The distance from phase mask to output screen.
    lam : double
        The vacuum wavelength of the laser radiation. in micron

    Returns
    -------
    Eout : array-like, 1d
        Out put Efield 
    """
    k0= 2*np.pi/lam
    
    kx=fftfreq(Ein.shape[0], X/(Ein.shape[0]))*2*np.pi
    kz= np.sqrt(k0**2-kx**2);
    
    E0=ifft(fft(Ein)*np.exp(-1j*kz*Zin))
    
    xreal=np.linspace(-X/2, X/2, Ein.shape[0])
    phase_mask= np.exp(1j*k0*(xreal**2)/F/2)
    
    E_pm= E0*phase_mask
    Eout= ifft(fft(E_pm)*np.exp(-1j*kz*Zout))
    return Eout

'''
def phase_lens_2d(X, Y, Ein, F, Zin, Zout, lam):
    """ Generates the intensity profile that is an initial E field goes through a phase only lens at distance Z after the lens

    Parameters
    ----------
    X : double
        Physical width of the grid in the x direction, the grid goes from [-X/2, X/2].
    Y : double
        Physical width of the grid in the y direction, the grid goes from [-Y/2, Y/2].
    Ein : array-like, 2d [x][y]
        E field at the lens
    F : double
        focal lenth of the length
    Zin : double
        The distance from the input E field to the phase mask
    Zout : double
        Intensity at Zout where the intensity output will be. The distance from phase mask to output screen.
    lam : double
        The vacuum wavelength of the laser radiation. in micron

    Returns
    -------
    Eout : array-like, 2d
        Out put Efield 
    """

    k0= 2*np.pi/lam
    
    kx=fftfreq(4*Ein.shape[1], X/(4*Ein.shape[1]))*2*np.pi
    ky=fftfreq(4*Ein.shape[0], Y/(4*Ein.shape[0]))*2*np.pi
    kX, kY= np.meshgrid(kx, ky);
    kz= np.sqrt(k0**2-kX**2-kY**2);
    
    E= np.zeros((Ein.shape[0]*4, Ein.shape[1]*4), dtype=complex)
    E[int(E.shape[0]*3/8):int(E.shape[0]*5/8), int(E.shape[1]*3/8):int(E.shape[1]*5/8)]= Ein
    
    E0=ifft2(fft2(E)*np.exp(-1j*kz*Zin))

    xreal=np.linspace(-X, X, Ein.shape[1]*4)
    yreal=np.linspace(-Y, Y, Ein.shape[0]*4)
    Xreal, Yreal= np.meshgrid(xreal, yreal)                 
    phase_mask= np.exp(1j*k0*(Xreal**2+Yreal**2)/F/2/8)
    
    E_pm= E0*phase_mask
    Eout= ifft2(fft2(E_pm)*np.exp(-1j*kz*Zout))
    
    Eout= Eout[int(E.shape[0]*3/8):int(E.shape[0]*5/8), int(E.shape[1]*3/8):int(E.shape[1]*5/8)]

    return Eout
'''
def phase_lens_2d(X, Y, Ein, F, Zin, Zout, lam):
    """ Generates the intensity profile that is an initial E field goes through a phase only lens at distance Z after the lens

    Parameters
    ----------
    X : double
        Physical width of the grid in the x direction, the grid goes from [-X/2, X/2].
    Y : double
        Physical width of the grid in the y direction, the grid goes from [-Y/2, Y/2].
    Ein : array-like, 2d [x][y]
        E field at the lens
    F : double
        focal lenth of the length
    Zin : double
        The distance from the input E field to the phase mask
    Zout : double
        Intensity at Zout where the intensity output will be. The distance from phase mask to output screen.
    lam : double
        The vacuum wavelength of the laser radiation. in micron

    Returns
    -------
    Eout : array-like, 2d
        Out put Efield 
    """
    k0= 2*np.pi/lam
    
    kx=fftfreq(2*Ein.shape[1], X/(2*Ein.shape[1]))*2*np.pi
    ky=fftfreq(2*Ein.shape[0], Y/(2*Ein.shape[0]))*2*np.pi
    kX, kY= np.meshgrid(kx, ky);
    kz= np.sqrt(k0**2-kX**2-kY**2);
    
    E= np.zeros((Ein.shape[0]*2, Ein.shape[1]*2), dtype=complex)
    E[int(E.shape[0]/4):int(E.shape[0]*3/4), int(E.shape[1]/4):int(E.shape[1]*3/4)]= Ein
    
    E0=ifft2(fft2(E)*np.exp(-1j*kz*Zin))

    xreal=np.linspace(-X, X, Ein.shape[1]*2)
    yreal=np.linspace(-Y, Y, Ein.shape[0]*2)
    Xreal, Yreal= np.meshgrid(xreal, yreal)                 
    phase_mask= np.exp(1j*k0*(Xreal**2+Yreal**2)/F/2/4)
    
    E_pm= E0*phase_mask
    Eout= ifft2(fft2(E_pm)*np.exp(-1j*kz*Zout))
    
    Eout= Eout[int(E.shape[0]/4):int(E.shape[0]*3/4), int(E.shape[1]/4):int(E.shape[1]*3/4)]
#    Iout= abs(Eout)**2
    return Eout


