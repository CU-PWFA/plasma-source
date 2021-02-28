# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:43:25 2020

@author: Keenan
"""

"""
Code to reproduce math/results from Alec Thomas' PR A&B 2010 paper. 

Definition of variables:
------------------------

chi0n     : array_like,
            Dot product of 0th order position interpolation coefficient (the 
            particle trajectory) and 4-wave vector     
chi1n     : array_like, 
            Dot product of 1st order position interpolation coefficient and 
            4-wave vector (Thomas Eq. 18). Array is along time, not particle 
            number.      
chi2n     : array_like
            Dot product of 2nd order position interpolation coefficient and
            4-wave vector (THomas Eq. 19). Array is along time, not particle 
            number.     
dtau      : float
            Non-interpolated proper time step.   
ImIn      : array_like
            imaginary part of spectral intensity integral given by Thomas 
            Eq. 13 or Eq. 22
ISx/y/z   : array_like
            Combination of integrals defined by Thomas Eq. 28
phi_p/m   : array_like
            Reduced variables phi_+/- from Thomas Eq. 17
psi_p/m   : array_like
            Reduced variables psi_+/- from Thomas Eqs. 15 and 16    
ReIn      : array_like
            real part of spectral intensity integral given by Thomas Eq. 12 or
            Eq. 21
RSx/y/z   : Array_like
            Combination of integrals defined by Thomax Eq. 27
s_hat     : array_like
            Unit vector of observation
sign      : int
            Sign of current frequency value
theta     : float
            Observation angle with respect to the z-axis in the y-z plane       
theta_p/m : array_like
            Reduced variables theta_+/- from Thomas Eq. 14    
v4        : array_like
            Particle 4 velocity array along time.         
v1n       : array_like
            1st order velocity interpolation coefficient array along time. 
"""

import numpy as np
from scipy.constants import c, e, mu_0
from scipy.special import erf
import matplotlib.pyplot as plt


def get_fresnel(z):
    """
    Function to compute fresnel integrals. Needed as scipy does not handle 
    imaginary input well. 
    
    Parameters:
    -----------
    z : float, or array_like
        Input into the fresnel integrals
        
    Returns:
    --------
    C : same type as z,
        The Fresnel integral C
    
    S : same type as z
        The Fresnel integral S
    """
    
    ip    = 1j + 1
    im    = 1 - 1j
    arg_p = ip * np.sqrt(np.pi) * z * 0.5
    arg_m = im * np.sqrt(np.pi) * z * 0.5
    S     = (ip/4) * (erf(arg_p) - 1j * erf(arg_m))
    C     = (im/4) * (erf(arg_p) + 1j * erf(arg_m))
    
    return S, C


def get_phis(chi1n, chi2n, dtau):
    """
    Function to compute reduced variables phi pluse/minus
    """
    phi_p = (0.25 * dtau**2) * chi2n + (0.5 * dtau) * chi1n
    phi_m = (0.25 * dtau**2) * chi2n - (0.5 * dtau) * chi1n
    
    return phi_p, phi_m

def get_psis(chi1n, chi2n, v0n, v1n):
    """
    Function to compute reduced variables Psi plus/minus
    """
    psi_p = np.zeros(len(chi2n), dtype = "complex128")
    psi_m = np.zeros(len(chi2n), dtype = "complex128")
    psi_base = np.zeros(len(chi2n), dtype = "complex128")
    psi_base += np.sqrt(2 * np.pi / abs(chi2n)) \
               * (2 * chi2n * v0n - chi1n * v1n)
    psi_base[chi2n < 0] = 1j * psi_base[chi2n <0]
    psi_p += psi_base * np.cos(0.25 * chi1n**2 / chi2n)
    psi_m +=  psi_base * np.sin(0.25 * chi1n**2 / chi2n)
    
    return psi_p, psi_m

def get_thetas(chi1n, chi2n, dtau):
    """
    Function to compute the reduced variables theta plus/minus 
    """
    
    theta_p = np.zeros(len(chi2n), dtype = "complex128")
    theta_m = np.zeros(len(chi2n), dtype = "complex128")
    theta_p += (chi1n + chi2n * dtau) / np.sqrt(2 * np.pi * abs(chi2n))
    theta_m += (chi1n - chi2n * dtau) / np.sqrt(2 * np.pi * abs(chi2n))
    
    theta_p[chi2n < 0] = 1j * theta_p[chi2n < 0]
    theta_m[chi2n < 0] = 1j * theta_m[chi2n < 0]
    return theta_p, theta_m

def fresnel_I(chi1n, chi2n, v0n, v1n, dtau, sign, plot = False):
    
    """
    Function to compute the spectral intensity integral (Eqs. 12 & 13) using
    Fresnel integrals. NOTE: This computes one vector component for one 
    particle. Must call for each vector component with a different velocity 
    input and loop over all particles. 
    """
    
    # Preallocate for complex values
    ReIn = np.zeros(len(chi2n), dtype = "complex128")
    ImIn = np.zeros(len(chi2n), dtype = "complex128")
    
    # Compute reduced variables
    phi_p, phi_m = get_phis(chi1n, chi2n, dtau)
    
    psi_p, psi_m = get_psis(chi1n, chi2n, v0n, v1n)
    
    theta_p, theta_m = get_thetas(chi1n, chi2n, dtau)
    
    
    Sp, Cp = get_fresnel(theta_p)
    Sm, Cm  = get_fresnel(theta_m)
    
    sin_p = np.sin(phi_p)
    sin_m = np.sin(phi_m)
    cos_p = np.cos(phi_p)
    cos_m = np.cos(phi_m)
    
    ReIn += psi_p * (Cp - Cm) + psi_m * (Sp - Sm) + 2 * v1n * (sin_p - sin_m)
    ReIn = (0.25 / chi2n) * ReIn
    
    ImIn += psi_p * (Sp - Sm) - psi_m * (Cp - Cm) - 2 * v1n * (cos_p - cos_m)
    ImIn = sign * (0.25 / chi2n) * ImIn
    
    if plot:
        # Phis
        fig = plt.figure(figsize = (8, 6))
        ax1  = fig.add_subplot(121)
        ax2  = fig.add_subplot(122)
        ax1.plot(phi_p, label = "phi +")
        ax1.legend()
        ax2.plot(phi_m, label = "phi -")
        ax2.legend()
        plt.show()
        # psis (real)
        fig = plt.figure(figsize = (8, 6))
        ax1  = fig.add_subplot(121)
        ax2  = fig.add_subplot(122)
        ax1.plot(np.real(psi_p), label = "Re psi +")
        ax1.legend()
        ax2.plot(np.real(psi_m), label = "Re psi -")
        ax2.legend()
        plt.show()
        # psi
        fig = plt.figure(figsize = (8, 6))
        ax1  = fig.add_subplot(121)
        ax2  = fig.add_subplot(122)
        ax1.plot(np.imag(psi_p), label = "Im psi +")
        ax1.legend()
        ax2.plot(np.imag(psi_m), label = "Im psi -")
        ax2.legend()
        plt.show()
        # Thetas (real)
        fig = plt.figure(figsize = (8, 6))
        ax1  = fig.add_subplot(121)
        ax2  = fig.add_subplot(122)
        ax1.plot(np.real(theta_p), label = "Re Theta +")
        ax1.legend()
        ax2.plot(np.real(theta_m), label = "Re Theta -")
        ax2.legend()
        plt.show()
        # Thetas (imag)
        fig = plt.figure(figsize = (8, 6))
        ax1  = fig.add_subplot(121)
        ax2  = fig.add_subplot(122)
        ax1.plot(np.imag(theta_p), label = "Im Theta +")
        ax1.legend()
        ax2.plot(np.imag(theta_m), label = "Im Theta -")
        ax2.legend()
        plt.show()
        
        #Re and Im
        # psi
        fig = plt.figure(figsize = (8, 6))
        ax1  = fig.add_subplot(121)
        ax2  = fig.add_subplot(122)
        ax1.plot(ReIn, label = "Real I")
        ax1.legend()
        ax2.plot(ImIn, label = "Imag I")
        ax2.legend()
        plt.show()
        
    return ReIn, ImIn

def taylor_I(chi1n, chi2n, v0n, v1n, dtau, sign):
    """
    Function to compute the spectral intensity integral (Eqs. 21 & 22), using
    a Taylor expansion method. NOTE: This computes one vector component for 
    one particle. Must call for each vector component with a different 
    velocity input and loop over all particles. 
    """
    
    def sinc(x):
        return np.sin(x) / x
    
    # Computed reduced variables
    I0n = sinc(0.5 * chi1n * dtau) * dtau
    I1n = (dtau / chi1n) \
           * (sinc(0.5 * chi1n * dtau) - np.cos(0.5 * chi1n * dtau))
    I2n = (0.25 * dtau**3) * sinc(0.5 * chi1n * dtau) - (2 / chi1n) * I1n
    
    ReIn = v0n * I0n
    ImIn = v1n * I1n + sign * (chi2n * v0n * I2n)
    
    return ReIn, ImIn

def get_chis(w, s_hat, x4, x41n, x42n):
    """
    Function to compute chi0n, chi1n, and chi2n
    """
    kappa = (w/c) * np.array([1, s_hat[0], s_hat[1], s_hat[2]])
    chi0n = -kappa[0] * x4[0] + kappa[1] * x4[1] + kappa[2] * x4[2]\
            + kappa[3] * x4[3]
    chi1n = -kappa[0] * x41n[0] + kappa[1] * x41n[1] + kappa[2] * x41n[2] \
            + kappa[3] * x41n[3]
    chi2n = -kappa[0] * x42n[0] + kappa[1] * x42n[1] + kappa[2] * x42n[2] \
            + kappa[3] * x42n[3]
    chi0n = np.longdouble(chi0n)
    chi1n = np.longdouble(chi1n)
    chi2n = np.longdouble(chi2n)
    
    return chi0n, chi1n, chi2n

def get_d2I(w, s_hat, x4, x41n, x42n, v4, v1n, dtau, theta):
    """
    Function to compute the spectral intensity from Thomas eq. (27)
    """
    
    # Looping for now, speedup later
    d2I = np.zeros(len(w))

    
    for i in range(len(w)):
        om = w[i]
        sign = np.sign(om)
        chi0n, chi1n, chi2n = get_chis(om, s_hat, x4, x41n, x42n)
        
        # Preallocate
        RIx = np.zeros(len(chi2n), dtype = "complex128")
        RIy = np.zeros(len(chi2n), dtype = "complex128")
        RIz = np.zeros(len(chi2n), dtype = "complex128")
        
        IIx = np.zeros(len(chi2n), dtype = "complex128")
        IIy = np.zeros(len(chi2n), dtype = "complex128")
        IIz = np.zeros(len(chi2n), dtype = "complex128")
        
        # Index where to use Fresnel solution or Taylor approximate
        f_ind = np.squeeze(np.argwhere(abs(chi2n * dtau**2) > 1e-3))
        t_ind = np.squeeze(np.argwhere(abs(chi2n * dtau**2) <= 1e-3))
        
        # Compute Fresnel solution
        fresRIx, fresIIx = fresnel_I(chi1n[f_ind], chi2n[f_ind], \
                                     v4[1][f_ind], v1n[0][f_ind], dtau, sign)
        fresRIy, fresIIy = fresnel_I(chi1n[f_ind], chi2n[f_ind], \
                                     v4[2][f_ind], v1n[1][f_ind], dtau, sign)
        fresRIz, fresIIz = fresnel_I(chi1n[f_ind], chi2n[f_ind], \
                                     v4[3][f_ind], v1n[2][f_ind], dtau, sign)
        
        RIx[f_ind] = fresRIx
        RIy[f_ind] = fresRIy
        RIz[f_ind] = fresRIz
        IIx[f_ind] = fresIIx
        IIy[f_ind] = fresIIy
        IIz[f_ind] = fresIIz
        # Compute Taylor approximation
        tayRIx, tayIIx = taylor_I(chi1n[t_ind], chi2n[t_ind], v4[1][t_ind],
                                     v1n[0][t_ind], dtau, sign)
        tayRIy, tayIIy = taylor_I(chi1n[t_ind], chi2n[t_ind], v4[2][t_ind],
                                     v1n[1][t_ind], dtau, sign)
        tayRIz, tayIIz = taylor_I(chi1n[t_ind], chi2n[t_ind], v4[3][t_ind],
                                     v1n[2][t_ind], dtau, sign)
        
        RIx[t_ind] = tayRIx
        RIy[t_ind] = tayRIy
        RIz[t_ind] = tayRIz
        IIx[t_ind] = tayIIx
        IIy[t_ind] = tayIIy
        IIz[t_ind] = tayIIz
        
        # Compute d2I (Thomas Equations 27-29)
        RSx = np.sum(RIx * np.cos(chi0n) - IIx * np.sin(chi0n))
        RSy = np.sum(RIy * np.cos(chi0n) - IIy * np.sin(chi0n))
        RSz = np.sum(RIz * np.cos(chi0n) - IIz * np.sin(chi0n))
        
        ISx = np.sum(IIx * np.cos(chi0n) + RIx * np.sin(chi0n))
        ISy = np.sum(IIy * np.cos(chi0n) + RIy * np.sin(chi0n))
        ISz = np.sum(IIz * np.cos(chi0n) + RIz * np.sin(chi0n))
        
        d2I[i] = RSx**2 + ISx**2 \
             + (RSy * np.cos(theta) - RSz * np.sin(theta))**2 \
             + (ISy * np.cos(theta) - ISz * np.sin(theta))**2
    d2I = d2I * (mu_0 * e**2 * (w)**2) / (16 * np.pi**3)  
        
    return d2I