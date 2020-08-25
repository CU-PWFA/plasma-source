# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:00:06 2020

Module containing functions to compute useful reduced variables for betatron
radiation calculations
@author: keenan
"""


import matplotlib.pyplot as plt
plt.style.use('huntstone')
import numpy as np
from scipy.constants import c, mu_0, epsilon_0, e
from scipy.special import fresnel, erf
import sys
# Avoid creating pyc files
sys.dont_write_bytecode = True
sys.path.insert(0, "../")
from plotting import makefig


def remove_spikes(arr, n_it = 20):
    # Removes spikes from an array
    
    for i in range(n_it):
        arr[arr == np.nanmax(arr)] = np.nan
        arr[arr == np.nanmin(arr)] = np.nan
    return arr
    
def dot_product(vec1, vec2):
    """
    Computes the 4-vector dot product
    
    Parameters:
    -----------
    vec1 : array_like
           1st vector
    vec2 : array_like
           2nd vector
    
    Returns:
    --------
    v1v2 : float
           The 4-vector dot product
    """
    
    return (vec1[0] * vec2[0]) - (vec1[1] * vec2[1]) - (vec1[2] * vec2[2]) \
            - (vec1[3] * vec2[3])

def cross_product(vec1, vec2):
    """
    Computes the cross product
    
    Parameters:
    -----------
    vec1 : array_like
           1st vector
    vec2 : array_like
           2nd vector
    
    Returns:
    --------
    v1v2 : array_like
           The vector product
    """
    
    v1v20 = vec1[1]*vec2[2] - vec1[2]*vec2[1]
    v1v21 = vec1[2]*vec2[0] - vec1[0]*vec2[2]
    v1v22 = vec1[0]*vec2[1] - vec1[1]*vec2[0]
    
    return (v1v20, v1v21, v1v22)

def get_SC(z):
    """
    Function to compute the Fresnel integrals S and C
    Parameters:
    ----------
    z : array_like
        Array where entries are purely imaginary or real (array itself can 
        be a mix)
    """
    ip    = 1j + 1
    im    = 1 - 1j
    arg_p = ip * np.sqrt(np.pi) * z * 0.5
    arg_m = im * np.sqrt(np.pi) * z * 0.5
    S     = (ip/4) * (erf(arg_p) - 1j * erf(arg_m))
    C     = (im/4) * (erf(arg_p) + 1j * erf(arg_m))
    return S, C
def get_chis(kappa, x0n, x1n, x2n, plot = False):
    """
    Function to compute reduced variables chi1n and chi2n  
    (Thomas eqns 18-19)
    
    Parameters:
    -----------
    kappa : array_like
            The 4-wave vector (omega / c * (1, s_hat)) where s_hat is the unit
            observation vector
    x0n   : array_like
            Array of the particles position
    x1n   : array_like
            Array of particles linear position interpolation coefficient
    x2n   : array_like
            Array of particles quadratic position interpolation coefficient
    plot  : bool, optional
            Whether or not to plot, default = False
    Returns:
    --------
    chi_n : array_like
             Reduced variable chi_1n (eqn 18)
    ch2n  : array_like 
             Reduced variable chi_2n (eqn 19)
    """
    
    chi0n = dot_product(kappa, x0n)
    chi1n = dot_product(kappa, x1n)
    chi2n = dot_product(kappa, x2n)
    
    if plot:
        fig = plt.figure(dpi = 200)
        
        ax1 = fig.add_subplot(311)
        ax1.plot(chi0n)
        ax1.set_xticklabels([])
        ax1.set_ylabel(r'$\chi_{0n}$')
        
        ax2 = fig.add_subplot(312)
        ax2.plot(chi1n)
        ax2.set_xticklabels([])
        ax2.set_ylabel(r'$\chi_{1n}$')
        
        ax2 = fig.add_subplot(313)
        ax2.plot(chi2n)
        ax2.set_xlabel('n')
        ax2.set_ylabel(r'$\chi_{2n}$')
    
    return chi0n, chi1n, chi2n
    
def get_phis(chi1n, chi2n, dtau, plot = False):
    """
    Function to get reduced variables Phi_+/- (eqn 17)
    
    Parameters:
    -----------
    chi1n : array_like
            chi1n as computed above
    chi2n : array_like
            chi2n as computed above
    dtau  : array_like
            Proper time step
    
    Returns:
    --------
    phi_p : array_like
            Reduced variable Phi_+
    phi_m : array_like
            Reduced variable Phi_+
    """
    
    Phi_p = (dtau**2 / 4) * chi2n + (dtau / 2) * chi1n
    Phi_m = (dtau**2 / 4) * chi2n - (dtau / 2) * chi1n
    
    if plot:
        fig, ax = makefig(xlab = 'n', ylab = r'$\Phi$ [rad]')
        ax.plot(Phi_p, label = r'$\Phi_{+}$')
        ax.plot(Phi_m, label = r'$\Phi_{-}$')
        ax.legend()
        plt.show()
    return Phi_p, Phi_m
    
def get_psis(chi1n, chi2n, v0n, v1n, plot = False):
    """
    Function to get component of psi vectors as defined in eqns (15-16)
    
    Parameters:
    -----------
    chi1n : array_like
            chi1n as computed above
    chi2n : array_like
            chi2n as computed above
    v0n   : array_like
            1D velocity array along tau
    v1n   : array_like
            1D velocity interpolation coefficient along tau
            
    Returns:
    --------
    psi_p : array_like
            psi as given in eqn 15
    psi_m : array_like
            psi as given in eqn 16
    """
    
    
    psi_base = np.sqrt(2 * np.pi / chi2n, dtype = 'complex128')
    psi_base = psi_base * (2*chi2n * v0n - chi1n * v1n) 
    psi_p    = psi_base * np.cos(chi1n**2 / (4 * chi2n))
    psi_m    = psi_base * np.sin(chi1n**2 / (4 * chi2n))
    
    if plot:
        fig, ax = makefig(xlab = 'n', ylab = r'$\Psi$')
        ax.plot(psi_p, label = r'$\Psi_{+}$')
        ax.plot(psi_m, label = r'$\Psi_{-}$')
        ax.legend()
    return psi_p, psi_m


def get_thetas(chi1n, chi2n, dtau, plot = False):
    """
    Function to get theta_+/m as defined in eqn 14
    
    Parameters:
    -----------
    chi1n : array_like
            As defined above
    chi2n : array_like
            As defined above
    dtau  : float
            As defined above
            
    Returns:
    --------
    theta_p : array_like
              As defined in eqn 14
    theta_m : array_like
              As defined in eqn 14
    """
    
    norm    = np.sqrt(2 * np.pi * chi2n, dtype = 'complex128') 
    theta_p = (chi1n + chi2n * dtau) / norm
    theta_m = (chi1n - chi2n * dtau) / norm

    
    if plot:
        fig, ax = makefig(xlab = 'n', ylab = r'$\Theta$')
        ax.plot(np.real(theta_p), label = r'$Re[\Theta_{+}]$')
        ax.plot(np.real(theta_m), label = r'$Re[\Theta_{-}]$')
        ax.plot(np.imag(theta_p), label = r'$Im[\Theta_{+}]$')
        ax.plot(np.imag(theta_m), label = r'$Im[\Theta_{-}]$')
        ax.legend()
    return theta_p, theta_m

def get_fresnels(theta_p, theta_m, plot = False):
    """
    Function to compute the fresnel integrals of theta_p and theta_m
    
    Parameters:
    -----------
    theta_p : array_like
              theta_p as defined above
    theta_m : array_like
              theta_m as defined above

    Returns:
    --------
    Sp : array_like
         FresnelS of theta_p
    Cp : array_like
         FresnelS of theta_p
    Sm : array_like
         FresnelS of theta_m
    Cm : array_like
         FresnelS of theta_m
    """
    
    
    Sp, Cp = get_SC(theta_p)
    Sm, Cm = get_SC(theta_m)
    
    if plot:
        fig, ax = makefig(xlab = 'n', ylab = 'Fresnel Integral')
        ax.plot(Sp, label = r'S($\theta_{+}$)')
        ax.plot(Cp, label = r'C($\theta_{+}$)')
        ax.plot(Sm, label = r'S($\theta_{-}$)')
        ax.plot(Cm, label = r'C($\theta_{-}$)')
        ax.legend()
        plt.show()
    
    return Sp, Cp, Sm, Cm

def fres_In(chi2n, psi_p, psi_m, Cp, Cm, Sp, Sm, v1n, phi_p, phi_m):
    """
    Function to compute the Larmor integral analytically using Fresnel integrals
    All inputs are defined above. 
    
    Returns:
    --------
    RI  : array_like
          Array of the real part of the integral
    ImI : array_like
          Array of the imaginary part of the integral
    """
    
    norm = 1 / (4 * chi2n)
    
    R1 = psi_p * (Cp - Cm)
    R2 = psi_m * (Sp - Sm)
    R3 = 2 * v1n * (np.sin(phi_p) - np.sin(phi_m))
    
    RI = norm * (R1 + R2 + R3)
    
    I1 = psi_p * (Sp - Sm)
    I2 = psi_m * (Cp - Cm)
    I3 = 2 * v1n * (np.cos(phi_p) - np.cos(phi_m))
    
    ImI  = norm * (I1 - I2 - I3)
    
    return RI, ImI 
    
def tay_In(chi1n, chi2n, v0n, v1n, dtau):
    """
    Function to compute the Larmor integral using a Taylor expansion approach
    All inputs/outputs are defined above
    """
    
    arg = chi1n * dtau * 0.5
    I0  = np.sinc(arg/np.pi) * dtau
    I1  = (dtau / chi1n) * (np.sinc(arg/np.pi) - np.cos(arg))
    I2  = (dtau**3 / 4) * np.sinc(arg/np.pi) - (2 * I1 / chi1n)
    
    RI  = v0n * I0
    ImI = v1n * I1 + chi2n * v0n * I2
    
    return RI, ImI
   
def get_single_d2I(x0n, x1n, x2n, th, v0n, v1n, dtau, w, cut = 1e-3):
    """
    Function to compute a single spectral intensity.
    
    Parameters:
    -----------
    x0n  : array_like
           Array of particle 4-displacement
    x1n  : array_like
           Array of 1st order position interpolation coefficients
    x2n  : array_like
           Array of 2nd order position interpolation coefficients
    v0n  : array_like
           Array of particle velocity
    v1n  : array_like
           Array of 1st order velocity interpolation coefficients
    th   : float
           Observation angle between x and z
    dtau : float
           Proper time step corresponding to particle trajectory
    w    : float
           Frequency
    cut : float, optional
          Point at which to switch between Fresnel/Taylor methods
    
    Returns:
    --------
    d2I : float
          The integrated spectral intensity for the given theta/phi pair
    """
    
    # Geometry and chis
    s    = np.array([np.cos(th), 0, np.sin(th)])
    kappa = w / c * np.array([1, s[0], s[1], s[2]])
    chi0n, chi1n, chi2n = get_chis(kappa, x0n, x1n, x2n)
    
   # Break into Fresnel approximation and Taylor expansion
    # Fresnel
    fi = np.argwhere(abs(chi2n * dtau**2) > cut)
    phi_p, phi_m = get_phis(chi1n[fi], chi2n[fi], dtau)
    
    psi_px, psi_mx = get_psis(chi1n[fi], chi2n[fi], v0n[0][fi], v1n[0][fi])
    psi_py, psi_my = get_psis(chi1n[fi], chi2n[fi], v0n[1][fi], v1n[1][fi])
    psi_pz, psi_mz = get_psis(chi1n[fi], chi2n[fi], v0n[2][fi], v1n[2][fi])
    
    theta_p, theta_m = get_thetas(chi1n[fi], chi2n[fi], dtau)
    
    Sp, Cp, Sm, Cm   = get_fresnels(theta_p, theta_m)
    
    Rfx, Ifx = fres_In(chi2n[fi], psi_px, psi_mx, Cp, Cm, Sp, Sm,\
                          v1n[0][fi], phi_p, phi_m)
    Rfy, Ify = fres_In(chi2n[fi], psi_py, psi_my, Cp, Cm, Sp, Sm,\
                          v1n[1][fi], phi_p, phi_m)
    Rfz, Ifz = fres_In(chi2n[fi], psi_pz, psi_mz, Cp, Cm, Sp, Sm,\
                          v1n[2][fi], phi_p, phi_m) 
    
    # Taylor
    ti = np.argwhere(abs(chi2n * dtau**2) <= cut)
    Rtx, Itx = tay_In(chi1n[ti], chi2n[ti], v0n[0][ti], v1n[0][ti], dtau)
    Rty, Ity = tay_In(chi1n[ti], chi2n[ti], v0n[1][ti], v1n[1][ti], dtau)
    Rtz, Itz = tay_In(chi1n[ti], chi2n[ti], v0n[2][ti], v1n[2][ti], dtau) 
    
    R_test = Rfx + Rfy + Rfz
    I_test = Ifx + Ify + Ifz
    if any(np.imag(R_test)) or any(np.imag(I_test)):
        print("Error: complex values found")
        return
    # Create full arrays
    Rx = np.zeros(len(chi2n)); Ix = np.zeros(len(chi2n));
    Ry = np.zeros(len(chi2n)); Iy = np.zeros(len(chi2n));
    Rz = np.zeros(len(chi2n)); Iz = np.zeros(len(chi2n));
    
    Rx[fi] = np.real(Rfx); Ix[fi] = np.real(Ifx); 
    Ry[fi] = np.real(Rfy); Iy[fi] = np.real(Ify); 
    Rz[fi] = np.real(Rfz); Iz[fi] = np.real(Ifz);
    
    Rx[ti] = Rtx; Ix[ti] = Itx; 
    Ry[ti] = Rty; Iy[ti] = Ity; 
    Rz[ti] = Rtz; Iz[ti] = Itz  
    
    
    # Compute vectors present in eqn 26
    chi0n = dot_product(kappa, x0n)
    c0n   = np.cos(chi0n)
    s0n   = np.sin(chi0n)
    RCIS  = (Rx*c0n - Ix*s0n, Ry*c0n - Iy*s0n, Rz*c0n - Iz * s0n)
    ICRS  = (Ix * c0n + Rx * s0n, Iy * c0n + Ry * s0n, Iz * c0n + Rz * c0n)
    # Compute cross product with observation vector
    sxRCIS = cross_product(s, RCIS)
    sxICRS = cross_product(s, ICRS)
    # Take square of the vectors
    sxRCIS2 = sxRCIS[0]**2 + sxRCIS[1]**2 + sxRCIS[2]**2
    sxICRS2 = sxICRS[0]**2 + sxICRS[1]**2 + sxICRS[2]**2
    # Sum over all timesteps to get d2I
    d2I = np.sum(sxRCIS2) + np.sum(sxICRS2)
    # Add units
    units = (mu_0 * e**2 * c / (16 * np.pi**3)) * w**2
    d2I = units * d2I   
    
    return d2I   
    
def get_d2I(x0n, x1n, x2n, th, v0n, v1n, dtau, w, cut = 1e-3):
    """
    Function to get d2I as an array along multiple values of th and w. Inputs
    are the same as get_single_d2I, but th and w are now arrays. 
    """                                    
    
    d2I = np.zeros((len(w), len(th)))
    for i in range(len(w)):
        if (i+1)%(len(w)/10)==0:
            prog = int((i+1)/len(w) * 100)
            print(str(prog) + "%")
        for j in range(len(th)):
            d2I[i,j] = get_single_d2I(x0n, x1n, x2n, th[j], v0n, \
                                      v1n, dtau, w[i], cut)
    return d2I
