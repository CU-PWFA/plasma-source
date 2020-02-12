#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:33:36 2019

Module containing functions for computing the spectral intensity from a
charged particle with an arbitrary trajectory

@author: keenan
"""

# Python code
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, mu_0
from scipy.special import fresnel
import sys
# Custom modules
sys.path.insert(0, "../")
from v_x_interp import position_interp, velocity_interp
from plotting import makefig
plt.style.use('huntstone')

def get_coeffs(x, v, tau, tau_int, plot = False):
    """
    Function to interpolate position and velocity and return the interpolation
    coefficients
    
    Parameters:
    -----------
    x : array_like
        Array of the particles 4-position (ct, x, y, z) for each time step
    v : array_like
        Array of the particles 4-velocity (c*gamma, vx, vy, vz) for each time
        step
    tau : array_like
        Proper time array corresponding to x and v
    tau_int : array_like
        Array of times to interpolate over
    plot : boolean
        Whether or not to plot the interpolation
        
    Returns:
    --------
    x1n : array_like
        Array of linear interpolation coefficients for the 4-position
    x2n : array_like
        Array of quadratic interpolation coefficients for the 4-position
    v1n : array_like
        Array of linear interpolation coefficients for the 4-velocity
    """
    
    # Interpolate position
    ct1, ct2, ct_int = position_interp(x[0], tau, tau_int)
    x1, x2, x_int    = position_interp(x[1], tau, tau_int)
    y1, y2, y_int    = position_interp(x[2], tau, tau_int)
    z1, z2, z_int    = position_interp(x[3], tau, tau_int)
    
    x4_interp = (ct_int, x_int, y_int, z_int)
    x1n       = (ct1, x1, y1, z1)
    x2n       = (ct2, x2, y2, z2)
    
    # Interpolate velocity
    v1, v_int   = velocity_interp(v[0], tau, tau_int)
    vx1, vx_int = velocity_interp(v[1], tau, tau_int)
    vy1, vy_int = velocity_interp(v[2], tau, tau_int)
    vz1, vz_int = velocity_interp(v[3], tau, tau_int)

    v4_interp = (v_int, vx_int, vy_int, vz_int)
    v1n       = (v1, vx1, vy1, vz1)    
        
    # Plot if requested
    def plot_interp(x, x_int, x1, x2, tau, tau_int, xlab, ylab):
        fig1, ax1 = makefig(xlab = xlab, ylab = ylab)
        ax1.plot(tau, x, label = 'WARGSim')
        ax1.plot(tau_int, x_int, '--', label = 'Inteprolation')
        ax1.legend()
        fig2, ax2 = makefig(xlab = xlab, ylab = 'Linear coefficient')
        ax2.plot(tau, x1)
        
        if isinstance(x2, np.ndarray):
            fig3, ax3 = makefig(xlab = xlab, ylab = 'Quadratic coefficient')
            ax3.plot(tau, x2)
            plt.show()
    if plot:
        ylabs = [r'ct [m]', r'x [$\mu$m]', r'y [$\mu$m]', \
                 r'z [$\mu$m]']
        vylabs = [r'$\gamma$c [m/s]', r'$v_x$ [$\mu$m/s]',\
                  r'$v_y$ [$\mu$m/s]', r'$v_z$ [$\mu$m/s]']
        xlab  = r'$\tau$ [fs]'
        for i in range(len(x4_interp)):
            if i == 0:
                x_use = x[i]
                x_int_use = x4_interp[i]
                v_use = v[i]
                v_int_use = v4_interp[i]
            else:
                x_use = x[i] * 1e6
                x_int_use = x4_interp[i] * 1e6
                v_use = v[i] * 1e6
                v_int_use = v4_interp[i] * 1e6
            plot_interp(x_use, x_int_use, x1n[i], x2n[i], tau*1e15, \
                        tau_int*1e15, xlab, ylabs[i])
            plot_interp(v_use, v_int_use, v1n[i], 0, tau*1e15, \
                        tau_int*1e15, xlab, vylabs[i])
    # Only return momentum part of 4-velocity
    return x1n, x2n, (vx1, vy1, vz1)
    
    
def get_chis(w, s, x1n, x2n, tau, plot = False):
    """
    Function to compute chi_1n and chi_2n the dot product of the 4-wave 
    vector and the 4-position interpolation coefficients
    
    Parameters:
    -----------
    w   : float
          Query frequency in rad/s
    s   : array_like
          unit vector in the direction of observation (sx, sy, sz)
    x1n : array_like
          Array of linear 4-position interpolation coefficients
    x2n : array_like
          Array of quadratic 4-position interpolation coefficients
    tau : array
          Array of proper times corresponding to the x1n and x2n
    plot : boolean, optional
            Whether or not to plot the chis
    
    Returns:
    --------
    chi1n : array_like
             Array of kappa * x1n
    chi2n : array_like
             Array of kappa * x2n
    """
    
    kappa = np.array([1, s[0], s[1], s[2]]) * (w / c)    
    chi1n = kappa[0] * x1n[0] - kappa[1] * x1n[1] - kappa[2] * x1n[2] \
            - kappa[3] * x1n[3]
    chi2n = kappa[0] * x2n[0] - kappa[1] * x2n[1] - kappa[2] * x2n[2] \
            - kappa[3] * x2n[3]
            
    if plot:
        fig1, ax1 = makefig(xlab = r'$\tau$ [fs]', ylab = '$\chi_{1n}$')
        fig2, ax2, = makefig(xlab = r'$\tau$ [fs]', ylab = '$\chi_{2n}$')
        ax1.plot(tau * 1e15, chi1n)
        ax2.plot(tau * 1e15, chi2n)
        plt.show()
    return chi1n, chi2n
    
def get_angles(chi1n, chi2n, dtau):
    """
    Function to get reduced variables theta and phi following Thomas' method
    Parameters:
    -----------
    chi1n : array_like
            Defined above
    chi2n : array_like
            Defined above
    dtau  : float or array of same length as chi2n
            Time-step of the intital particle tracking simulation
    Returns:
    --------
    theta_p : array_like
              Reduced variable theta_plus
    theta_m : array_like
              Reduced variable theta_minus
    phi_p   : array_like
              Reduced variable phi_plus
    phi_m   : array_like
              Reduced variable phi_minus
    """
    theta_p = (chi1n + chi2n * dtau) / np.sqrt(2 * np.pi * chi2n)
    theta_m = (chi1n + chi2n * dtau) / np.sqrt(2 * np.pi * chi2n)
    
    phi_p   = ((dtau**2 / 4) * chi2n) + ((dtau / 2) * chi1n)
    phi_m   = ((dtau**2 / 4) * chi2n) - ((dtau / 2) * chi1n)
    
    return theta_p, theta_m, phi_p, phi_m

def get_psi(chi1n, chi2n, v, v1):
    """
    Function to get reduced variables psi_plus and psi_minus following Thomas'.
    Note both psis are 4-vector quantities
    
    Parameters:
    -----------
    chi1n : array_like
            Defined above
    chi2n : array_like
            Defined above
    v     : array_like
            Array of particles 4-velocity (c * gamma, vx, vy, vz)
    v1    : array_like
            Array of the linear interpolation coefficient for v
    Returns:
    --------
    psi_p : array_like
            Reduced variable psi_plus
    psi_m : array_like
            Reduced variable psi_minus
    """    
    
    # Use functions as psi must be computed for each component of the velocity
    
    def get_velo_term(v0n, v1n):
        return 2 * chi2n * v0n - chi1n * v1n
        
    p_factor = np.sqrt(2 * np.pi / chi2n) * np.cos(chi1n**2 / (4 * chi2n))
    m_factor = np.sqrt(2 * np.pi / chi2n) * np.sin(chi1n**2 / (4 * chi2n))
    
    psi_px   = p_factor * get_velo_term(v[0], v1[0])
    psi_py   = p_factor * get_velo_term(v[1], v1[1])
    psi_pz   = p_factor * get_velo_term(v[2], v1[2])
    
    psi_mx   = m_factor * get_velo_term(v[0], v1[0])
    psi_my   = m_factor * get_velo_term(v[1], v1[1])
    psi_mz   = m_factor * get_velo_term(v[2], v1[2])
    
    return (psi_px, psi_py, psi_pz), (psi_mx, psi_my, psi_mz)
    
def get_I(chi1n, chi2n, dtau, v, v1):
    """
    Function to compute the spectral intensity integral following Thomas
    
    Parameters: all defined above
    -----------
    
    Returns:
    --------
    I : array_like
        Array of the spectral intensity integral values (note: I is a vector)
    """
    
    # Helper functions
    def get_RI_comp(chi2n, psi_p, psi_m, Sp, Sm, Cp, Cm, phi_p, phi_m, v1):
        t1 = psi_p * (Cp - Cm)
        t2 = psi_m * (Sp - Sm)
        t3 = 2 * v1 * (np.sin(phi_p) - np.sin(phi_m))
        return (1 / (4 * chi2n)) * (t1 + t2 + t3)
        
    def get_ImI_comp(chi2n, psi_p, psi_m, Sp, Sm, Cp, Cm, phi_p, phi_m, v1):
        t1 = psi_p * (Sp - Sm)
        t2 = psi_m * (Cp - Cm)
        t3 = 2 * v1 * (np.cos(phi_p) - np.cos(phi_m))
        return (1 / (4 * chi2n)) * (t1 - t2 - t3)
        
    # Search for when to use taylor expansion
    ind    = ~np.isnan(chi2n)
    ind_ex = np.argwhere(dtau[ind] * chi2n[ind] >= 0.1)
    ind_t  = np.argwhere(dtau[ind] * chi2n[ind] < 0.1)
    if isinstance(dtau, np.ndarray):
        dtau_ex = dtau[ind_ex]
        dtau_t  = dtau[ind_t]
    else:
        dtau_ex = dtau
        dtau_t  = dtau
    # Preallocate components of I both real and imaginary
    RIx = np.zeros(len(chi2n))
    RIy = np.zeros(len(chi2n))
    RIz = np.zeros(len(chi2n))
    
    ImIx = np.zeros(len(chi2n))
    ImIy = np.zeros(len(chi2n))
    ImIz = np.zeros(len(chi2n))

    
    theta_p, theta_m, phi_p, phi_m = get_angles(chi1n[ind_ex], chi2n[ind_ex], \
                                                dtau_ex)
                                                
    v_use = (v[0][ind_ex], v[1][ind_ex], v[2][ind_ex])
    v1_use = (v1[0][ind_ex], v1[1][ind_ex], v1[2][ind_ex])
    psi_p, psi_m = get_psi(chi1n[ind_ex], chi2n[ind_ex], v_use, v1_use)
    
    Sp, Cp = fresnel(theta_p)
    Sm, Cm = fresnel(theta_m)
    
    # Get real part where taylor is invalid    
    RIx[ind_ex] = get_RI_comp(chi2n[ind_ex], psi_p[0], psi_m[0], Sp, Sm, Cp,\
                              Cm, phi_p, phi_m, v1_use[0])
                              
    RIy[ind_ex] = get_RI_comp(chi2n[ind_ex], psi_p[1], psi_m[1], Sp, Sm, Cp,\
                              Cm, phi_p, phi_m, v1_use[1])
                              
    RIz[ind_ex] = get_RI_comp(chi2n[ind_ex], psi_p[2], psi_m[2], Sp, Sm, Cp,\
                              Cm, phi_p, phi_m, v1_use[2])
                              
    # Get imaginary part where taylor is invalid               
    ImIx[ind_ex] = get_ImI_comp(chi2n[ind_ex], psi_p[0], psi_m[0], Sp, Sm, Cp,\
                              Cm, phi_p, phi_m, v1_use[0])
                              
    ImIy[ind_ex] = get_ImI_comp(chi2n[ind_ex], psi_p[1], psi_m[1], Sp, Sm, Cp,\
                              Cm, phi_p, phi_m, v1_use[1])
                              
    ImIz[ind_ex] = get_ImI_comp(chi2n[ind_ex], psi_p[2], psi_m[2], Sp, Sm, Cp,\
                              Cm, phi_p, phi_m, v1_use[2])
                              
    # Get real part where taylor is valid
    if any(ind_t):                       
        RIx[ind_t] = v[0][ind_t] * np.sin(chi1n[ind_t] * dtau_t / 2) * dtau_t
        RIy[ind_t] = v[1][ind_t] * np.sin(chi1n[ind_t] * dtau_t / 2) * dtau_t
        RIz[ind_t] = v[2][ind_t] * np.sin(chi1n[ind_t] * dtau_t / 2) * dtau_t
        
        # Get imaginary part where taylor is valid
        I1 = (dtau_t / chi1n[ind_t]) 
        I1 = I1 * (np.sinc(chi1n[ind_t] * dtau_t / 2)\
                   - np.cos(chi1n[ind_t] * dtau_t / 2))
        I2 = (dtau_t**3 / 4) * np.sinc(chi1n[ind_t] * dtau / 2) \
              - (2 / chi1n[ind_t]) * I1
              
              
        ImIx[ind_t] = v1[0][ind_t] * I1 + chi2n[ind_t] * v[0][ind_t] * I2
        ImIy[ind_t] = v1[1][ind_t] * I1 + chi2n[ind_t] * v[1][ind_t] * I2
        ImIz[ind_t] = v1[2][ind_t] * I1 + chi2n[ind_t] * v[2][ind_t] * I2
    
    Ix = RIx + ImIx
    Iy = RIy + ImIy
    Iz = RIz + ImIz
    
    return Ix, Iy, Iz
    

def get_intensity(s, w, x, v, tau, tau_int, dtau):
    """
    Function to compute spectral intensity for a given frequency w
    
    Parameters: defined above
    
    Returns:
    --------
    Intensity : float
                Spectral intensity fordd the given s and w
    """
    
    x1n, x2n, v1n = get_coeffs(x, v, tau, tau_int)
    v             = (v[1], v[2], v[3])
    chi1n, chi2n  = get_chis(w, s, x1n, x2n, tau)
    chi2n[chi2n <= 0]  = np.nan # breaks for negative chi2n
    
    kappa = np.array([1, s[0], s[1], s[2]]) * (w / c)   
    kdx   = kappa[0] * x[0] - kappa[1] * x[1] - kappa[2] * x[2] \
            - kappa[3] * x[3]
    Ix, Iy, Iz  = get_I(chi1n, chi2n, dtau, v, v1n)
    
    Ix = Ix * np.exp(1j * kdx)
    Iy = Iy * np.exp(1j * kdx)
    Iz = Iz * np.exp(1j * kdx)
    # Compute cross product of integral with observation vector
    
    sxIx = s[1] * Iz - s[2] * Iy
    sxIy = s[2] * Ix - s[0] * Iz
    sxIz = s[0] * Iy - s[1] * Ix
    
    sxIx_sum = np.nansum(sxIx)
    sxIy_sum = np.nansum(sxIy)
    szIz_sum = np.nansum(sxIz)
    
    sxI2 = np.abs(sxIx_sum)**2 + np.abs(sxIy_sum)**2 + np.abs(szIz_sum)**2
    
    d2I = sxI2 * mu_0 * e**2 * c * w**2 / (16 * np.pi**3)
    
    return d2I
        
    
    
    


