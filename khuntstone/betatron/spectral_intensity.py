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


def get_single_intensity(theta, w, x0n, x1n, x2n, v0n, v1n, dtau):
    """
    Function to get the radiated intensity from a single particle for a given
    observation angle and frequency
    
    Parameters:
    -----------
    theta : float
            The observation angle (y-z plane), rad
    w     : float
            The query frequency, rad/s
    x0n   : array_like
            Array of particles 4-position throughout the simulation
    x1n   : array_like
            Array of the linear interpolation coeffient of the trajectory
    x2n   : array_like
            Array of the quadratic interpolation coefficient of the trajectory
    v0n   : array_like
            Array of the particles velocity throughout the simulation
    v1n   : array_like
            Array of the linear interpolation coefficient of the velocity
    dtau  : array_like
            Array of the proper time step throughout the simulation
    """
    #%% Assign observation angle (y-z) and frequency then compute wave-vector

    sx = 0; sy = np.sin(theta); sz = np.cos(theta) # unit observation vector
    kappa = np.array([1, sx, sy, sz]) * (w/c)
    
    # Get reduced chi variables
    chi0n = kappa[0]*x0n[0] - kappa[1]*x0n[1] - kappa[2]*x0n[2] \
            - kappa[3]*x0n[3]
    chi1n = kappa[0]*x1n[0] - kappa[1]*x1n[1] - kappa[2]*x1n[2] \
            - kappa[3]*x1n[3]
    chi2n = kappa[0]*x2n[0] - kappa[1]*x2n[1] - kappa[2]*x2n[2] \
            - kappa[3]*x2n[3]
            
    ind_ex   = np.argwhere(chi2n * dtau > 0.1)
    ind_t    = np.argwhere(chi2n * dtau <= 0.1)
    chi2n[chi2n <= 0] = np.nan
    
    #%% Get reduced variable angles (only where taylor expansion isn't valid)
    
    dtex     = dtau[ind_ex]
    chi1n_ex = chi1n[ind_ex]
    chi2n_ex = chi2n[ind_ex]
    
    
    phi_p  = (dtex**2 / 4) * chi2n_ex + (dtex / 2) * chi1n_ex
    phi_m  = (dtex**2 / 4) * chi2n_ex - (dtex / 2) * chi1n_ex
    
    theta_p = chi1n_ex + chi2n_ex * dtex / (np.sqrt(2 * np.pi * chi2n_ex))
    theta_m = chi1n_ex - chi2n_ex * dtex / (np.sqrt(2 * np.pi * chi2n_ex))
    #%% Get reduced vector psis (only where taylor expansion isn't valid)
    p_factor = np.sqrt(2 * np.pi / chi2n_ex) * np.cos(chi1n_ex**2 \
              / (4 * chi2n_ex))
    m_factor = np.sqrt(2 * np.pi / chi2n_ex) * np.sin(chi1n_ex**2 \
              / (4 * chi2n_ex))
    
    psi_p_x  = p_factor * (2 * chi2n_ex * v0n[0][ind_ex] \
                - chi1n_ex * v1n[0][ind_ex])
    psi_p_y  = p_factor * (2 * chi2n_ex * v0n[1][ind_ex] \
               - chi1n_ex * v1n[1][ind_ex])
    psi_p_z  = p_factor * (2 * chi2n_ex * v0n[2][ind_ex] \
               - chi1n_ex * v1n[2][ind_ex])
    
    psi_m_x  = m_factor * (2 * chi2n_ex * v0n[0][ind_ex] \
               - chi1n_ex * v1n[0][ind_ex])
    psi_m_y  = m_factor * (2 * chi2n_ex * v0n[1][ind_ex] \
               - chi1n_ex * v1n[1][ind_ex])
    psi_m_z  = m_factor * (2 * chi2n_ex * v0n[2][ind_ex] \
               - chi1n_ex * v1n[2][ind_ex])
    
    #%% Get real and imaginary parts of the integral 
    RIx = np.zeros(len(chi2n))
    RIy = np.zeros(len(chi2n))
    RIz = np.zeros(len(chi2n))
    
    ImIx = np.zeros(len(chi2n))
    ImIy = np.zeros(len(chi2n))
    ImIz = np.zeros(len(chi2n))
    
    # Compute useful factors, fresnel functions
    invchi2 = (1 / (4 * chi2n_ex))
    Cp, Sp  = fresnel(theta_p)
    Cm, Sm  = fresnel(theta_m)
    
    RIx[ind_ex] = invchi2 * ((psi_p_x * (Cp - Cm)) + (psi_m_x * (Sp-Sm)) \
          - (2*v1n[0][ind_ex]) * (np.sin(phi_p) - np.sin(phi_m)))
    
    RIy[ind_ex] = invchi2 * ((psi_p_y * (Cp - Cm)) + (psi_m_y * (Sp-Sm)) \
          - (2*v1n[1][ind_ex]) * (np.sin(phi_p) - np.sin(phi_m)))
          
    RIz[ind_ex] = invchi2 * ((psi_p_z * (Cp - Cm)) + (psi_m_z * (Sp-Sm)) \
          - (2*v1n[2][ind_ex]) * (np.sin(phi_p) - np.sin(phi_m)))
    
    ImIx[ind_ex] = invchi2 * ((psi_p_x * (Sp-Sm)) - (psi_m_x * (Cp-Cm)) \
                   - 2*v1n[0][ind_ex] * (np.cos(phi_p) - np.cos(phi_m)))
                   
    ImIy[ind_ex] = invchi2 * ((psi_p_y * (Sp-Sm)) - (psi_m_y * (Cp-Cm)) \
                   - 2*v1n[1][ind_ex] * (np.cos(phi_p) - np.cos(phi_m)))
                   
    ImIz[ind_ex] = invchi2 * ((psi_p_z * (Sp-Sm)) - (psi_m_z * (Cp-Cm)) \
                   - 2*v1n[2][ind_ex] * (np.cos(phi_p) - np.cos(phi_m))) 
                   
    # Imaginary parts have a +/- in them
                   
    #ImIx_m = -1.0 * ImIx
    #ImIy_m = -1.0 * ImIy
    #ImIz_m = -1.0 * ImIz
    #%% Compute using taylor approximation where appropriate
    
    dtt     = dtau[ind_t]
    chi1n_t = chi1n[ind_t]
    chi2n_t = chi2n[ind_t]             
    
    RIx[ind_t] = v0n[0][ind_t] * np.sinc(chi1n_t * dtt / 2) * dtt
    RIy[ind_t] = v0n[1][ind_t] * np.sinc(chi1n_t * dtt / 2) * dtt
    RIz[ind_t] = v0n[2][ind_t] * np.sinc(chi1n_t * dtt / 2) * dtt
    
    
    I1 = (dtt / chi1n_t) * (np.sinc(chi1n_t * dtt / 2) \
          - np.cos(chi1n_t * dtt / 2))
    I2 = (dtt**3 / 4) * np.sinc(chi1n_t * dtt / 2) - (2 / chi1n_t) * I1
    
    ImIx[ind_t] =  v1n[0][ind_t] * I1 + chi2n_t * v0n[0][ind_t] * I2
    ImIy[ind_t] =  v1n[1][ind_t] * I1 + chi2n_t * v0n[1][ind_t] * I2
    ImIz[ind_t] =  v1n[2][ind_t] * I1 + chi2n_t * v0n[2][ind_t] * I2
    #%% Compute spectral intensity
    # Add back the phase factor
    Ix = np.exp(1j*chi0n) * (RIx + 1j*ImIx)
    Iy = np.exp(1j*chi0n) * (RIy + 1j*ImIy)
    Iz = np.exp(1j*chi0n) * (RIz + 1j*ImIz)
    
    # Compute componets of s cross I
    scIx = sy*np.nansum(Iz) - sz*np.nansum(Iy)
    scIy = sz*np.nansum(Ix) - sx*np.nansum(Iz)
    scIz = sx*np.nansum(Iy) - sy*np.nansum(Ix)
    # Take magntiude squared of the cross product
    
    scI2 = scIx * np.conj(scIx) + scIy * np.conj(scIy) + scIz * np.conj(scIz)
    
    # Finally get the radiate intensity
    d2I = mu_0 * e**2 * c * w**2 / (16 * np.pi**3) * scI2
    #print(d2I)
    
    # Test comparison
    #RSx = np.nansum(RIx * np.cos(chi0n) - ImIx*np.sin(chi0n))
    #RSy = np.nansum(RIy * np.cos(chi0n) - ImIy*np.sin(chi0n))
    #RSz = np.nansum(RIz * np.cos(chi0n) - ImIz*np.sin(chi0n))
    
    #ImSx = np.nansum(ImIx * np.cos(chi0n) + RIx * np.sin(chi0n))
    #ImSy = np.nansum(ImIx * np.cos(chi0n) + RIx * np.sin(chi0n))
    #ImSz = np.nansum(ImIx * np.cos(chi0n) + RIx * np.sin(chi0n))
    
    #x_term = RSx**2 + ImSx**2
    #y_term = (RSy*np.cos(theta) - RSz*np.sin(theta))**2
    #z_term = (ImSy*np.cos(theta) - ImSz*np.sin(theta))**2
    #d2I_test = (mu_0 * e**2 * c * w**2 / (16 * np.pi**3)) \
    #           * (x_term + y_term + z_term)
               
    return d2I
