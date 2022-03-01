# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:43:25 2020

@author: Keenan
"""

"""
Code to reproduce math/results from Alec Thomas' PR A&B 2010 paper. 
"""
# Imports and constants
import numpy as np
from scipy.constants import c, e, mu_0
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.special import erf

# Other constants
inv_c        = 1/c
inv_pi       = 1/np.pi
sqrt_2pi     = np.sqrt(2*np.pi)
inv_sqrt_2pi = 1 / sqrt_2pi

# Math functions
def get_fresnel(z):
    """
    Function to compute fresnel integrals. Needed as scipy does not handle 
    imaginary input well. 
    
    Parameters:
    -----------
    z : complex float or array *scipy.special.erf does not support complex256*
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

def sincu(x):
    """
    Unnormalized sinc function
    """
    return np.sin(x)/x

# Interpolation coeff function
def get_dtau(tau, var):
    """
    Returns 1st and 2nd order derivatives w.r.t. proper time
    """
    s = ius(tau, var)
    xp = s.derivative(n=1)
    xpp = s.derivative(n=2)
    return xp(tau), xpp(tau)

"""
Below are functions for computing many of the equations/variables in Thomas'
paper. Note they are not called for get_d2I as that nearly doubles computation
time. They are however useful for debugging. 
"""

def get_chis(w, s_hat, x4, x41, x42):
    """
    Function to compute reduced chi variables
    """
    kappa = w * inv_c * np.array([1, s_hat[0], s_hat[1], s_hat[2]])
    chi0  = -kappa[0]*x4[0]+kappa[1]*x4[1]+kappa[2]*x4[2]+kappa[3]*x4[3]
    chi1  = -kappa[0]*x41[0]+kappa[1]*x41[1]+kappa[2]*x41[2]+kappa[3]*x42[3]
    chi2  = -kappa[0]*x42[0]+kappa[1]*x42[1]+kappa[2]*x42[2]+kappa[3]*x42[3]
    inv_chi1         = 1/chi1
    inv_chi2         = 1/chi2
    sqrt_chi2        = np.complex128(np.sqrt(abs(chi2)))
    z_ind            = np.argwhere(chi2<0)
    sqrt_chi2[z_ind] = 1j*sqrt_chi2[z_ind]
    inv_sqrt_chi2    = 1/sqrt_chi2
    return chi0, chi1, chi2, inv_chi1, inv_chi2, inv_sqrt_chi2

def get_phis(chi1, chi2, dtau):
    """
    Function to compute reduced variables phi plus/minus
    """
    phi_p = (0.25 * dtau*dtau) * chi2 + (0.5 * dtau) * chi1
    phi_m = (0.25 * dtau*dtau) * chi2 - (0.5 * dtau) * chi1
    
    return phi_p, phi_m

def get_psis(chi1, chi2, inv_chi2, inv_sqrt_chi2, v0, v1):
    """
    Function to compute reduced variables Psi plus/minus
    """
    psi_base = sqrt_2pi*inv_sqrt_chi2*(2*chi2*v0-chi1*v1)
    psi_p    = psi_base*np.cos(0.25*chi1*chi1*inv_chi2)
    psi_m    = psi_base*np.sin(0.25*chi1*chi1*inv_chi2)
    
    return psi_p, psi_m

def get_thetas(chi1, chi2, inv_sqrt_chi2, dtau):
    """
    Function to compute the reduced variables theta plus/minus 
    """
    theta_p = inv_sqrt_2pi*inv_sqrt_chi2*(chi1+chi2*dtau)
    theta_m  = inv_sqrt_2pi*inv_sqrt_chi2*(chi1-chi2*dtau)
    
    return theta_p, theta_m
    

def fresnel_I(chi1, chi2, inv_chi2, inv_sqrt_chi2, v0, v1, dtau, sign):
    
    """
    Function to compute the spectral intensity integral (Eqs. 12 & 13) using
    Fresnel integrals.
    """
    
    # Compute reduced variables
    phi_p, phi_m = get_phis(chi1, chi2, dtau)
    
    psi_p, psi_m = get_psis(chi1, chi2, inv_chi2, inv_sqrt_chi2, v0, v1)
    
    theta_p, theta_m = get_thetas(chi1, chi2, inv_sqrt_chi2, dtau)
    
    
    Sp, Cp = get_fresnel(np.complex128(theta_p))
    Sm, Cm  = get_fresnel(np.complex128(theta_m))
    
    sin_p = np.sin(phi_p)
    sin_m = np.sin(phi_m)
    cos_p = np.cos(phi_p)
    cos_m = np.cos(phi_m)
    
    ReIn = (0.25*inv_chi2) * (psi_p*(Cp-Cm) \
                              + psi_m*(Sp-Sm) \
                              +2*v1*(sin_p-sin_m))
    
    ImIn = (0.25*inv_chi2)*(psi_p*(Sp-Sm)\
                            -psi_m*(Cp-Cm)\
                            -2*v1*cos_p-cos_m)
    
    
        
    return ReIn, ImIn

def taylor_I(chi1, chi2, inv_chi1, v0, v1, dtau, sign):
    """
    Function to compute the spectral intensity integral (Eqs. 21 & 22), using
    a Taylor expansion method.
    """

    sinc_chi1dtau = sincu(0.5*chi1*dtau)
    cos_chi1dtau = sincu(0.5*chi1*dtau)
    I0  = sinc_chi1dtau*dtau
    I1  = (dtau*inv_chi1)*(sinc_chi1dtau - cos_chi1dtau)
    I2  = (0.25*dtau*dtau*dtau) * sinc_chi1dtau - 2*inv_chi1*I1
    
    ReIn = v0*I0
    ImIn = v1*I1+sign*(chi2*v0*I2)
    
    return ReIn, ImIn

def get_d2I(w, s_hat, gamma, x4, x41, x42, v0, v1, dtau, theta, small=1e-3):
    
    # Preallocate
    d2I = np.zeros(len(w), dtype = "complex128")
    #dchi2 = np.zeros((len(w), len(x4[0])))
    # Verbosity
    status = len(w) / 10
            
    for i in range(len(w)):
        if (i+1)%status==0:
            print(i+1, "of", len(w))
        om = w[i]
        kappa = om * inv_c * np.array([1, s_hat[0], s_hat[1], s_hat[2]])
        
        chi0  = -kappa[0]*x4[0] + kappa[1]*x4[1] + kappa[2]*x4[2]\
                + kappa[3]*x4[3]
                
        chi1 = -kappa[0]*x41[0] + kappa[1]*x41[1] + kappa[2]*x41[2] \
                + kappa[3]*x41[3]
                
        chi2 = -kappa[0]*x42[0] + kappa[1]*x42[1] + kappa[2]*x42[2] \
              + kappa[3]*x42[3]
        
        
        chi0 = np.longdouble(chi0)
        chi1 = np.longdouble(chi1)
        chi2 = np.longdouble(chi2)
        
        
        
        #dchi2[i, :] = chi2
        # Only divide once
        inv_chi1         = 1/chi1
        inv_chi2         = 1/chi2
        z_ind            = np.argwhere(chi2 < 0)
        sqrt_chi2        = np.complex128(np.sqrt(abs(chi2)))
        sqrt_chi2[z_ind] = 1j*sqrt_chi2[z_ind]
        inv_sqrt_chi2    = 1/sqrt_chi2
        
        # Fresnel approach
        phi_p = (0.25*dtau*dtau)*chi2 + (0.5*dtau*chi1)
        phi_m = (0.25*dtau*dtau)*chi2 + (0.5*dtau*chi1)
        
        theta_p = inv_sqrt_2pi*inv_sqrt_chi2*(chi1+chi2*dtau)
        theta_m = inv_sqrt_2pi*inv_sqrt_chi2*(chi1-chi2*dtau)
        
        
        psi_base_x = sqrt_2pi*inv_sqrt_chi2*(2*chi2*v0[0]-chi1*v1[0])
        psi_base_y = sqrt_2pi*inv_sqrt_chi2*(2*chi2*v0[1]-chi1*v1[1])
        psi_base_z = sqrt_2pi*inv_sqrt_chi2*(2*chi2*v0[2]-chi1*v1[2])
    
        psi_px     = psi_base_x*np.cos(0.25*chi1*chi1*inv_chi2)
        psi_py     = psi_base_y*np.cos(0.25*chi1*chi1*inv_chi2)
        psi_pz     = psi_base_z*np.cos(0.25*chi1*chi1*inv_chi2)
    
        psi_mx     = psi_base_x*np.sin(0.25*chi1*chi1*inv_chi2)
        psi_my     = psi_base_y*np.sin(0.25*chi1*chi1*inv_chi2)
        psi_mz     = psi_base_z*np.sin(0.25*chi1*chi1*inv_chi2)
        
        Sp, Cp = get_fresnel(np.complex128(theta_p))
        Sm, Cm = get_fresnel(np.complex128(theta_m))
        
        sin_p  = np.sin(phi_p)
        sin_m  = np.sin(phi_m)
        cos_p = np.cos(phi_p)
        cos_m = np.cos(phi_m)
        
        RIfx  = (0.25*inv_chi2) * (psi_px*(Cp-Cm) \
                                  + psi_mx*(Sp-Sm) \
                                  +2*v1[0]*(sin_p-sin_m))
        RIfy  = (0.25*inv_chi2) * (psi_py*(Cp-Cm) \
                                  + psi_my*(Sp-Sm) \
                                  +2*v1[1]*(sin_p-sin_m))
        RIfz  = (0.25*inv_chi2) * (psi_pz*(Cp-Cm) \
                                  + psi_mz*(Sp-Sm) \
                                  +2*v1[2]*(sin_p-sin_m))
        
        IIfx  = (0.25*inv_chi2)*(psi_px*(Sp-Sm)\
                                -psi_mx*(Cp-Cm)\
                                -2*v1[0]*cos_p-cos_m)
        
        IIfy  = (0.25*inv_chi2)*(psi_py*(Sp-Sm)\
                                -psi_my*(Cp-Cm)\
                                -2*v1[1]*cos_p-cos_m)
        
        IIfz  = (0.25*inv_chi2)*(psi_pz*(Sp-Sm)\
                                -psi_mz*(Cp-Cm)\
                                -2*v1[2]*cos_p-cos_m)
        
        # Taylor approximation
        
        sinc_chi1dtau = sincu(0.5*chi1*dtau)
        cos_chi1dtau = sincu(0.5*chi1*dtau)
        I0  = sinc_chi1dtau*dtau
        I1  = (dtau*inv_chi1)*(sinc_chi1dtau - cos_chi1dtau)
        I2  = (0.25*dtau*dtau*dtau) * sinc_chi1dtau - 2*inv_chi1*I1
        
        
        RItx = v0[0]*I0
        RIty = v0[1]*I0
        RItz = v0[2]*I0
        
        IItx = v1[0] * I1 + (chi2 * v0[0] * I2)
        IIty = v1[1] * I1 + (chi2 * v0[1] * I2)
        IItz = v1[2] * I1 + (chi2 * v0[2] * I2)
        
        # Apply taylor expansion filter
        comp_par = abs(chi2*dtau*dtau)
        t_ind    = np.argwhere(comp_par <= small)
        
        RIx        = RIfx
        RIx[t_ind] = RItx[t_ind]
        RIy        = RIfy
        RIy[t_ind] = RIty[t_ind]
        RIz        = RIfz
        RIz[t_ind] = RItz[t_ind]
        
        IIx        = IIfx
        IIx[t_ind] = IItx[t_ind]
        IIy        = IIfy
        IIy[t_ind] = IIty[t_ind]
        IIz        = IIfz
        IIz[t_ind] = IItz[t_ind]
        
        # Cartesian form 
        RSx = np.sum(RIx * np.cos(chi0) - IIx * np.sin(chi0))
        RSy = np.sum(RIy * np.cos(chi0) - IIy * np.sin(chi0))
        RSz = np.sum(RIz * np.cos(chi0) - IIz * np.sin(chi0))
    
        ISx = np.sum(IIx * np.cos(chi0) + RIx * np.sin(chi0))
        ISy = np.sum(IIy * np.cos(chi0) + RIy * np.sin(chi0))
        ISz = np.sum(IIz * np.cos(chi0) + RIz * np.sin(chi0))  
        
        d2I[i] = RSx**2 + ISx**2 \
                 + (RSy * np.cos(theta) - RSz * np.sin(theta))**2 \
                 + (ISy * np.cos(theta) - ISz * np.sin(theta))**2
    
  
    d2I = d2I * (mu_0*e*e*c*w*w*gamma*gamma)*(inv_pi*inv_pi*inv_pi)/(16)    
    
    
    return d2I