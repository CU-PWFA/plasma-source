#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 19:32:52 2021

@author: litos
"""

import numpy as np
from scipy.constants import c, e, m_e, mu_0
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.special import erf
from scipy.special import kv
import matplotlib.pyplot as plt

pi = np.pi
inv_pi = 1/pi
inv_c  = 1/c

# unnormalized sinc function
def sincu(x):
    return np.sin(x)/x

# calculate product of four-vectors
def four_prod(a,b):
    return a[0]*b[0] - np.dot(a[1:],b[1:])

def get_fresnel(z):
    ip    = 1 + 1j
    im    = 1 - 1j
    arg_p = np.complex128(ip * np.sqrt(np.pi) * z * 0.5)
    arg_m = np.complex128(im * np.sqrt(np.pi) * z * 0.5)
    erf_p = erf(arg_p)
    erf_m = erf(arg_m)
    S     = (0.25*ip) * (erf_p - 1j*erf_m)
    C     = (0.25*im) * (erf_p + 1j*erf_m)
    return S, C

# input parameters
w0   = 2.5e15
Nw   = 10000
w    = w0*np.linspace(1e0,1e7,Nw)

shat   = np.array([1,0,0],dtype="longdouble")
kaphat = np.concatenate(([inv_c],shat*inv_c))

gamma     = np.longdouble(1000)
inv_gamma = 1/gamma
beta      = np.sqrt(1-inv_gamma**2)

B       = np.longdouble(m_e*w0/e)
rho     = np.longdouble(gamma*beta*m_e*c/(B*e))
inv_rho = 1/rho

dtau = pi*(1e4)/w0
Ntau = Nw
tau  = np.linspace(0,Ntau*dtau,Ntau)
t    = gamma * tau

phi_obs = 90*pi/180

# position vectors
x  = rho * np.sin(beta*c*t * 0.5*inv_pi*inv_rho)
y  = rho * (1 - np.cos(beta*c*t * 0.5*inv_pi*inv_rho))
z  = np.zeros(len(x),dtype="longdouble")
x4 = np.array([c*t,x,y,z],dtype="longdouble")

# get quadratic interpolation coefficients
def get_interp_coeffs(tau,x,x_spl):
    xp_spl  = x_spl.derivative(n=1)
    xpp_spl = x_spl.derivative(n=2)
    v0 = xp_spl(tau)
    v1 = xpp_spl(tau)
    x2 = (0.5*v1)
    x1 = v1-2*c*tau
    x0 = x_spl(tau) - x1*tau - x2*tau*tau
    # x0 = x
    return x0, x1, x2, v0, v1

# create 2nd order splines
x_spl = ius(tau,x,k=2)
y_spl = ius(tau,y,k=2)
z_spl = ius(tau,z,k=2)

# quadratic coefficient vectors
x0,x1,x2,vx0,vx1 = get_interp_coeffs(tau,x,x_spl)
y0,y1,y2,vy0,vy1 = get_interp_coeffs(tau,y,y_spl)
z0,z1,z2,vz0,vz1 = get_interp_coeffs(tau,z,z_spl)

# quadratic coefficient position four-vectors
zero_array = np.zeros(len(tau),dtype="longdouble")
x40 = np.array([c*t,x0,y0,z0],dtype="longdouble")
x41 = np.array([zero_array + gamma * c,x1,y1,z1],dtype="longdouble")
x42 = np.array([zero_array,x2,y2,z2],dtype="longdouble")

# quadratic coefficient velocity three-vectors
v0 = inv_c*np.array([vx0,vy0,vz0])
v1 = inv_c*np.array([vx1,vy1,vz1])

# chi four-vector scalar products
chi0 = np.complex128(w*four_prod(kaphat,x40))
chi1 = np.complex128(w*four_prod(kaphat,x41))
chi2 = np.complex128(w*four_prod(kaphat,x42))

# theta scalars
sqrt_2pi      = np.sqrt(2*pi)
inv_sqrt_2pi  = 1/sqrt_2pi
sqrt_chi2     = np.sqrt(chi2)
inv_sqrt_chi2 = 1/sqrt_chi2
inv_chi2      = inv_sqrt_chi2*inv_sqrt_chi2

theta_p = (chi1+chi2*dtau)*inv_sqrt_2pi*inv_sqrt_chi2
theta_m = (chi1-chi2*dtau)*inv_sqrt_2pi*inv_sqrt_chi2

# psi three-vectors
psi_pre_fact = sqrt_2pi*inv_sqrt_chi2*(2*chi2*v0 - chi1*v1)
psi_arg      = 0.25*chi1*chi1*inv_chi2

psi_p = psi_pre_fact*np.cos(psi_arg)
psi_m = psi_pre_fact*np.sin(psi_arg)

# phi scalars
phi_p = 0.25*dtau*dtau*chi2 + 0.5*dtau*chi1
phi_m = 0.25*dtau*dtau*chi2 - 0.5*dtau*chi1

# Fresnel integrals
S_theta_p, C_theta_p = get_fresnel(theta_p)
S_theta_m, C_theta_m = get_fresnel(theta_m)

# Real and Imaginary integral solutions using Fresnel integrals
sin_phi_p = np.sin(phi_p)
sin_phi_m = np.sin(phi_m)
cos_phi_p = np.cos(phi_p)
cos_phi_m = np.cos(phi_m)

RIf = 0.25*inv_chi2*( psi_p*(C_theta_p - C_theta_m) \
                    + psi_m*(S_theta_p - S_theta_m) \
                    + 2*v1*(sin_phi_p - sin_phi_m))

IIf = 0.25*inv_chi2*( psi_p*(S_theta_p - S_theta_m) \
                    - psi_m*(C_theta_p - C_theta_m) \
                    - 2*v1*(cos_phi_p - cos_phi_m))

# Real and Imaginary integral solutions using Taylor expansion
sinc_chidtau = sincu(0.5*chi1*dtau)
cos_chidtau  = np.cos(0.5*chi1*dtau)

I0 = sinc_chidtau*dtau
I1 = (dtau/chi1)*(sinc_chidtau - cos_chidtau)
I2 = 0.25*dtau*dtau*dtau*sinc_chidtau - 2*I1/chi1

RIt = v0*I0
IIt = v1*I1 + chi2*v0*I2

# choose Fourier or Taylor integral solutions
comp_par = np.absolute(chi2*tau*tau)
small = 1e-3
RI = RIf
II = IIf
for i in range(3):
    RI[i,comp_par<=small] = RIt[i,comp_par<=small]
    II[i,comp_par<=small] = IIt[i,comp_par<=small]

# cartesian form of integral solutions
cos_chi0 = np.cos(chi0)
sin_chi0 = np.sin(chi0)

RS = np.sum(RI*cos_chi0 - II*sin_chi0,1)
IS = np.sum(II*cos_chi0 + RI*sin_chi0,1)

# spectral intensity
cos_theta_obs = np.cos(phi_obs)
sin_theta_obs = np.sin(phi_obs)
d2I_prefact = 0.065*mu_0*e*e*c*inv_pi*inv_pi*inv_pi
d2I_RS = np.complex128(RS[1]*cos_theta_obs - RS[2]*sin_theta_obs)
d2I_IS = np.complex128(IS[1]*cos_theta_obs - IS[2]*sin_theta_obs)

d2I_thomas = np.complex128(d2I_prefact * w*w * \
                           (RS[0]*RS[0] + IS[0]*IS[0] \
                            + d2I_RS*d2I_RS + d2I_IS*d2I_IS) )
    
def analytic_synch(w, gamma, rho, theta_obs):
    t1 = (mu_0 * e**2 * c * w**2 / (12 * np.pi**3)) * (rho / c)**2
    t2 = ((1 / gamma**2) + theta_obs**2)**2
    xi = w * rho / (3 * c * gamma**3 * (1 + theta_obs**2 * gamma**2)**(1.5))
    t3 = kv(2/3, xi)**2 + ((theta_obs**2 * kv(1/3, xi)**2) / ((1 / gamma**2) \
                           + theta_obs**2))
    return t1 * t2 * t3
    
theta_obs = 0*pi/180
d2I_ana = analytic_synch(w, gamma, rho, theta_obs)

# Plot & Compare
fig1 = plt.figure()
ax1  = fig1.add_subplot(121)
ax2  = fig1.add_subplot(122)
ax1.set_xlabel(r'$\omega / \omega_0$')
ax1.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$')
ax1.set_title("Thomas")
ax2.set_xlabel(r'$\omega / \omega_0$')
ax2.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$')
ax2.set_title("Analytical")
ax1.semilogx(w / w0, np.real(d2I_thomas))
ax2.semilogx(w / w0, d2I_ana)
plt.show()

fig2 = plt.figure()
ax3  = fig2.gca()
ax3.set_xlabel(r'$\omega / \omega_0$')
ax3.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$')
ax3.semilogx(w / w0, d2I_thomas / max(d2I_thomas), label = "Thomas")
ax3.semilogx(w / w0, d2I_ana / max(d2I_ana), label = "Analytic")
ax3.legend()
plt.show()