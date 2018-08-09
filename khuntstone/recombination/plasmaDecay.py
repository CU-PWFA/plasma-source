#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 14:51:46 2018
plasma decay based on coefficients from Simposon, Hinnov, and Gurachev
@author: keenan
"""
import coeffs as cf
from consts import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

# Flat top plasma density 
n0 = 5e16; #cm^-3;
# range of temperature
T = np.array([1000, 5000, 10000]);
# Recombination coefficients
def getAlphas(n0, T, kb = kb_eV):
    alpha_s  = cf.simpson(T)[0]; # m^6/s
    alpha_s0 = cf.simpson(T)[1]; # m^6/s
    alpha_si = cf.simpson(T)[2]; # m^6/s
    alpha_hh = cf.hh(kb * T, n0); # cm^3/s
    alpha_g  = cf.gurachev(T, n0); # cm^3/s
    return alpha_s, alpha_s0, alpha_si, alpha_hh, alpha_g

# get half-life for gurachev and HH
def getHalfLife(n0, alpha):
    return 1/(n0 * alpha);
# get density as a function of time for gurachev and HH
def getN(n0, alpha, t):
    n = np.zeros(np.shape(t), dtype = 'float')
    for i in range(len(t)):
        n[i] = ((1./n0) + alpha * t[i])**(-1.);
    return n
# get half_life for Simpson
def getHalfLife_s(n0, alpha):
    return 1/((n0*1e6)**2 * alpha);
def getN_s(n0, alpha, t):
    n = np.zeros(np.shape(t), dtype = 'float')
    for i in range(len(t)):
        n[i] = ((1./(n0 * 1e6)) + alpha * n0 * 1e6 * t[i])**(-1.);
    return n
# Plot decay for simpson cases
def plotDecayS(n0, T, alpha):
    for i in range(len(T)):
        t_end = getHalfLife_s(n0, alpha[i])
        t = np.linspace(0, t_end, 10000)
        n = getN_s(n0, alpha[i], t)
        plt.loglog(t, n*1e-6, label = 'T = ' + str(T[i]))
    plt.ylabel(r'$n_e$ [$cm^{-3}$]')
    plt.xlabel(r't [s]')
    if len(T) > 1:
        plt.legend()
    else:
        plt.title('Simpson Plasma Decay - T = ' + str(T[0]) + 'K')
    plt.show()

# plot decay for HH and Gurachev cases
def plotDecay(n0, T, alpha):
    for i in range(len(T)):
        t_end = getHalfLife(n0, alpha[i])
        t = np.linspace(0, t_end, 1000)
        n = getN(n0, alpha[i], t)
        plt.loglog(t, n, label = 'T = ' + str(T[i]))
    plt.ylabel(r'$n_e$ [$cm^{-3}$]')
    plt.xlabel(r't [s]')
    if len(T) > 1:
        plt.legend()
    else:
        plt.title('Simpson Plasma Decay - T = ' + str(T[0]) + 'K')
    plt.show()
# Plot half-lives for a variety of temperatures
def plotHalfT(n0, T):
    alphas    = getAlphas(n0, T);
    t_half_s  = [getHalfLife_s(n0, alphas[i]) for i in range(3)];
    t_half_hh = getHalfLife(n0, alphas[3])
    t_half_g  = getHalfLife(n0, alphas[4])
    labels = [r'$\alpha_s$', r'$\alpha_{s,0}$', r'$\alpha_{s,\infty}$']
    for i in range(3):
        plt.loglog(T, t_half_s[i], label = labels[i])
    plt.loglog(T, t_half_hh, label = r'$\alpha_{HH}$')
    plt.loglog(T, t_half_g, label = r'$\alpha_{G}$')
    plt.legend()
    plt.xlabel('T [K]')
    plt.ylabel(r'$\tau$ [s]')
    plt.ylim([1e-15, 1e0])
    plt.show()
#T = np.linspace(1e2, 15e3, num = 1000)
#plotHalfT(n0, T)
# rough estimate of temperature from ponderomotive energy
def getT(I, m = me_SI, eps = eps0, e = e_SI, kb = kb_SI):
    c = 3e8; # m/s
    l = 800e-9; # m
    f = c/l; # /s carrier frequency
    U = e**2 * I / (2 * c * eps0 * m * f**2)
    return (2./3.) * U / kb
I = np.linspace(1e18, 1e19, 1000)
T = getT(I);
T = T / 11600; # eV
plt.plot(I, T)
plt.xlabel('I [$W/m^2$]')
plt.ylabel('T [eV]')
plt.show()
    
    
    
    
    
    
    
    
    