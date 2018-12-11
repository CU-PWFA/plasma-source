#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 11:09:36 2018
Comparing recombination rates from Hinnov and Hirschberg, Hoffert and Lien, 
Simpson, and Gurachev
@author: keenan
"""
import sys
sys.path.insert(0, "../Constants")
import SI
from consts import *
import numpy as np
from matplotlib import pyplot as plt

kb_SI = SI.boltzmann
kb_eV = kb_SI * 6.242e18;
h_SI  = SI.planck
me    = SI.elecMass
# Simpson method calculate coeff in m^6/s as a polynomial of temperature

# Polynomial coefficients
def coeffs():
    C0 = np.array([2.647e-15, 3.499e-16, 2.817e-15, 1.878e-10]);
    C1 = np.array([1.622e-18, 1.407e-19, 2.476e-18, -9.932e-16]);
    C2 = np.array([-9.615e-23, 2.406e-23, -1.956e-22, -1.489e-18]);
    C3 = np.array([5.824e-27, -3.612e-29, 9.750e-27, 1.297e-22]);
    C4 = np.array([-1.373e-31, -2.347e-32, -2.026e-31, -3.429e-27]);
    C = np.array([C0, C1, C2, C3, C4]);
    return C
def getPolys(T):
    # T in K
    C = coeffs();
    S0i = C[0][0] + C[1][0]*T + C[2][0]*T**2 + C[3][0]*T**3 + C[4][0]*T**4;
    S0s = C[0][1] + C[1][1]*T + C[2][1]*T**2 + C[3][1]*T**3 + C[4][1]*T**4;
    P0s = C[0][2] + C[1][2]*T + C[2][2]*T**2 + C[3][2]*T**3 + C[4][2]*T**4;
    Psi = C[0][3] + C[1][3]*T + C[2][3]*T**2 + C[3][3]*T**3 + C[4][3]*T**4;
    return S0i, S0s, P0s, Psi;

# degeneracies
g0 = 1.;
gi = 6.;
# energy difference to continuum
es = 4.21; # eV
# Coefficent for alpha
def getC(T, kb = kb_SI, h = h_SI, me = me_SI):
    return .5 * (2 * np.pi * me * kb * T / h**2)**(-1.5);
def getBranch(T, g0 = g0, es = es, kb = kb_eV):
    # g0 - degeneracy, es - energy difference to continuum for Ar
    # Get polynomials needed (S0s, P0s, Psi)
    polys = getPolys(T)
    S0s = polys[1]; P0s  = polys[2]; Psi = polys[3];
    b = g0 * (S0s * np.exp(es/(kb * T)) - P0s)/Psi
    return b, polys
# Recombination coeff from Simpson 
def simpson(T, g0 = g0, gi = gi, es = es, kb = kb_eV ):
    [b, polys] = getBranch(T);
    C = getC(T);
    S0i = polys[0]; S0s = polys[1];
    # Three alphas for b = , b= 0 ,b = inf
    alpha   = (C * g0/gi) * (S0i + S0s * np.exp(es/(kb*T))/(1 + b));
    alpha_0 = (C * g0/gi) * (S0i + S0s * np.exp(es/(kb*T))/(1 + 0));
    alpha_i = (C * g0/gi) * (S0i + S0s * np.exp(es/(kb*T))/(1 + np.inf));
    return alpha, alpha_0, alpha_i
# Plot for temperature range 
def plotSimpson(T):
    [alpha, alpha_0, alpha_i] = simpson(T);
    b = getBranch(T)[0];
    plt.semilogy(T/1e3, alpha * 1e40,'-b', label = r'$\alpha$')
    plt.semilogy(T/1e3, alpha_0 * 1e40, '-r', label = r'$\alpha_0$ (b = 0)')
    plt.semilogy(T/1e3, alpha_i * 1e40, '-g', label = r'$\alpha_i$' +\
                                                     ' (b = $\infty$)')
    plt.semilogy(T/1e3, b,'-y', label = 'b')
    plt.legend()
    plt.xlabel(r'$T_e$ [$10^3$ K]')
    plt.ylabel(r'$\alpha$ [$10^{-40} m^6 s^{-1}$]')
    # Lines for simpson range
    xpos = [2, 16]
    for xc in xpos:
        plt.axvline(x = xc, color = 'k', linestyle = '--')
    plt.ylim([1e-3, 1e4])
    plt.show()
#T = np.arange(1, 25e3, 10);
#plotSimpson(T)

# Hinnov and Hirschberg approach, coeff is a function of density and temp
def hh(kT, ne):
    ar = 2.7e-13 * kT**(-3./4.);
    ac = 5.6e-27 * kT**(-9./2.) * ne;
    return ac + ar;

def plotHH(kT, ne):   
    for i in kT:
        alpha = hh(i, ne);
        plt.plot(np.log10(ne), np.log10(alpha), label = 'kT =' +str(i))
    plt.xlabel(r'log$_{10}$' + ' ' + '$n_e$')
    plt.ylabel(r'log$_{10}$' + ' ' + '$\\alpha$')
    plt.ylim([-13, -7])
    plt.legend()
    plt.show()
kT = np.array([0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0])
ne = np.linspace(1e7, 1e15, num = 1e5)
#plotHH(kT, ne)

# Gurachev approach, coeff is a function of density and temp
def gurachev(T, ne, kb = kb_cgs, Z = 1, e = e_cgs, m = me_cgs):
    num = 4 * np.sqrt(2) * np.pi**(3./2.) * e**10 * Z**3 * ne * \
          np.sqrt(Z**2 + 1);
    den = 9 * np.sqrt(m) * (kb * T)**(9./2.)
    return num/den;
def plotGurachev(T, ne, kb = kb_eV):
    for i in T:
        alpha = gurachev(i, ne)
        plt.plot(np.log10(ne), np.log10(alpha), label = 'kT = ' + str(i*kb))
    plt.xlabel(r'log$_{10}$' + ' ' + '$n_e$')
    plt.ylabel(r'log$_{10}$' + ' ' + '$\\alpha$')
    plt.legend()
    plt.show()
T = kT / kb_eV;
#plotGurachev(T, ne)
def plotHH_Gurachev(ne, kT = 0.2):
    alpha_g = gurachev(kT/kb_eV, ne)
    alpha_hh = hh(kT, ne);
    plt.plot(np.log10(ne), np.log10(alpha_g), label = 'Gurachev')
    plt.plot(np.log10(ne), np.log10(alpha_hh), label = 'Hinnov & Hirschberg')
    plt.xlabel(r'log$_{10}$' + ' ' + '$n_e$')
    plt.ylabel(r'log$_{10}$' + ' ' + '$\\alpha$')
    plt.legend()
    plt.show()
#plotHH_Gurachev(ne)

def desai(T, ne):
    return 1.28e5 * T**(-1.8) * 10**(-(3410/T)) * ne**(-0.64);

