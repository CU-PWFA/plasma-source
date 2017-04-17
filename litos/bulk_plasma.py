#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 16:35:13 2017

@author: litos
"""

## common libraries
#import sys
from numpy import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt

## custom functions
from fitramp import fitramp
from propbeam import propbeam
from plasma_ramp import plasma_ramp
from draw_ellipse import draw_ellipse
from calc_Bmag import calc_Bmag
from calc_M import calc_M

## define constants
c  = 3e8 # m/s, speed of light
me = 0.511e-3 # GeV, rest mass energy of electron

## set plasma parameters
# use values for main accelerator plasma
np0 = 5e16 # cm^-3, plasma number density
wp0 = (5.64e4)*sqrt(np0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
dEdm = 4.0 # GeV/m, accelerating gradient
dgdz = dEdm/me # change in gamma per unit length

## set initial beam parameters
# use values for nominal vacuum waist
gb0   = 20000 # Lorentz factor
eps0   = 5e-6 # m-rad
kb0 = kp0/sqrt(2*gb0) # m^-1, betatron wave number
beta0  = 1/kb0 # m
alpha0 = 0
gamma0 = (1+(alpha0**2))/beta0 # 1/m
T0     = [beta0,alpha0,gamma0] # Twiss vector
Ttarg  = [beta0,alpha0,gamma0]

# Lp : full length of bulk plasma
Lp0 = 1.00 # m
nz  = 100001
z   = linspace(0,Lp0,nz)
dz  = abs(z[1]-z[0])

beta = zeros(nz)
beta[0] = beta0
M = zeros(nz)
M[0] = 1
beam = zeros((nz,5))
beam[0][:] = [gb0,eps0,beta0,alpha0,gamma0]
for i in range(1,len(z)):
    ibeam = propbeam(beam[i-1][:],[0,dz],[np0,np0])
    beam[i][:] = ibeam[len(ibeam)-1][:]
    beam[i][0] = beam[i][0]+dgdz*dz
    beta[i]    = beam[i][2]
    iTbeam     = [beam[i][2],beam[i][3],beam[i][4]]
    M[i]       = calc_M(iTbeam,Ttarg)

plt.figure()
plt.plot(z,beta)
plt.title('Beam evolution in bulk plasma')
plt.xlabel(r'$z$ (m)')
plt.ylabel(r'$\beta$ (m)')
plt.show()

plt.figure()
plt.plot(z,M)
plt.title('Mismatch evolution in bulk plasma')
plt.xlabel(r'$z$ (m)')
plt.ylabel(r'$M$ (m)')
plt.show()

z_ind = argwhere(M > 1.1 )
z_lim = z[z_ind[0]]
print(z_lim)
beta_lim = beta[z_ind[0]]
print(beta_lim/beta0)
print(beam[z_ind[0]][0])