#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:04:48 2020

Propagating and ebeam through a plasma and dumping the 6D-phase space for 
calculating the betatron radiation.

@author: keenan
"""

# standard python code
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patches
import scipy.constants as const
from scipy import optimize
import psutil
import timeit
tic = timeit.default_timer() # time the process

# add WARGSim directory location to $PATH
sys.path.insert(0, "/home/keenan/WARGSim")
# import WARGSim packages
from beams import electronbeam
from beams import betatronbeam
import calc.electron as ecalc
from elements import pwfa as pwfa
from elements import drift as drift
from elements import quadrupole as quadrupole
from elements import dipole as dipole
from elements import beamline as beamline
#from optics import profile
from interactions import interactions

# Define output location. Don't forget final slash!
#path = '/home/mike/work/github_repos/CU-PWFA/plasma-source/litos/butterfly/'
path = '/home/keenan/Dumps/'

# Get the number of available cores
Ncores = psutil.cpu_count()
numThread  = Ncores # int, number of threads to run on

dumpPeriod = 1e9 # int, step period for electron beam dump

Nslice = 1000


#%%
# define electron beam

N       = int(1e6) # number of particles in the beam
# eV, electron rest mass energy:
me = const.physical_constants['electron mass energy equivalent in MeV'][0]*(1e6)
beamE   = 10e9 # beam energy in eV
gb0     = beamE/me # reltivistic Lorentz factor (centroid value)
eps_n0   = 3.0e-6 # m-rad, normalized emittance
beta_x0  = 0.25   # m, x beta function at z=0
beta_y0  = 0.25   # m, y beta function at z=0
alpha_x0 = 0.00   # x alpha function at z=0
alpha_y0 = 0.00   # y alpha function at z=0
#gamma_x0 = (1+alpha_x0**2)/beta_x0 # 1/m, x gamma function at z=0
#gamma_y0 = (1+alpha_y0**2)/beta_y0 # 1/m, y gamma function at z=0
rms_z0  = 0.00   # m, rms bunch length
rms_gb0 = 0.01 # relative energy spread


# Useful functions for checking emittance growth

def match_beta(np0,gb0):
    np0 = np0*(1e17)
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = 1.0/kbeta
    return beta

def calc_Bmag(beta,beta_m):
    return (1/2)*(beta/beta_m + beta_m/beta)

def calc_beta_fact(Bmag):
    return (Bmag+np.sqrt((Bmag**2) - 1))




# force matched Twiss params:
Bmag_target = 5.00



n0 = 0.34 # make sure this agrees with value used to define plasma!!
beta_m = match_beta(n0,gb0)
beta_x00 = calc_beta_fact(Bmag_target)*beta_m
beta_y00 = calc_beta_fact(Bmag_target)*beta_m
Bmag = calc_Bmag(beta_x0,beta_m)
print(f'beta factor = {beta_x0/beta_m :.2f}')
print(f'Bmag = {Bmag :.2f}')


electronParams = {
        'name'      : 'ElectronBeam',
        'path'      : path,
        'load'      : False,
        'N'         : N,
        'shape'     : 'Gauss',
        'gb0'       : gb0,
        'rms_gb'    : rms_gb0,
        'rms_z'     : rms_z0,
        'eps_nx'    : eps_n0,
        'eps_ny'    : eps_n0,
        'beta_x'    : beta_x00,
        'beta_y'    : beta_y00,
        'alpha_x'   : alpha_x0,
        'alpha_y'   : alpha_x0
    }
my_ebeam = electronbeam.ElectronBeam(electronParams)

eps_nx0, eps_ny0 = my_ebeam.get_emit_n()
print('init eps_nx = ',eps_nx0)
print('init eps_ny = ',eps_ny0)


#%%
# introduce offset at vacuum waist

sig_x = np.sqrt(eps_n0*beta_x0/gb0)
sig_xp = np.sqrt(eps_n0/(beta_x0*gb0))
sig_y = np.sqrt(eps_n0*beta_y0/gb0)
sig_yp = np.sqrt(eps_n0/(beta_y0*gb0))

#my_ebeam.ptcls[:,0] += 1.00*sig_x
#my_ebeam.ptcls[:,1] += 1.00*sig_xp
#my_ebeam.ptcls[:,2] += 1.00*sig_y
#my_ebeam.ptcls[:,3] += -1.00*sig_yp



#%%
# define the plasma

def match_hw(np0,gb0,beta0):
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = beta0*kbeta
    hw = (1.96e-3)*beta**2 + 1.49*beta - 3.08
    return hw/kbeta

def match_zw(np0,gb0,beta0):
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = beta0*kbeta
    zw = (-1.47e-2)*beta**2 - 4.34*beta + 17.9
    return zw/kbeta

#n0    = 0.34  # 1/(10^17 cm^-3), flat-top plasma density
shape = 'Gauss' # ramp shape
L_ft  = 0.40 # m, length of flat-top
hw_up = 0 # m, half-width of up-ramp
L_up  = 8*hw_up # m, full length of up-ramp
hw_dn = 1.00*match_hw(n0*(1e17),2*gb0,beta_x0) # m, half-width of down-ramp
L_dn  = min(1.600,8*hw_dn) #L_up  # m, full length of down-ramp
L_p   = L_up + L_ft + L_dn # m, full length of plasma region

pwfaParams ={
    'name'  : 'PWFA0',
    'L'     : L_p,
    'n0'    : n0,
    'autoNz': True,
    'gb0'   : gb0,
    'shape' : shape,
    'L_ft'  : L_ft,
    'hw_up' : hw_up,
    'hw_dn' : hw_dn
}
pwfa0 = pwfa.PWFA(pwfaParams)

dumpPeriod = int(len(pwfa0.dz) / 100)
#%%
# propagate the beam through the plasma

# propagate the beam from the waist back to the start of the simulation
#waist      = match_zw(n0*(1e17),gb0,beta_x0) # -0.3884 # m, waist distance upstream of plasma flat-top start
#Z_back     = -(L_up+waist)
#interactions.electron_vacuum(my_ebeam, np.linspace(0,Z_back,2), dumpPeriod, numThread)

# propagate the beam through the PWFA
#interactions.ebeam_bbeam_pwfa(my_ebeam, my_bbeam, pwfa0, dumpPeriod, numThread)
#interactions.ebeam_pwfa(my_ebeam, pwfa0, dumpPeriod, numThread)

ecalc.ebeam_prop_pwfa_dumps(my_ebeam.ptcls,
                          pwfa0.dz, pwfa0.ne, pwfa0.dgb, 
                          dumpPeriod, my_ebeam.save_ptcls,
                          numThread)