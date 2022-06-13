#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 09:23:31 2022

@author: chris
"""

# Add python modules to path and import
import numpy as np
import sys
# Be sure to include your own path to the python modules
sys.path.insert(0, "/home/chris/Desktop/plasma-source/plasma-source/cedoss/eos_doss/python/")
import currents as cp
import eo_signal as eos
from plotting import plot_signal#, plot_gamma
from scipy.constants import c, epsilon_0, elementary_charge
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import currents as cp
from crystal import crystal
from laser import laser
import phase_retard as pr
from plotting import makefig
import thz

def las_int(x, y, sigr = 0.25e-3):
    I_prof = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            r = np.sqrt(x[i]*x[i]+y[j]*y[j])
            if r > 1.5e-3:
                I_prof[i,j] = 0
            else:
                I_prof[i,j] = np.exp(-r*r/(2*sigr*sigr))
    return I_prof
def get_E(Q, sigz, r):
    sigt = sigz/c
    dt   = 0.1*sigt
    N    = 8000
    t    = np.linspace(-N*dt*0.5, 0.5*dt*N, N)
    E0   = Q/((2*np.pi)**(1.5)*epsilon_0*r*c*sigt)
    E    = E0*np.exp(-t*t/(2*sigt*sigt))
    return E, t

# Figure: EOS-BPM vs BPM comparison
# Beam paramters
Q    = 1.200e-9 
sigz = 20e-6
r0   = 1.3e-3

# Create setup dictionary
# NOTE: remember to change fpath to your local machine's directory
# cpath = "/media/keenan/Data_Storage/Data/currents/"
# cpath = "/mnt/md0/Data/currents/"
cpath = "/home/chris/Desktop/DataLoads/EOSBeamProfiles/"
setup = {
    "ctype"   : "gap",
    "d"       : 100e-6,
    "y0"      : 800e-9,
    "tp"      : 40e-15,
    "r0"      : r0,
    "fpath"   : cpath,
    "nslice"  : 100,
    "tau"     : np.linspace(-1000, 1000, 1000)*1e-15,
    "method"  : "spatial",
    "process" : "cross",
    "angle"   : 15,
    "theta"   : 25 # You can see this causes phase wrapping at certain angles
    
}

# Get Electric field for plotting
E, te  = get_E(Q, sigz, r0)
setup["tau"] = te
sig, tsig, gamma, tgamma = eos.E_signal(E, te, setup)
plot_signal(E, te, sig, tsig)
#plot_gamma(E, te, gamma, tgamma)