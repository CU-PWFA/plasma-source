#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 10:36:53 2022

@author: chris
"""



# Add python modules to path and import
import numpy as np
import sys
# Be sure to include your own path to the python modules
sys.path.insert(0, "/home/chris/Desktop/plasma-source/plasma-source/khuntstone/eos_bpm/python/")
import currents as cp
import eo_signal as eos
from plotting import plot_signal

from scipy.constants import c
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import currents as cp
from crystal import crystal
from laser import laser
import phase_retard as pr
from plotting import makefig
import thz

# Create setup dictionary
# NOTE: remember to change fpath to your local machine's directory
cpath = "/home/chris/Desktop/DataLoads/EOSBeamProfiles/"
setup = {
    "ctype"   : "gap",
    "d"       : 100e-6,
    "y0"      : 800e-9,
    "tp"      : 30e-15,
    "r0"      : 2.5e-3,
    "fpath"   : cpath,
    "nslice"  : 100,
    "tau"     : np.linspace(-1000, 1000, 1000)*1e-15,
    "method"  : "spatial",
    "process" : "near",
    "angle"   : 15,
    "theta"   : 25 # You can see this causes phase wrapping at certain angles
    
}
ind = 0 # Current profile index to use
sig, t_sig = eos.get_signal(ind, setup)

# Get Electric field for plotting
I, ti, p2p = cp.get_current(ind, setup["fpath"])
E, te      = cp.get_E(I, ti, setup["r0"])
plot_signal(E, te, sig, t_sig)