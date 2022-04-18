#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 11:40:27 2022

@author: chris
"""



# Notebook for developing new code, changes often
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, epsilon_0
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
try:
    plt.style.use("huntstone")
except:
    plt.style.use("seaborn-colorblind")
inv_c = 1/c
from scipy.interpolate import interp1d
import sys
sys.path.insert(0, "/home/chris/Desktop/plasma-source/plasma-source/cedoss/eos_doss/python/")
# Custom modules
from crystal import crystal
import currents as cp
import eo_signal as eos
from laser import laser
import phase_retard as pr
import plotting as pl
import thz

#import ipywidgets as widgets
from IPython.display import clear_output
#outs = widgets.Output()
import time
cpath = "/home/chris/Desktop/DataLoads/EOSBeamProfiles/"
#peak_sigs      = np.load("peakSignals.npy")
#peak_drive_sig = peak_sigs[0]
#peak_wit_sig   = peak_sigs[1]

# Testing 1mm ZnTe for first optomizeing signal
ctype  = "znte"
d      = 100e-6
y0     = 800e-9
tp     = 45e-15
nslice = 100
dtau   = 30
tau    = np.arange(-500, 1000+dtau, dtau)*1e-15
setup  = {"ctype"   : ctype,
         "d"       : d,
         "y0"      : y0,
         "tp"      : tp,
         "nslice"  : nslice,
         "method"  : "spatial",
         "process" : "cross", 
         "r0"      : 2.5e-3,
         "fpath"   : cpath,
         "tau"     : tau,
         "angle"   : 15}

sig, tsig = eos.get_signal(10, setup)
I, ti, p2p = cp.get_current(10, cpath)
"""
fig = plt.figure(figsize = (4,4), dpi = 300)
ax1  = fig.gca()
ax1.set_xlabel("t [ps]")
ax1.set_ylabel("Signal [AU]", color = "blue")
ax1.tick_params(axis = "y", color = "b", labelcolor = "blue")
ax1.plot(tsig*1e12, sig/max(sig), '-b')
ax2 = ax1.twinx()
ax2.set_ylabel("I [kA]", color = "r")
ax2.spines["left"].set_color("b")
ax2.spines["right"].set_color("r")
ax2.tick_params(axis = "y", color = "r", labelcolor = "r")
ax2.plot(ti*1e12, I, '-r')
ax2.set_xlim([-0.5, 0.5])
ax2.set_title(r'100 $\mu$m  ZnTe')
plt.show()
"""
# Test ZnTe for single bunch

eps0 = epsilon_0
Q    = 0.4e-12
sigz = 50e-6
r0   = 2.5e-3
sigt = sigz/c
E0   = Q / ((2*np.pi)**(1.5)*sigt*epsilon_0*r0*c)
dt   = sigt*0.1
N    = 8000
te   = np.linspace(-N*dt*0.5, N*dt*0.5, N)
E    = E0 * np.exp(-te*te/(2*sigt*sigt))
setup["d"] = 0.5e-3
setup["ctype"] = "znte"
print(setup["ctype"], setup["d"])
setup["tau"] = te
sig_z, tsig_z, gamma, tgamma = eos.E_signal(E, te, setup)
setup["d"] = 100e-6
setup["ctype"] = "gap"
print(setup["ctype"], setup["d"])
sig_g, tsig_g, gamma, tgamma = eos.E_signal(E, te, setup)

"""plt.plot(E)
plt.show()
plt.plot(tsig_g*1e12, sig_g)
plt.show()
plt.plot(tsig_z*1e12, sig_z)
plt.show()"""

fig = plt.figure(figsize = (2,2), dpi = 300)
ax1  = fig.gca()
ax1.set_xlabel("t [ps]")
ax1.set_ylabel("Signal [AU]", color = "blue")
ax1.tick_params(axis = "y", color = "b", labelcolor = "blue")
ax1.plot(tsig_z*1e12, sig_z/max(sig_z), '-b', label = "ZnTe")
#ax1.plot(tsig_g*1e12, sig_g/max(sig_z), '-.b', label = "GaP")
ax1.legend()
ax2 = ax1.twinx()
ax2.set_ylabel("E [V/m]", color = "r")
ax2.spines["left"].set_color("b")
ax2.spines["right"].set_color("r")
ax2.tick_params(axis = "y", color = "r", labelcolor = "r")
ax2.plot(te*1e12, E, '-r')
ax2.set_xlim([-0.5, 0.5])
ax2.set_title(r'1 mm Crystal')
plt.show()