# -*- coding: utf-8 -*-
# For testing and running code

# Standard python imports
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, epsilon_0
from scipy.interpolate import interp1d

# Custom imports 
import eo_signal as eos
from plotting import makefig
import tilt as tl
# Config
with open("/home/keenan/Dropbox/notebooks/eos/conf.json") as json_conf: 
    CONF = json.load(json_conf)
fpath = CONF["cp_lnx"]
def I_sfel(I1, I2, I3, zi, zsep, s12, s3, ojit = 0.05):
    # Getting a swiss-xfel current profile
    tsep = zsep / c
    dt   = 0.5 * tsep
    ti   = zi / c
    cur1 = I1 * np.exp(-(ti + dt)**2 / (2 * s12**2))
    cur2 = I2 * np.exp(-(ti - dt)**2 / (2 * s12**2))
    cur3 = I3 * np.exp(-ti**2 / (2 * s3**2))
    cur  = cur1 + cur2 + cur3
    # Add jitter
    sign = np.random.randint(2, size = len(cur))
    sign[sign == 0] = -1
    jitter = sign * np.random.random(size = len(cur)) * ojit
    return cur * (1 + jitter), ti


def main():
    I1        = 2.5e3
    I2        = 1.9e3
    I3        = 3e3
    zi        = np.linspace(-2.5e-5, 2.5e-5, 1000)
    zsep      = 14.1e-6
    s12       = 0.03e-5 / c
    s3        = 0.86e-5 / c
    I, ti     = I_sfel(I1, I2, I3, zi, zsep, s12, s3, ojit = 0.0)
    E, te, dx = tl.getE(I, ti, 2.5e-3, 0, 0)
    # Simulation parameters
    y0      = 800e-9
    tp      = 30e-15
    angle   = 15
    r0      = 2.5e-3
    ctype   = "GaP"
    d       = 50e-6
    method  = "cross"
    setup   = {"y0"     : y0,
               "tp"     : tp,
               "angle"  : angle,
               "fpath"  : fpath,
               "th"     : 0,
               "plot"   : False,
               "tilt"   : 0,
               "nslice" : 100,
               "r0"     : r0, 
               "ctype"  : ctype,
               "d"      : d,
               "method" : method,
               "plot"   : True
              }
    setup["tau"] = np.linspace(-2, 2, 1000)*1e-12
    #sig, t_sig, gamma, t_gamma = eos.E_signal(E, te, setup)
    #fig1 = plt.figure()
    #ax1  = fig1.gca()
    #E = 1e9 * np.exp(-te**2 / (2 * s3**2))
    #ax1.plot(te*1e15, E)
    #plt.show()
    #fig2 = plt.figure()
    #ax2  = fig2.gca()
    #ax2.plot(t_gamma * 1e15, gamma)
    #plt.show()
    dummy = eos.get_signal(0, setup)
main()