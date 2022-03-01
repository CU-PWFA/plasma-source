#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 13:42:32 2020

@author: valentinalee
"""


"""
This is a file to establish the initial laser pulse parameters for ionization simulation.
*Species: Helium*
*Linear polarized*
*All the units are in SI unit*
*only plotting unit are as labeled for comprehensing purpose*
*Wavelength: 800nm*

The spatial profile fitting equation is: 19.777e18*exp(-(r^2/(2*(29.417*1e-6)^2))^1.076)
The temporal profile fitting equation is: 19.777e18*exp((-1/2)*(t/12.8e-15)^2)
which is a gaussian of a FWHM=30fs


Both spatial and temporal are plotted in this document for visualization
"""

import numpy as np
import matplotlib.pyplot as plt

I0= 19.777e18
r= np.linspace(-200e-6, 200e-6, 1000)
Ir= np.exp(-(r**2/(2*(29.417*1e-6)**2))**1.076)

t= np.linspace(-100e-15, 100e-15, 1000)
It= np.exp((-1/2)*(t/12.8e-15)**2)

plt.figure(1)
plt.plot(r*1e6, Ir*I0*1e-4)
plt.xlabel('r ($\mu$m)')
plt.ylabel('I (W/cm^2)')

plt.figure(2)
plt.plot(t*1e15, It*I0*1e-4)
plt.xlabel('Time (fs)')
plt.ylabel('I (W/cm^2)')

IR, IT= np.meshgrid(Ir, It)
Irt= IR*IT

plt.figure(3)
plt.pcolormesh(t*1e15, r*1e6, Irt*I0*1e-4)
plt.xlabel('Time (fs)')
plt.ylabel('r ($\mu$m)')
plt.colorbar(label='I (W/cm^2)')


