#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 21:34:52 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)
import sys

#%%
p= np.linspace(0.5, 12, 30)
FirstSpacing= p*0.34+2
AGWidth= p*0.001+4
#%%
NorFS= FirstSpacing/np.amax(FirstSpacing)
NorAGW= AGWidth/np.amax(AGWidth)
#%%
plt.plot(p, AGWidth, '.')
plt.text(5, 80, "FAKE PLOT", size= 40, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))
plt.xlabel('Gas Pressure')

#%%
Pinmbar=14.2
PinPa= Pinmbar*100
kB=1.3807e-23
T=293
NoverV= PinPa/kB/T
DenIncm3=NoverV/(100**3)
#%%
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(p, NorFS, '.', label= 'First Peak Spacing')
ax.plot(p, NorAGW, '.', label= 'Afterglow Width')
ax.set_xlabel('Pressure (mbar)')
ax.set_ylabel('Normalized Parameter')
plt.text(0.5, 0.8, "FAKE PLOT", size= 40, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))


def PtoN(Pinmbar):
    PinPa= Pinmbar*100
    NoverV= PinPa/kB/T
    DenIncm3=NoverV/(100**3)
    return DenIncm3

def NtoP(NoverV):
    return NoverV*kB*T

secax = ax.secondary_xaxis('top', functions=(PtoN, NtoP))
secax.set_xlabel('Density ($cm^{-3}$)')
plt.legend()
plt.show()


#%%
#%%
I= np.linspace(300, 500, 30)
FirstSpacing= I*0.0005+2
AGWidth= p*3+4

#%%
NorFS= FirstSpacing/np.amax(FirstSpacing)
NorAGW= AGWidth/np.amax(AGWidth)

#%%
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(I, NorFS, '.', label= 'First Peak Spacing')
ax.plot(I, NorAGW, '.', label= 'Afterglow Width')
ax.set_xlabel('Laser Energy (mJ)')
ax.set_ylabel('Normalized Parameter')
plt.text(300, 0.78, "FAKE PLOT", size= 40, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))


def ItoW(I):
    return I*0.2

def WtoI(W):
    return W/0.2 

secax = ax.secondary_xaxis('top', functions=(ItoW, WtoI))
secax.set_xlabel('Calculated Plasma Width ($\mu m$)')
plt.legend()
plt.show()


#%%
#%%
t= np.linspace(0, 50, 30)
n_exp= 1e17*np.exp(-(t**2/15**2)**3)
n_th= 1e17*np.exp(-(t**2/14**2)**2.5)
#%%
plt.plot(t, n_exp, '.', label= 'Experiment')
plt.plot(t, n_th, '-.', label= 'Simulation')
plt.xlabel('time (ns)')
plt.ylabel('Central Density')
plt.text(20, 0.6e17, "FAKE PLOT", size= 35, bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))
plt.legend()
#%%







#%%
#%%
x = p
y1 = FirstSpacing
y2 = AGWidth
for line in sys.stdin:
  vals = line.split(',')
  x.append(float(vals[0]))
  y1.append(float(vals[1]))
  y2.append(float(vals[2]))

# Plot y1 vs x in blue on the left vertical axis.
plt.xlabel("x")
plt.ylabel("Blue", color="b")
plt.tick_params(axis="y", labelcolor="b")
plt.plot(x, y1, "b-", linewidth=2)

# Plot y2 vs x in red on the right vertical axis.
plt.twinx()
plt.ylabel("Red", color="r")
plt.tick_params(axis="y", labelcolor="r")
plt.plot(x, y2, "r-", linewidth=2)


