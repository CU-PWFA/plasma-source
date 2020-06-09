# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:18:48 2020
Main module for betatron radiation calculations. 
@author: keenan
"""

# Python modules
import numpy as np
from scipy.constants import c, mu_0, e

# Custom modules
import d2IMath as d2M
import trajectory as traj



fstring = ["x_noGain.npy", "y_noGain.npy", "gb_noGain.npy"]
# Load single particle trajectory
x_load  = np.load(fstring[0])
y_load  = np.load(fstring[1])
gb_load = np.load(fstring[2])


# Compute t, ct, tau, and z
L_plasma = 0.4
z        = np.linspace(0, L_plasma, len(x_load))
t        = z / c
# Compute tau (dtau is uniform with no energy gain)
dt = t[1] - t[0]
dtau = dt / gb_load
tau = np.zeros(len(dtau))
for i in range(1, len(dtau)):
    tau[i] = tau[i-1] + dtau[i-1]
# Create 4-position and interpolate trajectory
x = (c*t, x_load, y_load, z)
N_int = 1000
tau_int = np.linspace(tau[0], tau[-1], N_int)
x1n, x2n, x_int = traj.x_interp(x, tau, tau_int)
# Create 4-velocity and interpolate 
vx = np.diff(x[1]) / np.diff(tau)
vx= np.append(vx, vx[-1])

vy = np.diff(x[2]) / np.diff(tau)
vy= np.append(vy, vy[-1])

vz = np.diff(x[3]) / np.diff(tau)
vz= np.append(vz, vz[-1])

v = (vx, vy, vz)
v1n, v_int = traj.v_interp(v, tau, tau_int)
    

    

