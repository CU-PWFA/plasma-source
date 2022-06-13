# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:41:48 2020

@author: keenan
"""
import numpy as np

def read_data(path, nDumps, N):
    """
    Function to read in the phase space from a WARGSim run
    
    Parameters:
    -----------
    path : string
        The path to the dump files
    nDumps: int
        The number of dumps to read
    N : int
        The number of particles in the simulation
    Returns:
    --------
    x : array_like
        Particle x trajectory
    y : array_like
        Particle y trajectory
    gb : array_like
         Particle lorentz factor
    """
    x      = np.zeros((nDumps, N))
    xp     = np.zeros((nDumps, N))
    y      = np.zeros((nDumps, N))
    yp     = np.zeros((nDumps, N))
    z      = np.zeros((nDumps, N))
    gb     = np.zeros((nDumps, N))    
    
    for i in range(nDumps):
        ptcls   = np.load(path + "_ptcls_" + str(i) + ".npy")
        x[i,:]  = ptcls[:, 0]
        xp[i,:] = ptcls[:, 1]
        y[i,:]  = ptcls[:, 2]
        yp[i,:] = ptcls[:, 3]
        z[i, :] = ptcls[:, 4]
        gb[i,:] = ptcls[:, 5] 
        
        
    return x, xp, y, yp, z, gb
