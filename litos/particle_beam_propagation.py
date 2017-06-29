#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:38:19 2017

@author: mike
"""

import numpy as np
from collections import defaultdict
import nat_consts as nc
import particle_beam as pb

def prop_ebeam_drift(ebeam,s,last_only=False):
    """function to propagate ebeam through a drift"""
    plasma = defaultdict(dict)
    plasma["s"]    = s
    plasma["npl"]  = np.zeros(len(s))
    plasma["dgds"] = np.zeros(len(s))
    ebeam = prop_ebeam_plasma(ebeam,plasma,last_only)

    return ebeam

def prop_ebeam_plasma(ebeam,plasma,last_only=False):
    """function to propagate ebeam through plasma"""
    s    = plasma["s"]
    npl  = plasma["npl"]
    dgds = plasma["dgds"]
     
    twiss = ebeam[len(ebeam)-1]["twiss"]
    parts = ebeam[len(ebeam)-1]["parts"]
    
    # propagate beam
    for i in range(1,len(s)):      
        twiss = prop_twiss_plasma_step(twiss,s[i]-s[i-1],npl[i-1],dgds[i-1])
        parts = prop_parts_plasma_step(parts,s[i]-s[i-1],npl[i-1],dgds[i-1])

        if (last_only):
            continue
        else:
            step  = len(ebeam)
            ebeam = pb.append_ebeam_step(ebeam,step,s[i],twiss,parts)

    return ebeam

def prop_twiss_plasma_step(twiss,ds,npl=0,dgds=0):
    """propagate Twiss parameters through plasma for a single step"""
    beta  = twiss["beta"]
    alpha = twiss["alpha"]
    gamma = twiss["gamma"]
    eps   = twiss["eps"]
    gbC   = twiss["gbC"]
    dgb   = twiss["dgb"]
    dz    = twiss["dz"]
    T     = [beta,alpha,gamma]

    # calculate kb
    wp = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
    kp = wp/nc.c # m^-1, plasma wave nutarget_betaer
    kb = kp/np.sqrt(2*gbC) # m^-1, betatron wave number
    
    # create beam phase space transfer matrix
    if kb>0: # if plasma density is non-zero
        R = [ [np.cos(kb*ds)    , np.sin(kb*ds)/kb], \
              [-kb*np.sin(kb*ds), np.cos(kb*ds)]  ]
    else: # treat like drift
        R = [ [1, ds], \
              [0,  1]  ]

    # create Twiss parameter transfer matrix
    RTwiss = [ [  R[0][0]**2, -2*R[0][0]*R[0][1], R[0][1]**2 ], \
               [ -R[0][0]*R[1][0], \
                  R[0][1]*R[1][0]+R[0][0]*R[1][1], \
                 -R[0][1]*R[1][1] ], \
               [  R[1][0]**2, -2*R[1][0]*R[1][1], R[1][1]**2 ] ]

    # perform beam transport
    T = np.dot(RTwiss,T)
    [beta,alpha,gamma] = T

    # add energy gain/loss
    dgb = dgds*ds
    gbC = gbC + dgb

    # reduce angle from energy gain/loss
    dgb = dgds*ds
    if dgb!=0:
        beta  = beta*(1+dgb/gbC)
        gamma = gamma*(1-dgb/gbC)
        alpha = np.sign(alpha)*np.sqrt(np.abs(beta*gamma-1))

    # return transported beam params
    twiss["beta"]  = beta
    twiss["alpha"] = alpha
    twiss["gamma"] = gamma
    twiss["eps"]   = eps
    twiss["gbC"]   = gbC
    twiss["dgb"]   = dgb
    twiss["dz"]    = dz
    
    return twiss

def prop_parts_plasma_step(parts,ds,npl=0,dgds=0):
    """propagate macro particles through plasma for a single step"""
    x     = parts["x"]
    xp    = parts["xp"]
    y     = parts["y"]
    yp    = parts["yp"]
    z     = parts["z"]
    gb    = parts["gb"]
    npart = parts["npart"]
    
    # loop over particles
    for i in range(0,npart):
    
        # calculate kb
        wp = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
        kp = wp/nc.c # m^-1, plasma wave nutarget_betaer
        kb = kp/np.sqrt(2*gb[i]) # m^-1, betatron wave number
     
        # beam phase space transfer matrix
        if kb>0: # if plasma density is non-zero
            R = [ [np.cos(kb*ds)    , np.sin(kb*ds)/kb], \
                  [-kb*np.sin(kb*ds), np.cos(kb*ds)   ]  ]
        else: # treat like drift
            R = [ [1, ds], \
                  [0, 1]  ]
    
        # perform beam transport
        X = [x[i],xp[i]]
        Y = [y[i],yp[i]]
        
        X = np.dot(R,X)
        Y = np.dot(R,Y)
    
        # add energy gain/loss
        Dgb = dgds*ds
        gb[i]  = gb[i] + Dgb
        
        # reduce angle from energy gain/loss
        if Dgb!=0:
            R = [ [1,        0], \
                  [0, 1-Dgb/gb[i]]  ]
            X = np.dot(R,X)
            Y = np.dot(R,Y)
        
        [x[i],xp[i]] = X
        [y[i],yp[i]] = Y

    # return transported beam params
    parts["x"]  = x
    parts["xp"] = xp
    parts["y"]  = y
    parts["yp"] = yp
    parts["z"]  = z
    parts["gb"] = gb
    parts["npart"] = npart
 
    return parts