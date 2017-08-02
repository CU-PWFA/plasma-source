#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:38:19 2017

@author: mike
"""

import numpy as np
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import nat_consts as nc
import particle_beam as pb

def prop_ebeam_drift(ebeam,s,last_only=False):
    """function to propagate ebeam through a drift"""
    vacuum  = defaultdict(dict)
    vacuum["s"]    = s
    vacuum["npl"]  = np.zeros(len(s))
    vacuum["dgds"] = np.zeros(len(s))
    prop_ebeam_plasma(ebeam,vacuum,last_only)
    return

def prop_ebeam_plasma(ebeam,plasma,last_only=False):
    """function to propagate ebeam through plasma"""
    s     = plasma["s"]
    npl   = plasma["npl"]
    dgds  = plasma["dgds"]
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    
    # propagate beam
    prop_twiss_plasma(twiss,s,npl,dgds)
    prop_parts_plasma(parts,s,npl,dgds)

    nstep = len(s)
    for i in range(1,nstep):
        step = len(ebeam)
        pb.append_ebeam_step(ebeam,step,s[i],twiss[i],parts[i])

    return

def prop_twiss_plasma(twiss,s,npl,dgds):
    nstep = len(s)
    i_twiss = twiss[len(twiss)-1].copy()
    for i in range(1,nstep):
        prop_twiss_plasma_step(i_twiss,s[i]-s[i-1],npl[i-1],dgds[i-1])
        twiss[len(twiss)] = i_twiss.copy()
    return

def prop_twiss_plasma_step(i_twiss,ds=0,npl=0,dgds=0):
    """propagate Twiss parameters through plasma for a single step"""
    beta  = i_twiss["beta"]
    alpha = i_twiss["alpha"]
    gamma = i_twiss["gamma"]
    eps   = i_twiss["eps"]
    gbC   = i_twiss["gbC"]
    dgb   = i_twiss["dgb"]
    dz    = i_twiss["dz"]
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
    i_twiss["beta"]  = beta
    i_twiss["alpha"] = alpha
    i_twiss["gamma"] = gamma
    i_twiss["eps"]   = eps
    i_twiss["gbC"]   = gbC
    i_twiss["dgb"]   = dgb
    i_twiss["dz"]    = dz
    
    return

def prop_parts_plasma(parts,s,npl,dgds):
    nstep = len(s)
    npart = parts[len(parts)-1]["npart"]
    
    if npart>0:
        phase6D0 = []
        for i in range(0,npart):
            x  = parts[len(parts)-1]["x"][i]
            xp = parts[len(parts)-1]["xp"][i]
            y  = parts[len(parts)-1]["y"][i]
            yp = parts[len(parts)-1]["yp"][i]
            z  = parts[len(parts)-1]["z"][i]
            gb = parts[len(parts)-1]["gb"][i]
            phase6D0 = np.append(phase6D0,[x,xp,y,yp,z,gb])
   
        # propagate individual particles in parallel
        num_cores = multiprocessing.cpu_count()
        phase6D = Parallel(n_jobs=num_cores)\
            (delayed(prop_part_plasma)(phase6D0[i*6:(i+1)*6],s,npl,dgds)\
             for i in range(0,npart))
        phase6D = np.reshape(phase6D,[npart,6*nstep])
    
        # loop over steps to create parts dictionary
        for i in range(1,nstep):
            j = len(parts)
            parts[j]["x"]  = phase6D[:,i*6+0]
            parts[j]["xp"] = phase6D[:,i*6+1]
            parts[j]["y"]  = phase6D[:,i*6+2]
            parts[j]["yp"] = phase6D[:,i*6+3]
            parts[j]["z"]  = phase6D[:,i*6+4]
            parts[j]["gb"] = phase6D[:,i*6+5]
            parts[j]["npart"] = parts[j-1]["npart"]
            parts[j]["dist"]  = parts[j-1]["dist"]

    else:
        for i in range(1,nstep):
            j = len(parts)
            parts[j] = parts[j-1].copy()

    return

def prop_part_plasma(phase6D0,s,npl,dgds):
    """propagate single macro particle through multiple steps"""
    """needs to return unique object for effective use in parallel"""
    i_phase6D = phase6D0.copy()
    phase6D   = []
    phase6D   = np.append(phase6D,i_phase6D)
    for i in range(1,len(s)):
        i_phase6D = prop_part_plasma_step(i_phase6D,s[i]-s[i-1],npl[i],dgds[i])
        phase6D = np.append(phase6D,i_phase6D)
    return phase6D

def prop_part_plasma_step(i_phase6D,ds=0,npl=0,dgds=0):
    """propagate single macro particle through plasma for a single step"""
    """needs to return unique object for effective use in parallel"""
    [x,xp,y,yp,z,gb] = i_phase6D

    # calculate kb
    wp = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
    kp = wp/nc.c # m^-1, plasma wave number
    kb = kp/np.sqrt(2*gb) # m^-1, betatron wave number
 
    # beam phase space transfer matrix
    if kb>0: # if plasma density is non-zero
        R = [ [np.cos(kb*ds)    , np.sin(kb*ds)/kb], \
              [-kb*np.sin(kb*ds), np.cos(kb*ds)   ]  ]
    else: # treat like drift
        R = [ [1, ds], \
              [0, 1]  ]

    # perform beam transport
    X = [x,xp]
    Y = [y,yp]
    
    X = np.dot(R,X)
    Y = np.dot(R,Y)

    # add energy gain/loss
    Dgb = dgds*ds
    gb  = gb + Dgb
    
    # reduce angle from energy gain/loss
    if Dgb!=0:
        R = [ [1,        0], \
              [0, 1-Dgb/gb]  ]
        X = np.dot(R,X)
        Y = np.dot(R,Y)

    [x,xp] = X
    [y,yp] = Y
        
    i_phase6D = [x,xp,y,yp,z,gb]
    return i_phase6D