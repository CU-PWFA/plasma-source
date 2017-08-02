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

def prop_ebeam_drift(ebeam0,s,last_only=False):
    """function to propagate ebeam through a drift"""
    plasma = defaultdict(dict)
    plasma["s"]    = s
    plasma["npl"]  = np.zeros(len(s))
    plasma["dgds"] = np.zeros(len(s))
    ebeam = prop_ebeam_plasma(ebeam0,plasma,last_only)
    return ebeam

def prop_ebeam_plasma(ebeam,plasma,last_only=False):
    """function to propagate ebeam through plasma"""
    s     = plasma["s"]
    npl   = plasma["npl"]
    dgds  = plasma["dgds"]
    twiss0 = pb.get_twiss(ebeam,len(ebeam)-1)
    parts0 = pb.get_parts(ebeam,len(ebeam)-1)
    
    # propagate beam
    twiss = prop_twiss_plasma(twiss0,s,npl,dgds)
    parts = prop_parts_plasma(parts0,s,npl,dgds)
    
    nstep = len(s)
    for i in range(1,nstep):
        step = len(ebeam)
        ebeam = pb.append_ebeam_step(ebeam,step,s[i],twiss[i],parts[i])

    return ebeam

def prop_twiss_plasma(twiss0,s,npl,dgds):
    twiss = defaultdict(dict)
    itwiss = twiss0
    for i in range(0,len(s)):
        if i>0:
            itwiss = prop_twiss_plasma_step(itwiss,s[i]-s[i-1],\
                                            npl[i-1],dgds[i-1])
        twiss[i]["beta"]  = itwiss["beta"]
        twiss[i]["alpha"] = itwiss["alpha"]
        twiss[i]["gamma"] = itwiss["gamma"]
        twiss[i]["eps"]   = itwiss["eps"]
        twiss[i]["gbC"]   = itwiss["gbC"]
        twiss[i]["dgb"]   = itwiss["dgb"]
        twiss[i]["dz"]    = itwiss["dz"]
    return twiss

def prop_twiss_plasma_step(twiss,ds=0,npl=0,dgds=0):
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

def prop_parts_plasma(parts0,s,npl,dgds):
    parts = defaultdict(dict)
    
    nstep = len(s)
    npart = parts0["npart"]
    
    if npart>0:
        parts[0] = parts0
        phase6D0 = []
        for i in range(0,npart):
            x  = parts[0]["x"][i]
            xp = parts[0]["xp"][i]
            y  = parts[0]["y"][i]
            yp = parts[0]["yp"][i]
            z  = parts[0]["z"][i]
            gb = parts[0]["gb"][i]
            phase6D0 = np.append(phase6D0,[x,xp,y,yp,z,gb])
    
        # propagate individual particles in parallel
        num_cores = multiprocessing.cpu_count()
        phase6D = Parallel(n_jobs=num_cores)\
            (delayed(prop_part_plasma)(phase6D0[i*6:(i+1)*6],s,npl,dgds)\
             for i in range(0,npart))
        phase6D = np.reshape(phase6D,[npart,6*nstep])
    
        # loop over steps to create parts dictionary
        for i in range(0,nstep):
            parts[i]["x"]  = phase6D[:,i*6+0]
            parts[i]["xp"] = phase6D[:,i*6+1]
            parts[i]["y"]  = phase6D[:,i*6+2]
            parts[i]["yp"] = phase6D[:,i*6+3]
            parts[i]["z"]  = phase6D[:,i*6+4]
            parts[i]["gb"] = phase6D[:,i*6+5]
            parts[i]["npart"] = parts[0]["npart"]
            parts[i]["dist"]  = parts[0]["dist"]

    else:
        for i in range(0,nstep):
            parts[i] = parts0

    return parts

def prop_part_plasma(phase6D0,s,npl,dgds):
    phase6D = phase6D0
    # loop over steps
    for i in range(1,len(s)):
        iphase6D = prop_part_plasma_step(phase6D[(i-1)*6:(i-1+1)*6],\
                                         s[i]-s[i-1],npl[i],dgds[i])
        phase6D = np.append(phase6D,iphase6D)
    return phase6D

def prop_part_plasma_step(phase6D0,ds=0,npl=0,dgds=0):
    """propagate single macro particle through plasma for a single step"""
    [x,xp,y,yp,z,gb] = phase6D0

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
        
    phase6D = [x,xp,y,yp,z,gb]
    
    return phase6D



#def prop_parts_plasma_step(parts,ds=0,npl=0,dgds=0):
#    """propagate all macro particles through plasma for a single step"""
#
##    # propagate individual particles in parallel
##    npart = parts["npart"]
###    num_cores = multiprocessing.cpu_count()
##    num_cores = 1
##    phase6D = Parallel(n_jobs=num_cores)\
##        (delayed(prop_part_plasma_step)(parts,ds,npl,dgds,i)\
##         for i in range(0,npart))
##    phase6D = np.reshape(phase6D,[npart,6])
#
#
#    # propagate individual particles in parallel
#    npart = parts["npart"]
#    nprop = 100
#    num_cores = multiprocessing.cpu_count()
##    num_cores = 1
#    phase6D = Parallel(n_jobs=num_cores)\
#        (delayed(prop_part_plasma_step)(parts,ds,npl,dgds,nprop,i0)\
#         for i0 in range(0,int(npart/nprop)))
#    phase6D = np.reshape(phase6D,[npart,6])
#
#
#
#    # return transported beam params
#    parts["x"]  = phase6D[:,0]
#    parts["xp"] = phase6D[:,1]
#    parts["y"]  = phase6D[:,2]
#    parts["yp"] = phase6D[:,3]
#    parts["z"]  = phase6D[:,4]
#    parts["gb"] = phase6D[:,5]
#    return parts






#def prop_part_plasma_step(parts,ds=0,npl=0,dgds=0,nprop=1,i0=0):
#    """propagate subset of macro particles through plasma for a single step"""
#    phase6D = []
#    x  = np.zeros(nprop)
#    xp = np.zeros(nprop)
#    y  = np.zeros(nprop)
#    yp = np.zeros(nprop)
#    z  = np.zeros(nprop)
#    gb = np.zeros(nprop)
#    for i in range(0,nprop):
#        x[i]  = parts["x"][i0+i]
#        xp[i] = parts["xp"][i0+i]
#        y[i]  = parts["y"][i0+i]
#        yp[i] = parts["yp"][i0+i]
#        z[i]  = parts["z"][i0+i]
#        gb[i] = parts["gb"][i0+i]
#    
#        # calculate kb
#        wp = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
#        kp = wp/nc.c # m^-1, plasma wave nutarget_betaer
#        kb = kp/np.sqrt(2*gb[i]) # m^-1, betatron wave number
#     
#        # beam phase space transfer matrix
#        if kb>0: # if plasma density is non-zero
#            R = [ [np.cos(kb*ds)    , np.sin(kb*ds)/kb], \
#                  [-kb*np.sin(kb*ds), np.cos(kb*ds)   ]  ]
#        else: # treat like drift
#            R = [ [1, ds], \
#                  [0, 1]  ]
#    
#        # perform beam transport
#        X = [x[i],xp[i]]
#        Y = [y[i],yp[i]]
#        
#        X = np.dot(R,X)
#        Y = np.dot(R,Y)
#    
#        # add energy gain/loss
#        Dgb   = dgds*ds
#        gb[i] = gb[i] + Dgb
#        
#        # reduce angle from energy gain/loss
#        if Dgb!=0:
#            R = [ [1,           0], \
#                  [0, 1-Dgb/gb[i]]  ]
#            X = np.dot(R,X)
#            Y = np.dot(R,Y)
#    
#        [x[i],xp[i]] = X
#        [y[i],yp[i]] = Y
#
#        phase6D = np.append(phase6D,[x[i],xp[i],y[i],yp[i],z[i],gb[i]])
#
#    return phase6D

