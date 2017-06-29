#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:38:19 2017

@author: mike
"""

import numpy as np
import nat_consts as nc
import mike_math as mm

def ebeam_propagation(ebeam,plasma,last_only=False):
    """function to propagate ebeam through plasma"""
    ebeam = prop_twiss(ebeam,plasma,last_only)
    ebeam = prop_parts(ebeam,plasma,last_only)
    return ebeam

def prop_twiss(ebeam,plasma,last_only=False):
    """propagate twiss parameters through entire plasma"""
    s    = plasma["s"]
    npl  = plasma["npl"]
    dgds = plasma["dgds"]
    
    """
    # buffer npl with zeros or
    # make zero everywhere if not given
    if len(npl)<len(s):
        dl = len(s)-len(npl)
        npl = np.append(npl,np.zeros((1,dl-1)))
        
    # buffer dgds with zeros or
    # make zero everywhere if not given
    if len(dgds)<len(s):
        dl = len(s)-len(dgds)
        dgds = np.append(dgds,np.zeros((1,dl-1)))    
    """
    
    Twiss = 
    
    
    
    # propagate beam ellipse
    for i in range(0,len(s)):
        if (last_only):
            Twiss = \
                prop_twiss_step(Twiss,s[i+1]-s[i],npl[i],dgds[i])
        else:
            Twiss[i+1][:] = \
                prop_twiss_step(Twiss[i][:],s[i+1]-s[i],npl[i],dgds[i])
    
    
    
    beam = np.zeros([len(s),len(beam0)])
    beam[0][:] = beam0

    # propagate beam ellipse
    for i in range(0,len(npl)):
        if (last_only):
            beam[:] = propbeamstep(beam[:],s[i+1]-s[i],npl[i],dgds[i])
        else:
            beam[i+1][:] = propbeamstep(beam[i][:],s[i+1]-s[i],npl[i],dgds[i])

    return beam

def prop_twiss_step(ebeam,ds,npl,dgds=0):
    """propagate Twiss parameters for a single step"""
    c = 3e8 # m/s, speed of light
    
    # define beam params
    [gb,eps,beta,alpha,gamma,dE,npart] = beam
    T = [beta,alpha,gamma]

    # calculate kb
    wp = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
    kp = wp/c # m^-1, plasma wave nutarget_betaer
    kb = kp/np.sqrt(2*gb) # m^-1, betatron wave number
    
    # beam phase space transfer matrix
    if kb>0: # if plasma density is non-zero
        R = [ [np.cos(kb*ds)    , np.sin(kb*ds)/kb], \
              [-kb*np.sin(kb*ds), np.cos(kb*ds)]  ]
    else: # treat like drift
        R = [ [1, ds], \
              [0,  1]  ]

    # Twiss parameter transfer matrix
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
    gb  = gb + dgb

    # reduce angle from energy gain/loss
    dgb = dgds*ds
    if dgb!=0:
        beta  = beta*(1+dgb/gb)
        gamma = gamma*(1-dgb/gb)
        alpha = np.sign(alpha)*np.sqrt(np.abs(beta*gamma-1))

    # return transported beam params
    beam = [gb,eps,beta,alpha,gamma,dE,npart]
    return beam
