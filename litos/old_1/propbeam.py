#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 28 16:44:55 2017

@author: litos
"""

import numpy as np

## define function
def propbeamlast(beam0,s,npl=[0],dgds=[0]):
    beam = propbeam(beam0,s,npl,dgds)
    beam = beam[len(beam)-1][:]
    
    return beam
    
## define function
def propbeam(beam0,s,npl=[0],dgds=[0]):

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
    
    beam = np.zeros([len(s),len(beam0)])
    beam[0][:] = beam0

    # propagate beam ellipse
    for i in range(0,len(npl)):
        beam[i+1][:] = propbeamstep(beam[i][:],s[i+1]-s[i],npl[i],dgds[i])

    return beam

## single step propagation
def propbeamstep(beam,ds,npl,dgds=0):

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
