#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 11:22:06 2017

@author: litos
"""

import numpy as np

## define function
def proppartslast(parts0,s,npl=[0],dgds=[0]):
    parts = propparts(parts0,s,npl,dgds)
    parts = parts[len(parts)-1][:]
    
    return parts

## propagate a group of individual particles
def propparts(parts0,s,npl=[0],dgds=[0]):
    
    # buffer npl with zeros or
    # make zero everywhere if not given
    if len(npl)<len(s):
        dn = len(s)-len(npl)
        npl = np.append(npl,np.zeros((1,dn-1)))
        
    # buffer dgds with zeros or
    # make zero everywhere if not given
    if len(dgds)<len(s):
        dn = len(s)-len(dgds)
        dgds = np.append(dgds,np.zeros((1,dn-1)))
        
    # step over particles
    npart = parts0.shape[0]
    
    parts = np.zeros((len(s),npart,6))
    
    for ipart in range(0,npart):
        parts[0][ipart][:] = parts0[ipart][:]
        # propagate individual particles
        for istep in range(0,len(npl)):
            parts[istep+1][ipart][:] =\
                proppartstep(parts[istep][ipart][:],\
                             s[istep+1]-s[istep],npl[istep],dgds[istep])
                
    return parts

## single step propagation
def proppartstep(part,ds,npl,dgds):

    c = 3e8 # m/s, speed of light
    
    # define particle phase space params
    [x,xp,y,yp,z,gb] = part

    # calculate kb
    wp = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
    kp = wp/c # m^-1, plasma wave nutarget_betaer
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
    dgb = dgds*ds
    gb = gb + dgb
    
    # reduce angle from energy gain/loss
    if dgb!=0:
        R = [ [1,        0], \
             [0, 1-dgb/gb]  ]
        X = np.dot(R,X)
        Y = np.dot(R,Y)
    
    [x,xp] = X
    [y,yp] = Y

    # return transported beam params
    part = [x,xp,y,yp,z,gb]
    return part
