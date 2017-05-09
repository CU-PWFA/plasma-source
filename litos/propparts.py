#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 11:22:06 2017

@author: litos
"""

from numpy import *

## define function
def proppartslast(parts0,s,npl=[0],dgds=0):
    parts = propparts(parts0,s,npl,dgds)
    parts = parts[len(parts)-1][:]
    
    return parts

## propagate a group of individual particles
def propparts(parts0,s,npl=[0],dgds=0):
    
    # buffer npl with zeros or
    # make zero everywhere if not given
    if len(npl)<len(s):
        dl = len(s)-len(npl)
        npl = append(npl,zeros((1,dl-1)))
        
    # step over particles
    npart = parts0.shape[0]
    
    parts = zeros((len(s),npart,6))
    
    for ipart in range(0,npart):
        parts[0][ipart][:] = parts0[ipart][:]
        # propagate individual particles
        for istep in range(0,len(npl)):
            parts[istep+1][ipart][:] =\
                proppartstep(parts[istep][ipart][:],\
                             s[istep+1]-s[istep],npl[istep],dgds)
                
    return parts

## single step propagation
def proppartstep(part,ds,npl,dgds=0):

    c = 3e8 # m/s, speed of light
    
    # define particle phase space params
    [x,xp,y,yp,z,gb] = part

    # calculate kb
    wp = (5.64e4)*sqrt(npl) # rad/s, plasma ang. freq.
    kp = wp/c # m^-1, plasma wave nutarget_betaer
    kb = kp/sqrt(2*gb) # m^-1, betatron wave number
 
    # beam phase space transfer matrix
    if kb>0: # if plasma density is non-zero
        R = [ [cos(kb*ds)    , sin(kb*ds)/kb], \
              [-kb*sin(kb*ds), cos(kb*ds)   ]  ]
    else: # treat like drift
        R = [ [1, ds], \
              [0, 1]  ]

    # perform beam transport
    X = [x,xp]
    Y = [y,yp]
    
    X = dot(R,X)
    Y = dot(R,Y)
    
    [x,xp] = X
    [y,yp] = Y

    # add energy
    gb = gb + dgds*ds

    # return transported beam params
    part = [x,xp,y,yp,z,gb]
    return part
