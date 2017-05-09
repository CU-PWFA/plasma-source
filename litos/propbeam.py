#!/usr/bin/python

#===============
# propagate beam
#===============

## common libraries
from numpy import *

## define constants
c = 3e8 # m/s, speed of light

## define function
def propbeamlast(beam0,s,npl=[0],dgds=0):
    beam = propbeam(beam0,s,npl,dgds)
    beam = beam[len(beam)-1][:]
    
    return beam
    
## define function
def propbeam(beam0,s,npl=[0],dgds=0):

    beam = zeros([len(s),len(beam0)])
    beam[0][:] = beam0

    if len(npl)<len(s):
        dl = len(s)-len(npl)
        npl = append(npl,zeros((1,dl-1)))

    for i in range(0,len(npl)):
        beam[i+1][:] = propstep(beam[i][:],s[i+1]-s[i],npl[i],dgds)

    return beam

## single step propagation
def propstep(beam,ds,npl,dgds=0):

    # define beam params
    [gb,eps,beta,alpha,gamma,dE,npart] = beam
    T = [beta,alpha,gamma]

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

    # Twiss parameter transfer matrix
    RTwiss = [ [ R[0][0]**2, -2*R[0][0]*R[0][1], R[0][1]**2 ], \
               [ -R[0][0]*R[1][0], \
                R[0][1]*R[1][0]+R[0][0]*R[1][1], \
                -R[0][1]*R[1][1] ], \
               [ R[1][0]**2, -2*R[1][0]*R[1][1], R[1][1]**2 ] ]

    # perform beam transport
    T = dot(RTwiss,T)
    [beta,alpha,gamma] = T
    
    # add energy
    gb = gb + dgds*ds

    # return transported beam params
    beam = [gb,eps,beta,alpha,gamma,dE,npart]
    return beam
