#!/usr/bin/python

#===============
# propagate beam
#===============

## common libraries
from numpy import *

## define constants
c = 3e8 # m/s, speed of light

## define function
def propbeam(beam0,z,npl=[0]):

    beam = zeros([len(z),len(beam0)])
    beam[0][:] = beam0

    for i in range(0,len(npl)):
        beam[i+1][:] = propstep(beam[i][:],z[i+1]-z[i],npl[i])

    return beam

## single step propagation
def propstep(beam,dz,npl):

    # define beam params
    [gb,eps,beta,alpha,gamma] = beam
    T = [beta,alpha,gamma]

    # calculate kb
    wp = (5.64e4)*sqrt(npl) # rad/s, plasma ang. freq.
    kp = wp/c # m^-1, plasma wave nutarget_betaer
    kb = kp/sqrt(2*gb) # m^-1, betatron wave number

    # beam phase space transfer matrix
    if kb>0: # if plasma density is non-zero
        R = [ [cos(kb*dz)    , sin(kb*dz)/kb], \
              [-kb*sin(kb*dz), cos(kb*dz)   ]  ]
    else: # treat like drift
        R = [ [1, dz], \
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

    # return transported beam params
    beam = [gb,eps,beta,alpha,gamma]
    return beam
