#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 17:02:43 2017

@author: litos
"""

from plasma_ramp import plasma_ramp
from propbeam import propbeamlast
from calc_M import calc_M

## define function
def fitupramp(x,args):

    # fit params
    [waist,lp,P] = x

    # define plasma
    npl0  = args[5]
    shape = args[6]
    s     = args[7]
    Lp    = args[8]
    pargs = [Lp,lp,P,1e-3]
    npl = plasma_ramp(npl0,shape,s[0:len(s)-1],pargs)

    # "waist" given w.r.t. z=L; need to calc. absolute position
    abs_waist = Lp+waist

    # define beam
    [gb,eps,beta,alpha,gamma] = args[0:5]
    beam0 = [gb,eps,beta,alpha,gamma,0,0]
    beam0 = propbeamlast(beam0,[0,-1*abs_waist])

    # propagate beam
    beam = propbeamlast(beam0,s,npl)
    [gb,eps,beta,alpha,gamma,dE,npart] = beam

    # calculate proximity to target
    targ_beta  = args[9]
    targ_alpha = args[10]
    targ_gamma = args[11]
    
    Ttarg = [targ_beta,targ_alpha,targ_gamma]
    Tbeam = [beta,alpha,gamma]
    M = calc_M(Tbeam,Ttarg)
    diff = M - 1
    d2 = diff**2

    return d2