#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 19:49:51 2017

@author: litos
"""

import numpy as np
from plasma_ramp import plasma_ramp
from propbeam import propbeam
from calc_Bmag import calc_Bmag

## define function
def fitdownramp(x,args):

    # fit params
    [waist,lp,P] = x

    # define plasma
    np0   = args[5]
    shape = args[6]
    z     = args[7]
    Lp    = args[8]
    pargs = [Lp,lp,P,1e-3]
    npl = plasma_ramp(np0,shape,z[0:len(z)-1],pargs,'down')

    # waist given w.r.t. z=L; need to convert
    abs_waist = Lp+waist

    # define beam
    [gb,eps,beta,alpha,gamma] = args[0:5]
    beam0 = args[0:5]

#    # propagate beam to desired waist location
#    iz_w = np.where(z<=abs_waist)[0]
#    z_w = z[iz_w]
#  
#    beam = propbeam(beam0,z_w,npl[0:len(z_w)-1])
#    [gb,eps,beta,alpha,gamma] = beam[len(beam)-1][:]
#
#    # calculate proximity to target
#    targ_beta  = args[9]
#    targ_alpha = args[10]
#    targ_gamma = args[11]
#    
#    Ttarg = [targ_beta,targ_alpha,targ_gamma]
#    Tbeam = [beta,alpha,gamma]
#    Bmag = calc_Bmag(Tbeam,Ttarg)
#    diff = Bmag - 1
#    d2 = diff**2

    # minimize gamma approach
    iz_w = np.where(z<=(Lp+3*lp))[0]
    z_w = z[iz_w]
    beam = propbeam(beam0,z_w,npl[0:len(z_w)-1])
    [gb,eps,beta,alpha,gamma] = beam[len(beam)-1][:]
    d2 = gamma

    print(lp)
    print(gamma)

    return d2