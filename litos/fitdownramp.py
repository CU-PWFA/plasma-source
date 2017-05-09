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

    # define initial beam
    [gb,eps,beta,alpha,gamma] = args[0:5]
    beam0 = args[0:5]

    # propagate virtual beam from virtual waist
    # to end of z through vacuum for matching targets
    v_beta  = args[9]
    v_alpha = args[10]
    v_gamma = args[11]
    vbeam0 = [gb,eps,v_beta,v_alpha,v_gamma]
    
    iz_w = np.where(z>=waist)[0]
    z_w = z[iz_w]
    
    vbeam = propbeam(vbeam0,z_w)
    [vgb,veps,targ_beta,targ_alpha,targ_gamma] = \
        vbeam[len(vbeam)-1][:]
    
    # propagate real beam through plasma
    # from z=0 to end
    beam = propbeam(beam0,z,npl[0:len(z)-1])
    [gb,eps,beta,alpha,gamma] = beam[len(beam)-1][:]
    
    # calculate proximity to target
    Ttarg = [targ_beta,targ_alpha,targ_gamma]
    Tbeam = [beta,alpha,gamma]
    
#    print(Ttarg)
#    print(Tbeam)
    
    Bmag = calc_Bmag(Tbeam,Ttarg)
    diff = Bmag - 1
    d2 = diff**2
    
#    print(d2)
 
    return d2