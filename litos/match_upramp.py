#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 17:02:43 2017

@author: litos
"""

import numpy as np
from plasma_ramp import plasma_ramp
from mc_beam import mc_beam
from propparts import propparts
from propparts import proppartslast
from calc_twiss import calc_twiss
from calc_M import calc_M

## define function
def match_upramp(x,args):

    # constants
    c  = 3e8 # m/s
    me = 0.511e-3 # GeV

    # fraction of particles to use
    frac  = 1.0
    
    # fit params
    [waist,hw_up] = x

    # define particles
    parts0 = args[0]

    # define plasma
    npl0  = args[1]
    dgds0 = args[2]
    shape = args[3]
    s     = args[4]
    L_up  = args[5]
    
    # make plasma up ramp
    s_up = s[(s<=L_up)]
    pargs = [L_up,hw_up]
    [npl_up,dgds_up] = plasma_ramp(npl0,shape,s_up[0:len(s_up)-1],pargs,'up',dgds0)

    # make bulk flat-top plasma
    s_ft = s[(s>L_up)]
    npl_ft  = npl0*np.ones(len(s_ft)-1)
    dgds_ft  = dgds0*np.ones(len(s_ft)-1)
    
    # combine plasma sections
    npl  = np.hstack((npl_up,npl_ft))
    dgds = np.hstack((dgds_up,dgds_ft))

    # "waist" given w.r.t. s=L_up; need to calc. absolute waist
    abs_waist = L_up+waist

    # set beam waist
    s_w    = abs_waist # m
    parts0 = proppartslast(parts0,[0,-s_w])

    # propagate particles
    parts = propparts(parts0,s,npl,dgds)
    par_x = np.zeros(4)
    par_y = np.zeros(4)
    navg  = 20
    for i in range(0,navg):
        iparts = parts[len(parts)-navg+i]
        [rms_parts,rms_X,rms_Y] = calc_twiss(iparts,frac)
        par_x += rms_X
        par_y += rms_Y
    par_avg = (par_x+par_y)/(2*navg)
    eps_avg = par_avg[0]

#    Tbeam = par_avg[1:]

#    beta_avg = par_avg[1]
#    gb = np.mean(parts[len(parts)-1,:,5])
#    wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
#    kp0 = wp0/c # m^-1, plasma wave number
#    kb = kp0/np.sqrt(2*gb) # m^-1, betatron wave number
#    dbeta = beta_avg-1/kb
#    dbeta2 = dbeta**2
    
#    parts = proppartslast(parts0,s,npl,dgds)
#    [rms_parts,rms_X,rms_Y] = calc_twiss(parts,frac)
#    Tx = np.array(rms_X[1:])
#    Ty = np.array(rms_Y[1:])
#    Tbeam = (Tx+Ty)/2
#    gb = np.mean(parts[:,5])
#    
#    wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
#    kp0 = wp0/c # m^-1, plasma wave number
#    kb = kp0/np.sqrt(2*gb) # m^-1, betatron wave number
#    Tmatch = [1/kb,0,kb]
#    M = calc_M(Tbeam,Tmatch)
     

#    parts = proppartslast(parts0,s,npl,dgds)
#    frac = 1.00
#    [rms_parts,rms_X,rms_Y] = calc_twiss(parts,frac)
#    [eps_x,beta_x,alpha_x,gamma_x] = rms_X
#    [eps_y,beta_y,alpha_y,gamma_y] = rms_Y
#    eps_avg = (eps_x+eps_y)/2

    print(eps_avg)

    return eps_avg