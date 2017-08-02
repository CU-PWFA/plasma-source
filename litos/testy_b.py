#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:10:13 2017

@author: mike
"""

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import mike_math as mm
from calc_M import calc_M

if __name__ == '__main__':
    
    # define plasma bulk (flat-top) properties
    npl0   = 1e17 # cm^-3
    dEds0  = 6e9 # eV/m
    dgds0  = dEds0/nc.me
    L_ft   = 0.20 # m
    
    # define plasma up-ramp
    shape_up = 'xu5'
    hw_up    = 0.002 #0.147 #0.05 # m
    L_up     = 0.20 #1.00 # m
    top_up   = L_up # m
    
    # define plasma down-ramp
    shape_dn = shape_up
    hw_dn    = hw_up # m
    L_dn     = L_up # m
    top_dn   = 0  # m
    
    # define longitudinal steps
    ds   = 250.0e-6 # m
    s_ft = np.linspace(0,L_ft,round(L_ft/ds+1))
    s_up = np.linspace(0,L_up,round(L_up/ds+1))
    s_dn = np.linspace(0,L_dn,round(L_dn/ds+1))
    
    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
    plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp)
    
    # define beam parameters
    gbC    = 20000 # relativistic lorentz factor
    eps    = 5.0e-6  # m-rad, normalized emittance
    beta   = 0.10  # m
    
#    # match beam to plasma
#    wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
#    kp0    = wp0/nc.c # m^-1, plasma wave number
#    kb     = kp0/np.sqrt(2*gbC)
#    beta   = 1.0/kb
#
#    beta   = 100*beta
#
#    eps    = 1.0/(4*kb*gbC)

    alpha  = 0.00
    gamma  = (1.0+alpha**2)/beta # 1/m
    dgb    = 0.01
    dz     = 0
    npart  = 100
    dist   = 'gauss'
    
    # make beam
    s0     = 0.0
    twiss = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts = pb.make_parts(twiss,npart,dist)
    ebeam = pb.make_ebeam(s0,twiss[0],parts[0])

#    # add offset to bunch
#    ebeam0[0]["x"] += 10e-6 # m
    
    # set beam waist position
    waist = 0.00 #-0.430 #-0.122 # m, waist location w.r.t L_up
    s_w   = L_up + waist # m
    
    # propagate beam backward from waist to start of plasma
    pbp.prop_ebeam_drift(ebeam,[0,-s_w],last_only=True)
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam = pb.make_ebeam(s0,twiss[len(ebeam)-1],parts[len(ebeam)-1])
    vbeam = ebeam.copy()

    # propagate beam through plasma
    pbp.prop_ebeam_plasma(ebeam,plasma,last_only=False)

    # propagate beam through vacuum
    pbp.prop_ebeam_drift(vbeam,plasma["s"],last_only=False)