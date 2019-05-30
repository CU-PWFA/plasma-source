#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 19:05:39 2017

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
import beam_ana as ba

if __name__ == '__main__':

    # define vacuum propagation distance
    L_vac    = 1.5-0.28
    
    # define beam parameters
    gbC    = 20000 # relativistic lorentz factor
    eps    = 5.0e-6  # m-rad, normalized emittance
    beta   = 0.10  # m

    alpha  = 0.00
    gamma  = (1.0+alpha**2)/beta # 1/m
    dgb    = 0.01
    dz     = 0
    npart  = 1000
    dist   = 'gauss'
    
    # make beam
    s0     = 0.0
    twiss = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts = pb.make_parts(twiss[0],npart,dist)
    ebeam = pb.make_ebeam(s0,twiss[0],parts[0])

    # add offset to bunch
    ebeam[0]["x"] += 10e-6 # m
    
    # set beam waist position
    waist = 0 #-0.28 #-0.105 # m, waist location w.r.t L_up
    s_w   = L_vac + waist # m
    
    # define longitudinal steps
    ds   = 0.01 # m
    s_vac = np.linspace(0,L_vac,int(L_vac/ds+1))

    # propagate beam backward from waist to start of vacuum
    pbp.prop_ebeam_drift(ebeam,[0,-s_w],last_only=True)
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam = pb.make_ebeam(s0,twiss[len(ebeam)-1],parts[len(ebeam)-1])
    
    # propagate beam through vacuum
    pbp.prop_ebeam_drift(ebeam,s_vac,last_only=False)
    
    # get beam centroid
    ebeam_cent = ba.calc_ebeam_cent(ebeam,len(ebeam)-1,frac=1.00)
    print(ebeam_cent)