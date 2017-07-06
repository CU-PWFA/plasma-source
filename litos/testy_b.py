#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:10:13 2017

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp

if __name__ == '__main__':
    
    # define plasma bulk (flat-top) properties
    npl0   = 1e17 # cm^-3
    dEds0  = 20.00e9 # eV/m
    dgds0  = dEds0/nc.me
    L_ft   = 0.50 #0.50 # m
    
    # define plasma up-ramp
    shape_up = 'gauss'
    hw_up    = 0.15 # m
    L_up     = 0.30 # m
    top_up   = L_up # m
    
    # define plasma down-ramp
    shape_dn = shape_up
    hw_dn    = hw_up # m
    L_dn     = L_up # m
    top_dn   = 0  # m
    
    # define longitudinal steps
    ds   = 0.01 # m
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
    eps    = 5e-6  # m-rad, normalized emittance
    beta   = 0.10 # m
    alpha  = 0.00
    gamma  = (1.0+alpha**2)/beta # 1/m
    dgb    = 0.01
    dz     = 0
    npart = 10000
    dist  = 'gauss'
    
    # make beam
    s0     = 0.0
    twiss0 = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts0 = pb.make_parts(twiss0,npart,dist)
    ebeam  = pb.make_ebeam(s0,twiss0,parts0)
    
    # set beam waist position
    s_w   = L_up -0.067 #-0.350 #-0.963 # m
    
    # propagate beam backward from waist to start of plasma
    ebeam = pbp.prop_ebeam_drift(ebeam,[0,-s_w],last_only=True)

    # propagate beam through plasma
    ebeam = pbp.prop_ebeam_plasma(ebeam,plasma,last_only=False)

    
    # analyze results
    
    
    # plot results
    
    
    
    #x      = np.zeros(len(parts[0]))
    #xp     = np.zeros(len(parts[0]))
    #y      = np.zeros(len(parts[0]))
    #yp     = np.zeros(len(parts[0]))
    #z      = np.zeros(len(parts[0]))
    #gb     = np.zeros(len(parts[0]))
    #
    #for j in range(0,len(parts[i])):
    #        [x[j],xp[j],y[j],yp[j],z[j],gb[j]] = parts[len(parts)-1,j,:]
    fig = plt.figure()
    plt.hist(ebeam[0]["x"],25)
    fig = plt.figure()
    plt.hist(ebeam[0]["y"],25)
    
    s    = np.zeros(len(ebeam))
    beta = np.zeros(len(ebeam))
    for i in range(0,len(ebeam)-1):
        s[i] = ebeam[i]["s"]
        beta[i] = ebeam[i]["beta"]
    
    fig = plt.figure()
    plt.scatter(s,beta)

