#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 20:46:26 2017

@author: litos
"""


import numpy as np
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    # define plasma bulk (flat-top) properties
    npl0   = 10*5e16                      # cm^-3, plasma density
    dEds0  = np.sqrt(npl0/(5e16))*16.67e9 # eV/m, energy gain rate
    dgds0  = dEds0/nc.me               # 1/m, energy gain rate for rel. gamma
    L_ft   = 0.00                      # m, length of flat-top
    
    # define plasma up-ramp
    shape_up = 'adiabatic' # shape of ramp
    hw_up    = 0.01  # m, half-width of ramp
    L_up     = 0.00     # m, full length of ramp
    top_up   = L_up    # m, relative location of ramp top
    
    # define plasma down-ramp
    shape_dn = 'adiabatic' # shape of ramp
    hw_dn    = 0.01    # m, half-width of ramp
    L_dn     = 10.00    # m, full length of ramp
    top_dn   = 0       # m, relative location of ramp top
    
    # define beam parameters
    npart  = 0    # number of macro particles
    dist   = 'gauss' # distribution shape in trace space
    gbC    = (100e9)/nc.me # centroid relativistic lorentz factor
    dgb    = 0.001    # relative energy spread (HWHM)
    dz     = 0       # spread in z (HWHM)
    eps    = 5.0e-6  # m-rad, normalized emittance
    beta   = 0.10    # m, Twiss at vac. waist
    alpha  = 0.00    # Twiss at vac. waist
    gamma  = (1.0+alpha**2)/beta # 1/m, Twiss at vac. waist
    auto_match = True # auto-match beam to plasma flat-top
    
    # calculate betatron wave number in flat-top
    wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq. (flat-top)
    kp0    = wp0/nc.c               # 1/m, plasma wave number (flat-top)
    kb     = kp0/np.sqrt(2*gbC)     # 1/m, betatron wave number (flat-top)
    
    # auto-match beam to plasma flat-top
    if (auto_match):
        beta   = 1.0/kb
        alpha  = 0.00    # Twiss at vac. waist
        gamma  = (1.0+alpha**2)/beta # 1/m, Twiss at vac. waist

    # make beam at vac. waist
    s0     = 0.0 # m
    twiss  = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts  = pb.make_parts(twiss[0],npart,dist)
    ebeam  = pb.make_ebeam(s0,twiss[0],parts[0])
    
#    # set beam waist position
#    waist = 0.00        # m, waist location w.r.t L_up
#    s_w   = L_up + waist # m, absolute wasit location
    
    # define longitudinal steps
    ds   = (1.0/kb)*(1./10.)                 # m, step size
    s_ft = np.linspace(0,L_ft,int(L_ft/ds+1)) # m, steps for flat-top
    s_up = np.linspace(0,L_up,int(L_up/ds+1)) # m, steps for up-ramp
    s_dn = np.linspace(0,L_dn,int(L_dn/ds+1)) # m, steps for down-ramp
       
    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0,gbC)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0,gbC)
    plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp) # output: plasma dict.

#    # propagate beam backward from vac. waist to start of simulation
#    pbp.prop_ebeam_drift(ebeam,[0,-s_w])
#    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
#    parts = pb.get_parts(ebeam,len(ebeam)-1)
#    ebeam = pb.make_ebeam(s0,twiss[len(ebeam)-1],parts[len(ebeam)-1])
    vbeam = ebeam.copy()

    # propagate beam through plasma
    pbp.prop_ebeam_plasma(ebeam,plasma) # output: ebeam dict.

    # propagate beam through vacuum
    pbp.prop_ebeam_drift(vbeam,plasma["s"]) # output: vbeam dict.
    
    #%%
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.semilogy(plasma["s"],plasma["npl"])
    plt.xlabel('z (m)')
    plt.ylabel(r'log10($n_p$)')
    
    #%%
    chop = 0
    plot_s = np.zeros(len(ebeam)-chop)
    plot_g = np.zeros(len(ebeam)-chop)
    plot_b = np.zeros(len(ebeam)-chop)
    plot_n = np.zeros(len(ebeam)-chop)
    for i in range(0,len(ebeam)-chop):
        
        plot_s[i] = ebeam[i]["s"]
        plot_g[i] = ebeam[i]["gbC"]*nc.me*1e-9
        plot_b[i] = ebeam[i]["beta"]
        plot_n[i] = plasma["npl"][i]
    
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.plot(plot_s,plot_b)
    plt.xlabel('z (m)')
    plt.ylabel(r'$\beta$')
    
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.plot(plot_s,plot_g)
    plt.xlabel('z (m)')
    plt.ylabel(r'$\gamma$')
    
    #%%
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.plot(plot_s,plot_n)
    plt.xlabel('z (m)')
    plt.ylabel(r'$n_p$')

    #%%
#    s_plot = np.zeros(np.floor(len(ebeam)/100))
#    beta_plot = np.zeros(np.floor(len(ebeam)/100))
#    n_plot = 0
#    for i in range(0,int(np.floor(len(ebeam)/100))):
#        j = int(i*100)
#        s_plot[n_plot] = ebeam[j]["s"]
#        beta_plot[n_plot] = ebeam[j]["beta"]
#        n_plot += 1
#    
#    fig, axes = plt.subplots(1,1, sharey=True)
#    plt.plot(s_plot,beta_plot)
#    plt.xlabel('z (m)')
#    plt.ylabel(r'$\beta$ (m)')

#    