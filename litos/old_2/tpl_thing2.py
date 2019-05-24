#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 17:20:48 2018

@author: mike
"""

import numpy as np
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp

if __name__ == '__main__':
    
    # define beam parameters
    npart  = 100   # number of macro particles
    dist   = 'gauss' # distribution shape in trace space
    gbC    = (10e9)/nc.me   # centroid relativistic lorentz factor
    dgb    = 0.001    # relative energy spread (HWHM)
    dz     = 0       # spread in z (HWHM)
    eps    = 5.0e-6  # m-rad, normalized emittance
    beta   = 0.10    # m, Twiss at vac. waist
    alpha  = 0.00    # Twiss at vac. waist
    gamma  = (1.0+alpha**2)/beta # 1/m, Twiss at vac. waist
    auto_match = False # auto-match beam to plasma flat-top
    
    # define plasma up-ramp
    shape_up = 'gauss' # shape of ramp
    hw_up    = 0.1403 # m, half-width of ramp
    L_up     = 0.00 #0.75 # m, full length of ramp
    top_up   = L_up    # m, relative location of ramp top
    
    # define plasma down-ramp
    shape_dn = shape_up # shape of ramp
    hw_dn    = 0.00 #hw_up    # m, half-width of ramp
    L_dn     = 0.10 #L_up     # m, full length of ramp
    top_dn   = 0        # m, relative location of ramp top
    
    # define thin plasma lens properties

    # choose either the focal length or the magnification 
    sigx   = np.sqrt(eps*beta/gbC)
    f      = gbC*sigx/5 #beta/2 # m, must be < beta/2 at lens 
    Kl     = (1+np.sqrt(1-4*((f/beta)**2)))/(2*f) # m 
    M      = f/(Kl*(beta**2)) 
    
#    M      = 0.5 # = waist beta / beta at lens 
#    Kl     = np.sqrt((1/M)-1)/beta 
#    f      = M*Kl*(beta**2) 
    
    print('M: ',M)
    print('f: ',f)

    npl0   = 1e17 # cm^-3 
    lp     = (5.65e13)*Kl*gbC/(npl0*(1e6)) # m 
    dEds0  = 0 # np.sqrt(npl0/(5e16))*16.67e9  # eV/m, energy gain rate
    dgds0  = dEds0/nc.me               # 1/m, energy gain rate for rel. gamma
    L_ft   = lp
    
    # calculate betatron wave number in flat-top
    wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq. (flat-top)
    kp0    = wp0/nc.c               # 1/m, plasma wave number (flat-top)
    kb0    = kp0/np.sqrt(2*gbC)     # 1/m, betatron wave number (flat-top)
    
    # auto-match beam to plasma flat-top
    if (auto_match):
        beta   = 1.0/kb0
        alpha  = 0.00    # Twiss at vac. waist
        gamma  = (1.0+alpha**2)/beta # 1/m, Twiss at vac. waist

    # make beam at vac. waist
    s0     = 0.0 # m
    twiss  = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts  = pb.make_parts(twiss[0],npart,dist)
    ebeam  = pb.make_ebeam(s0,twiss[0],parts[0])
    
    
    # set beam offset
    sigx = np.sqrt(eps*beta/gbC)
    x_offset = 2.5*sigx #1./kp0-3*sigx
    print('sigx: ',sigx)
    print('kp-1: ',1/kp0)
    print('x offset: ',x_offset)
    pb.add_offset(ebeam,'x',x_offset)
#    pb.add_offset(ebeam,'xp',-1e-4)
    
    
    # set beam waist position
    waist = 0.00 #(-0.3884)      # m, waist location w.r.t L_up
    s_w   = L_up + waist # m, absolute wasit location
    
    # define longitudinal steps
    ds   = (1.0/kb0)*(1./10.)                 # m, step size
    s_ft = np.linspace(0,L_ft,int(L_ft/ds+1)) # m, steps for flat-top
#    ds   = beta/10.
    s_up = np.linspace(0,L_up,int(L_up/ds+1)) # m, steps for up-ramp
    s_dn = np.linspace(0,L_dn,int(L_dn/ds+1)) # m, steps for down-ramp
       
    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
    plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp) # output: plasma dict.

    # propagate beam backward from vac. waist to start of simulation
    pbp.prop_ebeam_drift(ebeam,[0,-s_w])
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam = pb.make_ebeam(s0,twiss[len(ebeam)-1],parts[len(ebeam)-1])
    vbeam = ebeam.copy()

    # propagate beam through plasma
    pbp.prop_ebeam_plasma(ebeam,plasma) # output: ebeam dict.

    # propagate beam through vacuum
    pbp.prop_ebeam_drift(vbeam,plasma["s"]) # output: vbeam dict.