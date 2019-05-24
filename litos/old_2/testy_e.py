#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:07:38 2017

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
from matplotlib import cm
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import mike_math as mm
from calc_M import calc_M

import time as time

def scan_waist_lens(ebeam0,plasma0,waist,lens_srel,k):

    # set beam waist position
    s_w = plasma0["up_ramp"]["L"] + waist
    
    # make local copy of initial ebeam and plasma
    s0     = 0
    twiss0 = pb.get_twiss(ebeam0,len(ebeam0)-1)
    parts0 = pb.get_parts(ebeam0,len(ebeam0)-1)
    ebeam  = pb.make_ebeam(s0,twiss0,parts0)
    plasma = plasma0
    
    # propagate beam backward from waist to start of plasma
    ebeam = pbp.prop_ebeam_drift(ebeam,[0,-s_w],last_only=False)
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam = pb.make_ebeam(s0,twiss,parts)
    
    # insert plasma lens
    lens_npl0 = 1e17
    lens_L    = 100e-6
    lens_s0   = plasma0["up_ramp"]["L"] + lens_srel
    plasma    = ps.insert_lens(plasma0,lens_npl0,lens_L,lens_s0,add='no')

    # propagate beam through plasma
    ebeam = pbp.prop_ebeam_plasma(ebeam,plasma,last_only=False)

#    s     = np.zeros(len(ebeam))
#    beta  = np.zeros(len(ebeam))
#    for i in range(0,len(ebeam)):
#        s[i] = ebeam[i]["s"]
#        beta[i] = ebeam[i]["beta"]
#
#    ax2.scatter(s,beta)
#    ax1.scatter(s,plasma["npl"])
 

    # calculate mismatch parameter
    Tbeam  = [ebeam[len(ebeam)-1]["beta"],\
              ebeam[len(ebeam)-1]["alpha"],\
              ebeam[len(ebeam)-1]["gamma"]]
    wp0    = (5.64e4)*np.sqrt(plasma["npl"][-1]) # rad/s, plasma ang. freq.
    kp0    = wp0/nc.c # m^-1, plasma wave number
    kb     = kp0/np.sqrt(2*ebeam[len(ebeam)-1]["gbC"])
    Tmatch = [1.0/kb,0,kb]
    M      = calc_M(Tbeam,Tmatch)
    return M


if __name__ == '__main__':
    
    # define plasma bulk (flat-top) properties
    npl0   = 1e17 # cm^-3
    dEds0  = 6.00e9 # eV/m
    dgds0  = dEds0/nc.me
    L_ft   = 0.01 # m
    
    # define plasma up-ramp
    shape_up = 'gauss'
    hw_up    = 0.05 # m
    L_up     = 2.00 # m
    top_up   = L_up # m
    
    # define plasma down-ramp
    shape_dn = shape_up
    hw_dn    = hw_up # m
    L_dn     = 0 #L_up # m
    top_dn   = 0  # m
    
    # define longitudinal steps
    ds   = 0.0001 # m
    s_ft = np.linspace(0,L_ft,round(L_ft/ds+1))
    s_up = np.linspace(0,L_up,round(L_up/ds+1))
    s_dn = np.linspace(0,L_dn,round(L_dn/ds+1))
    
    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
    plasma0 = ps.make_plasma(bulk,up_ramp,dn_ramp)
    
    # define beam parameters
    gbC    = 20000 # relativistic lorentz factor
    eps    = 5e-6  # m-rad, normalized emittance
    beta   = 0.10 # m
    alpha  = 0.00
    gamma  = (1.0+alpha**2)/beta # 1/m
    dgb    = 0.01
    dz     = 0
    npart  = 0
    dist   = 'gauss'
    
    # make beam
    s0     = 0.0
    twiss0 = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts0 = pb.make_parts(twiss0,npart,dist)
    ebeam0 = pb.make_ebeam(s0,twiss0,parts0)
    
    # specify waist scan values
    nwaist = 51
    waist  = np.linspace(-0.20,-0.06,nwaist) # m, waist location w.r.t. L_up

    # specify relative position of thin plasma lens
    nlens_srel = 51
    lens_srel  = np.linspace(-0.22,-0.12,nlens_srel)

#    fig = plt.figure()
#    ax1 = fig.add_subplot(111)
#    
#    fig = plt.figure()
#    ax2 = fig.add_subplot(111)

    # perform scan
    num_cores = multiprocessing.cpu_count()
#    num_cores = 1
    M = Parallel(n_jobs=num_cores)\
        (delayed(scan_waist_lens)(ebeam0,plasma0,\
         waist[int(k/nlens_srel)],lens_srel[k%nlens_srel],k)\
         for k in range(nwaist*nlens_srel))

    M = np.reshape(M,[nwaist,nlens_srel])
#    print(M)
    
    # analyze results
    
    # plot results
    
    X = np.tile(waist.reshape(-1,1),(1,nlens_srel))
    Y = np.tile(lens_srel.T,(nwaist,1))
    
    
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,np.log10(M),100,\
                cmap=cm.Vega20c,\
                linewidth=2.0)
#    plt.scatter(-0.43,0.147,color='k')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$log_{10}$(M)')
    plt.ylabel(r'lens position [m]')
    plt.xlabel(r'waist position [m]')
    plt.title(r'beam matching optimization')
    
    
    
    
    
    
    
    
    levels = np.array([1.0,1.1,1.2,1.3,1.4,1.5,\
                       2.0,3.0,4.0,5.0,])
    labels = np.array([1.1,1.5,2.0,3.0,4.0,5.0])
    
    fig, axes = plt.subplots(1,1, sharey=True)
    
#    plt.contourf(X,Y,np.log10(M),100,\
#                cmap=cm.Vega20c,\
#                linewidth=2.0)
#    im = plt.imshow(M,cmap=cm.gray,origin='lower',extent=(-0.25,0,-0.25,0))

#    plt.contourf(X,Y,M,levels,cmap=cm.tab20b)

    CS = plt.contour(X,Y,M,levels,cmap=cm.tab20b)
    plt.clabel(CS,labels,fontsize=9, inline=1,fmt='%1.1f')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'M')
    cbar.set_ticks(levels)
#    cbar.set_ticklabels(levels)
    plt.ylabel(r'lens position [m]')
    plt.xlabel(r'waist position [m]')
    plt.title(r'beam matching optimization')
    
    
#    
#    min_y = np.zeros(M.shape[0])
#    min_x = np.zeros(M.shape[0])
#    for i in range(0,M.shape[0]):
#        min_y[i] = Y[0,i]
#        min_x[i] = X[np.argmin(M[:,i]),0]
#    
#    plt.scatter(min_x,min_y, s=10, c='b', marker=".", label='first')
#    
#    plt.xlabel(r'$waist$')
#    plt.ylabel(r'ramp half-length (m)')
#    plt.title("Ramp Type: %s" % shape_up)
#    plt.show()
