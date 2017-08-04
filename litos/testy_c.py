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

def scan_waist_hw_up(ebeam0,plasma0,waist,hw_up,k):

    # set beam waist position
    s_w = plasma0["up_ramp"]["L"].copy() + waist
    
    # make local copy of initial ebeam and plasma
#    s0     = 0
#    twiss0 = pb.get_twiss(ebeam0,len(ebeam0)-1)
#    parts0 = pb.get_parts(ebeam0,len(ebeam0)-1)
#    ebeam  = pb.make_ebeam(s0,twiss0,parts0)
#    plasma = plasma0
    
    ebeam  = ebeam0.copy()
    plasma = plasma0.copy()
    
    # propagate beam backward from waist to start of plasma
#    ebeam = pbp.prop_ebeam_drift(ebeam,[0,-s_w],last_only=False)
#    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
#    parts = pb.get_parts(ebeam,len(ebeam)-1)
#    ebeam = pb.make_ebeam(s0,twiss,parts)
    
    pbp.prop_ebeam_drift(ebeam,[0,-s_w],last_only=False)
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam = pb.make_ebeam(s0,twiss[len(ebeam)-1],parts[len(ebeam)-1])
    
    
    # modify up-ramp
    up_ramp  = plasma["up_ramp"]
    s_up     = up_ramp["s"]
    shape_up = up_ramp["shape"]
    top_up   = up_ramp["top_loc"]
    npl0     = up_ramp["npl0"]
    dgds0    = up_ramp["dgds0"]
    up_ramp  = ps.make_ramp(s_up,"up",shape_up,hw_up,top_up,npl0,dgds0)

    # make plasma
    bulk     = plasma["bulk"]
    dn_ramp  = plasma["dn_ramp"]
    plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp)
    
    ps.deform_plasma(plasma,'sin',0.01,0.05)

    # propagate beam through plasma
    pbp.prop_ebeam_plasma(ebeam,plasma,last_only=False)

    # calculate mismatch parameter
    eps    = ebeam[len(ebeam)-1]["eps"]
    beta   = ebeam[len(ebeam)-1]["beta"]
    alpha  = ebeam[len(ebeam)-1]["alpha"]
    gamma  = ebeam[len(ebeam)-1]["gamma"]
    wp0    = (5.64e4)*np.sqrt(plasma["npl"][-1]) # rad/s, plasma ang. freq.
    kp0    = wp0/nc.c # m^-1, plasma wave number
    kb     = kp0/np.sqrt(2*gbC)
    beta_m = 1.0/kb
    Tbeam  = [beta,alpha,gamma]
    Tmatch = [beta_m,0,1.0/beta_m]
    R_x    = np.sqrt(eps*beta/gbC+(beta_m**2)*eps*gamma/gbC)
    
    M      = calc_M(Tbeam,Tmatch)
    K      = gbC*kb*R_x
    
    return [M,K]


if __name__ == '__main__':
    
    # define plasma bulk (flat-top) properties
    npl0   = 1e17 # cm^-3
    dEds0  = 6e9 # eV/m
    dgds0  = dEds0/nc.me
    L_ft   = 0.00 # m
    
    # define plasma up-ramp
    shape_up = 'gauss'
    hw_up    = 0.01 # m
    L_up     = 2.00 # m
    top_up   = L_up # m
    
    # define plasma down-ramp
    shape_dn = shape_up
    hw_dn    = hw_up # m
    L_dn     = 0 #L_up # m
    top_dn   = 0  # m
    
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
    twiss  = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts  = pb.make_parts(twiss[0],npart,dist)
    ebeam0 = pb.make_ebeam(s0,twiss[0],parts[0])
    
    # define longitudinal steps
    wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
    kp0    = wp0/nc.c # m^-1, plasma wave number
    kb     = kp0/np.sqrt(2*gbC)
    
    ds   = (np.pi/kb)*(1./10.) # m
    s_ft = np.linspace(0,L_ft,int(L_ft/ds+1))
    s_up = np.linspace(0,L_up,int(L_up/ds+1))
    s_dn = np.linspace(0,L_dn,int(L_dn/ds+1))
    
    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
    plasma0 = ps.make_plasma(bulk,up_ramp,dn_ramp)
    
    # specify waist scan values
    nwaist = 101
    waist  = np.linspace(-0.55,-0.300,nwaist) # m, waist location w.r.t. L_up
#    waist  = np.linspace(-0.35,-0.55,nwaist) # m, waist location w.r.t. L_up
    # specify ramp half-width scan values
    nhw_up = 101
    hw_up  = np.linspace(0.18,0.12,nhw_up) # m, HWHM of up-ramp
#    hw_up  = np.linspace(0.12,0.18,nhw_up) # m, HWHM of up-ramp

    # perform scan
    num_cores = multiprocessing.cpu_count()
    MK = Parallel(n_jobs=num_cores)\
            (delayed(scan_waist_hw_up)(ebeam0,plasma0,\
             waist[int(k/nhw_up)],hw_up[k%nhw_up],k)\
             for k in range(nwaist*nhw_up))

    M = np.zeros(len(MK))
    K = np.zeros(len(MK)) 
    for i in range(len(MK)):
        M[i] = MK[i][0]
        K[i] = MK[i][1]
    
    M = np.reshape(M,[nwaist,nhw_up])
#    print(M)

    K = np.reshape(K,[nwaist,nhw_up])
#    print(K)

    # analyze results
    
    # find location of min(M)
    i_Mmin_x = np.argmin(np.min(M,1))
    i_Mmin_y = np.argmin(np.min(M,0))
    Mmin_x = waist[i_Mmin_x]
    Mmin_y = hw_up[i_Mmin_y]
    
    print('min waist: ',Mmin_x)
    print('min ramp: ',Mmin_y)
    
    # plot results
    
    # filled color contour map of M
    X = np.tile(waist.reshape(-1,1),(1,nhw_up))
    Y = np.tile(hw_up.T,(nwaist,1))
    
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,np.log10(M),100,\
                cmap=cm.Vega20c,\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Mmin_x,Mmin_y,color='k')
    cbar.ax.set_ylabel(r'$log_{10}$(M)')
    plt.ylabel(r'ramp half-width [m]')
    plt.xlabel(r'waist position [m]')
    plt.title(r'beam matching for %s ramp'%shape_up)
    
#    # filled color contour map of log10(K)
#    fig, axes = plt.subplots(1,1, sharey=True)
#    plt.contourf(X,Y,K,100,\
#                cmap=cm.Vega20c,\
#                linewidth=2.0)
##    plt.scatter(-0.43,0.147,color='k')
#    cbar = plt.colorbar()
#    cbar.ax.set_ylabel(r'K')
#    plt.ylabel(r'ramp width [m]')
#    plt.xlabel(r'waist position [m]')
#    plt.title(r'beam matching optimization')
    
    # thin line contour map of M
    levels = np.array([1.0,1.1,1.2,1.3,1.4,1.5,\
                       2.0,3.0,4.0,5.0,])
    labels = np.array([1.1,1.5,2.0,3.0,4.0,5.0])
    
    fig, axes = plt.subplots(1,1, sharey=True)
    CS = plt.contour(X,Y,M,levels,cmap=cm.tab20b)
    plt.clabel(CS,labels,fontsize=9, inline=1,fmt='%1.1f')
    cbar = plt.colorbar()
    plt.scatter(Mmin_x,Mmin_y,color='k')
    cbar.ax.set_ylabel(r'M')
    cbar.set_ticks(levels)
#    cbar.set_ticklabels(levels)
    plt.ylabel(r'ramp half-width [m]')
    plt.xlabel(r'waist position [m]')
    plt.title(r'beam matching for %s ramp'%shape_up)
    
    
#    # thin line contour map of K
##    levels = np.array([1.0,2.0,3.0,4.0,5.0,\
##                       6.0,7.0,8.0,9.0,\
##                       10,15,20,25])
##    labels = np.array([1.0,2.0,3.0,4.0,5.0,\
##                       6.0,7.0,8.0,9.0,\
##                       10,15,20,25])
#    
#    levels = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
#    labels = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
#    
#    fig, axes = plt.subplots(1,1, sharey=True)
#    CS = plt.contour(X,Y,K,levels,cmap=cm.tab20b)
#    plt.clabel(CS,labels,fontsize=9, inline=1,fmt='%1.1f')
#    cbar = plt.colorbar()
#    cbar.ax.set_ylabel(r'K')
#    cbar.set_ticks(levels)
##    cbar.set_ticklabels(levels)
#    plt.ylabel(r'ramp width [m]')
#    plt.xlabel(r'waist position [m]')
#    plt.title(r'beam matching for gauss ramp')
    
    
    
    
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
