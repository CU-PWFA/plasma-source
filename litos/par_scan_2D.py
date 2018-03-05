#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 10:52:37 2017

@author: litos
"""

import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import beam_ana as ba

def scan_waist_hw_up(ebeam0,plasma0,waist,hw_up,k):

    # set beam waist position
    s_w = plasma0["up_ramp"]["L"].copy() + waist
    
    # make local copy of initial ebeam and plasma
    ebeam  = ebeam0.copy()
    plasma = plasma0.copy()
    
    # propagate beam backward from waist to start of plasma
    pbp.prop_ebeam_drift(ebeam,[0,-s_w])
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

    # propagate beam through plasma
    pbp.prop_ebeam_plasma(ebeam,plasma)

    # calculate mismatch parameter
    eps    = ebeam[len(ebeam)-1]["eps"]
    beta   = ebeam[len(ebeam)-1]["beta"]
    alpha  = ebeam[len(ebeam)-1]["alpha"]
    gamma  = ebeam[len(ebeam)-1]["gamma"]
    gbC    = ebeam[len(ebeam)-1]["gbC"]
    wp0    = (5.64e4)*np.sqrt(plasma["npl"][-1]) # rad/s, plasma ang. freq.
    kp0    = wp0/nc.c # m^-1, plasma wave number
    kb     = kp0/np.sqrt(2*gbC)
    beta_m = 1.0/kb
    Tbeam  = [beta,alpha,gamma]
    Tmatch = [beta_m,0,1.0/beta_m]
    
    B      = ba.calc_Bmag(Tbeam,Tmatch)
    
    return [B]


if __name__ == '__main__':
    
    # define plasma bulk (flat-top) properties
    npl0   = 5e16                      # cm^-3, plasma density
    dEds0  = np.sqrt(npl0/(5e16))*16.67e9 # eV/m, energy gain rate
    dgds0  = dEds0/nc.me               # 1/m, energy gain rate for rel. gamma
    L_ft   = 0.00                      # m, length of flat-top
    
    # define plasma up-ramp
    shape_up = 'gauss' # shape of ramp
    hw_up    = 0.132  # m, half-width of ramp
    L_up     = 5*hw_up     # m, full length of ramp
    top_up   = L_up    # m, relative location of ramp top
    
    # define plasma down-ramp
    shape_dn = shape_up # shape of ramp
    hw_dn    = hw_up    # m, half-width of ramp
    L_dn     = 0.00     # m, full length of ramp
    top_dn   = 0        # m, relative location of ramp top
    
    # define beam parameters
    npart  = 0       # number of macro particles
    dist   = 'gauss' # distribution shape in trace space
    gbC    = (10e9)/nc.me # centroid relativistic lorentz factor
    dgb    = 0.01    # relative energy spread (HWHM)
    dz     = 0       # spread in z (HWHM)
    eps    = 5.0e-6  # m-rad, normalized emittance
    beta   = 0.20    # m, Twiss at vac. waist
    alpha  = 0.00    # Twiss at vac. waist
    gamma  = (1.0+alpha**2)/beta # 1/m, Twiss at vac. waist
    
    # calculate betatron wave number in flat-top
    wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq. (flat-top)
    kp0    = wp0/nc.c               # 1/m, plasma wave number (flat-top)
    kb     = kp0/np.sqrt(2*gbC)     # 1/m, betatron wave number (flat-top)

    # make beam at vac. waist
    s0     = 0.0 # m
    twiss  = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts  = pb.make_parts(twiss[0],npart,dist)
    ebeam0 = pb.make_ebeam(s0,twiss[0],parts[0])
    
    # define longitudinal steps
    ds   = (1.0/kb)*(1./40.)                  # m, step size
    s_ft = np.linspace(0,L_ft,int(L_ft/ds+1)) # m, steps for flat-top
    s_up = np.linspace(0,L_up,int(L_up/ds+1)) # m, steps for up-ramp
    s_dn = np.linspace(0,L_dn,int(L_dn/ds+1)) # m, steps for down-ramp

    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
    plasma0 = ps.make_plasma(bulk,up_ramp,dn_ramp) # output: plasma dict.
    
    # specify waist scan values
    nwaist = 200 # number of values
    waist  = np.linspace(-0.70,-0.20,nwaist) # m, waist location w.r.t. L_up

    # specify ramp half-width scan values
    nhw_up = 200 # number of values
    hw_up  = np.linspace(0.08,0.24,nhw_up) # m, HWHM of up-ramp

    # perform scan
    num_cores = multiprocessing.cpu_count() # get number of available cores
#    num_cores = 1 # use 1 core if number of particles > 0
    B  = Parallel(n_jobs=num_cores)\
            (delayed(scan_waist_hw_up)(ebeam0,plasma0,\
             waist[int(k/nhw_up)],hw_up[k%nhw_up],k)\
             for k in range(nwaist*nhw_up))

    # reshape scan output into 2D array
    B = np.reshape(B,[nwaist,nhw_up])

    #%%
    # analyze results
    
    # find location of min(B)
    i_Bmin_x = np.argmin(np.min(B,1))
    i_Bmin_y = np.argmin(np.min(B,0))
    Bmin_x = waist[i_Bmin_x]
    Bmin_y = hw_up[i_Bmin_y]
    Bmin   = np.min(np.min(B))
    
    print('matching waist pos. (m): ',Bmin_x)
    print('matching ramp HWHM (m): ',Bmin_y)
    print('matched B-mag: ',Bmin)
    print('emittance growth (%): ',100*(Bmin-1)/Bmin)
    
    #%%
    # plot results
    
    X = np.tile(waist.reshape(-1,1),(1,nhw_up))
    Y = np.tile(hw_up.T,(nwaist,1))
    
#    # filled color contour map of log10(B)
#    fig, axes = plt.subplots(1,1, sharey=True)
#    plt.contourf(X,Y,np.log10(B),100,\
#                cmap=plt.get_cmap('Vega20c'),\
#                linewidth=2.0)
#    cbar = plt.colorbar()
#    plt.scatter(Bmin_x,Bmin_y,color='k')
#    cbar.ax.set_ylabel(r'$log_{10}(B_m)$')
#    plt.ylabel(r'$\sigma_{\rm hw}$ [m]')
#    plt.xlabel(r'$z_{\beta^{*}}$ [m]')
##    plt.title(r'beam matching for %s ramp'%shape_up)
#
#    # thin line contour map of B
#    levels = np.array([1.01,1.05,1.1,1.5,2.0,3.0,4.0,5.0])
#    labels = np.array([1.01,1.05,1.1,1.5,2.0,3.0,4.0,5.0])
#    fig, axes = plt.subplots(1,1, sharey=True)
#    CS = plt.contour(X,Y,B,levels,cmap=plt.get_cmap('Vega20b'))
#    plt.clabel(CS,labels,fontsize=9,inline=1,fmt='%1.2f')
#    cbar = plt.colorbar()
#    plt.scatter(Bmin_x,Bmin_y,color='k')
#    cbar.ax.set_ylabel(r'$B_m$')
#    cbar.set_ticks(levels)
#    plt.ylabel(r'$\sigma_{\rm hw}$ [m]')
#    plt.xlabel(r'$z_{\beta^{*}}$ [m]')
##    plt.title(r'beam matching for %s ramp'%shape_up)

    # thin line contour map of B with flipped axes
    levels = np.array([1.01,1.05,1.1,1.5,2.0,3.0,4.0,5.0])
    labels = np.array([1.01,1.05,1.1,1.5,2.0,3.0,4.0,5.0])
    fig, axes = plt.subplots(1,1, sharey=True)
    CS = plt.contour(Y,X,B,levels,cmap=plt.get_cmap('Vega20b'))
    plt.clabel(CS,labels,fontsize=9,inline=1,fmt='%1.2f')
#    cbar = plt.colorbar()
    plt.scatter(Bmin_y,Bmin_x,color='k')
#    cbar.ax.set_ylabel(r'$B_m$')
#    cbar.set_ticks(levels)
    plt.xlabel(r'$\sigma_{\rm hw}$ [m]',fontsize=16)
    plt.ylabel(r'$z_v^*$ [m]',fontsize=16)
#    plt.ylabel(r'$z_{\beta^{*}}$ [m]',fontsize=16)
#    plt.title(r'beam matching for %s ramp'%shape_up)