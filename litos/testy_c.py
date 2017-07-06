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

def scan_waist_hw_up(ebeam0,plasma,waist,hw_up,ax1,ax2,k):
    # get scan indices
    nhw_up = len(hw_up)
    i = int(k/nhw_up)
    j = k%nhw_up
            
    print(i*nhw_up+j)

    # get scan variables
    iwaist = waist[i]
    jhw_up = hw_up[j]

    # set beam waist position
    s_w = plasma["up_ramp"]["L"] + iwaist
    
    print(s_w)
    
    # propagate beam backward from waist to start of plasma
    s0     = 0
    ebeam  = pbp.prop_ebeam_drift(ebeam0,[0,-s_w],last_only=True)
    twiss0 = pb.get_twiss(ebeam,len(ebeam)-1)
    parts0 = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam  = pb.make_ebeam(s0,twiss0,parts0)
    
    print(len(ebeam))
    
    
#    # modify up-ramp
#    up_ramp  = plasma["up_ramp"]
#    s_up     = up_ramp["s"]
#    shape_up = up_ramp["shape"]
#    top_up   = up_ramp["top_loc"]
#    npl0     = up_ramp["npl0"]
#    dgds0    = up_ramp["dgds0"]
#    up_ramp  = ps.make_ramp(s_up,"up",shape_up,jhw_up,top_up,npl0,dgds0)
#
#    # make plasma
#    bulk     = plasma["bulk"]
#    dn_ramp  = plasma["dn_ramp"]
#    plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp)

    # propagate beam through plasma
    ebeam = pbp.prop_ebeam_plasma(ebeam,plasma,last_only=False)


    s     = np.zeros(len(ebeam))
    beta  = np.zeros(len(ebeam))
    for i in range(0,len(ebeam)):
        s[i] = ebeam[i]["s"]
        beta[i] = ebeam[i]["beta"]

    ax2.scatter(s,beta)
    ax1.scatter(s,plasma["npl"])
 

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
    npl0   = 1e-17 #1e17 # cm^-3
    dEds0  = 6.00e9 # eV/m
    dgds0  = dEds0/nc.me
    L_ft   = 0.50 # m
    
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
    
    # define longitudinal steps
    ds   = 0.001 # m
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
    npart  = 0
    dist   = 'gauss'
    
    # make beam
    s0     = 0.0
    twiss0 = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
    parts0 = pb.make_parts(twiss0,npart,dist)
    ebeam0 = pb.make_ebeam(s0,twiss0,parts0)
    
    # specify waist scan values
    nwaist = 3 #61
    waist  = np.linspace(-1.2,0.0,nwaist) # m, waist location w.r.t. L_up
    # specify ramp half-width scan values
    nhw_up = 3 #50
    hw_up  = np.linspace(0.01,0.50,nhw_up) # m, HWHM of up-ramp
    
#    # initialize mismatch matrix
#    M = np.zeros([nwaist,nhw_up])
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    fig = plt.figure()
    ax2 = fig.add_subplot(111)


    # perform scan
    num_cores = multiprocessing.cpu_count()
    num_cores = 1
    M = Parallel(n_jobs=num_cores)\
        (delayed(scan_waist_hw_up)(ebeam0,plasma,waist,hw_up,ax1,ax2,k) for k in range(nwaist*nhw_up))

    M = np.reshape(M,[nwaist,nhw_up])
    print(M)
    
    # analyze results
    
    # plot results
    
    X = np.tile(waist.reshape(-1,1),(1,nhw_up))
    Y = np.tile(hw_up.T,(nwaist,1))
    
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,np.log10(M),100,\
                cmap=cm.Vega20c,\
                linewidth=2.0)
#    plt.pcolor(X,Y,M)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$log_{10}$(M)')
    
    min_y = np.zeros(M.shape[0])
    min_x = np.zeros(M.shape[0])
    for i in range(0,M.shape[0]):
        min_y[i] = Y[0,i]
        min_x[i] = X[np.argmin(M[:,i]),0]
    
    plt.scatter(min_x,min_y, s=10, c='b', marker=".", label='first')
    
    plt.xlabel(r'$waist$')
    plt.ylabel(r'ramp half-length (m)')
    plt.title("Ramp Type: %s" % shape_up)
    plt.show()
