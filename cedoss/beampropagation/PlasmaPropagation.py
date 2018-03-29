#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:30:12 2017

Functions for beam propagation using Mike's code

@author: chris
"""

import sys
import matplotlib.pyplot as plt
import scipy.stats as stats

sys.path.insert(0, "../../litos")

import numpy as np
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import beam_ana as ba

def_beta = 0.10
def_hwup = 0.14
def_waist = -0.387
def_ramp = 'gauss'
def_npl0 = 5e16
def_gbC = 19569.5
def_L_up = 1.5
def_dEds0 = 1

option_energyscale = 16.67e9 #To match Robert's measured 16.67 GeV/m energy gain in 5e16

def ReturnDefaultParams(beta_change=def_beta, hwup_change=def_hwup, waist_change=def_waist, 
                        ramp_change=def_ramp, npl0_change=def_npl0, gbC_change = def_gbC,
                        L_up_change=def_L_up, dEds0_change=def_dEds0):
    npl0   = npl0_change                      # cm^-3, plasma density
    dEds0  = np.sqrt(npl0/(5e16))*option_energyscale * dEds0_change # eV/m, energy gain rate
    shape_up = ramp_change # shape of ramp
    hw_up    = hwup_change  # m, half-width of ramp
    L_up     = L_up_change     # m, full length of ramp, 1.5 for the Gaussian ramps
    alpha  = 0.00    # Twiss at vac. waist
    beta   = beta_change    # m, Twiss at vac. waist
    gbC    = gbC_change   # centroid relativistic lorentz factor
    wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq. (flat-top)
    kp0    = wp0/nc.c               # 1/m, plasma wave number (flat-top)
    waist  = waist_change        # m, waist location w.r.t L_up
    def_params = {
            # define plasma bulk (flat-top) properties
            'npl0'   : npl0,                      # cm^-3, plasma density
            'dEds0'  : dEds0, # eV/m, energy gain rate
            'dgds0'  : dEds0/nc.me,               # 1/m, energy gain rate for rel. gamma
            'L_ft'   : 0.50,                      # m, length of flat-top
            
            # define plasma up-ramp
            'shape_up' : shape_up, # shape of ramp
            'hw_up'    : hw_up,  # m, half-width of ramp
            'L_up'     : L_up,     # m, full length of ramp
            'top_up'   : L_up,    # m, relative location of ramp top
            
            # define plasma down-ramp
            'shape_dn' : shape_up, # shape of ramp
            'hw_dn'    : hw_up,    # m, half-width of ramp
            'L_dn'     : L_up,     # m, full length of ramp
            'top_dn'   : 0,        # m, relative location of ramp top
            
            # define beam parameters
            'npart'  : 1000,    # number of macro particles
            'dist'   : 'gauss', # distribution shape in trace space
            'gbC'    : gbC,   # centroid relativistic lorentz factor
            'dgb'    : 0.01,    # relative energy spread (HWHM)
            'dz'     : 0,       # spread in z (HWHM)
            'eps'    : 7.0e-6,  # m-rad, normalized emittance
            'beta'   : beta,    # m, Twiss at vac. waist
            'alpha'  : alpha,    # Twiss at vac. waist
            'gamma'  : (1.0+alpha**2)/beta, # 1/m, Twiss at vac. waist
            'auto_match' : False, # auto-match beam to plasma flat-top
            
            # calculate betatron wave number in flat-top
            'wp0'    : wp0, # rad/s, plasma ang. freq. (flat-top)
            'kp0'    : kp0,               # 1/m, plasma wave number (flat-top)
            'kb'     : kp0/np.sqrt(2*gbC),     # 1/m, betatron wave number (flat-top)
            
            # set beam waist position
            'waist' : waist,        # m, waist location w.r.t L_up
            's_w'   : L_up + waist, # m, absolute wasit location
            's0'    : 0 #m
            }
    return def_params

"""
# auto-match beam to plasma flat-top
if (auto_match):
    beta   = 1.0/kb
    alpha  = 0.00    # Twiss at vac. waist
    gamma  = (1.0+alpha**2)/beta # 1/m, Twiss at vac. waist
"""

# make beam at vac. waist
def CallMakeTwiss(params):
    twiss = pb.make_twiss(params['beta'],params['alpha'],
                          params['gamma'],params['eps'],
                          params['gbC'],params['dgb'],params['dz'])
    return twiss

def CallMakeParts(twiss, params):
    parts  = pb.make_parts(twiss[0],params['npart'],params['dist'])
    return parts

def CallMakeBeam(twiss, parts, params):
    ebeam  = pb.make_ebeam(params['s0'],twiss[0],parts[0])
    return ebeam

def MakeBulkPlasma(params):
    # define longitudinal steps
    ds   = (1.0/params['kb'])*(1./10.)#10.)                  # m, step size
    s_ft = np.linspace(0,params['L_ft'],int(params['L_ft']/ds+1)) # m, steps for flat-top
    s_up = np.linspace(0,params['L_up'],int(params['L_up']/ds+1)) # m, steps for up-ramp
    s_dn = np.linspace(0,params['L_dn'],int(params['L_dn']/ds+1)) # m, steps for down-ramp
    
    npl0 = params['npl0']; dgds0 = params['dgds0']; shape_up=params['shape_up']
    hw_up = params['hw_up']; top_up = params['top_up']; hw_dn = params['hw_dn']
    shape_dn = params['shape_dn']; top_dn = params['top_dn']
    
    # make plasma
    bulk    = ps.make_bulk(s_ft,npl0,dgds0)
    up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
    dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
    plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp) # output: plasma dict.
    return plasma

def InsertPlasmaLens(tpl_n, tpl_l, tpl_offset, plasma):
    tpl_s0   = plasma["up_ramp"]["L"] + tpl_offset
    ps.insert_lens(plasma,tpl_n,tpl_l,tpl_s0,add='')
    return plasma

def PropagateBackwards(ebeam, params):
    # propagate beam backward from vac. waist to start of simulation
    pbp.prop_ebeam_drift(ebeam,[0,-1*params['s_w']])
    twiss = pb.get_twiss(ebeam,len(ebeam)-1)
    parts = pb.get_parts(ebeam,len(ebeam)-1)
    ebeam = pb.make_ebeam(params['s0'],twiss[len(ebeam)-1],parts[len(ebeam)-1])
    return ebeam

def PropagatePlasma(ebeam, plasma):
    # propagate beam through plasma
    pbp.prop_ebeam_plasma(ebeam,plasma) # output: ebeam dict.
    return ebeam
    
def PropagateVirtual(ebeam, plasma):
    vbeam = ebeam.copy()
    # propagate beam through vacuum
    pbp.prop_ebeam_drift(vbeam,plasma["s"]) # output: vbeam dict.
    return vbeam

def CalcBmag(ebeam, plasma):
    i_flat_start = np.argwhere(plasma["s"]>=plasma["up_ramp"]["top_loc"])[0][0]
    Tbeta   = ebeam[i_flat_start]["beta"]
    Talpha  = ebeam[i_flat_start]["alpha"]
    Tgamma  = ebeam[i_flat_start]["gamma"]
    TgbC    = ebeam[i_flat_start]["gbC"]
    #TgbC    = TgbC-0.734*plasma["bulk"]["dgds0"]*plasma["up_ramp"]["hw"]
    Twp0    = (5.64e4)*np.sqrt(plasma["bulk"]["npl0"]) # rad/s, plasma ang. freq.
    Tkp0    = Twp0/nc.c # m^-1, plasma wave number
    Tkb     = Tkp0/np.sqrt(2*TgbC)
    Tbeta_m = 1.0/Tkb
    TTbeam  = [Tbeta,Talpha,Tgamma]
    TTmatch = [Tbeta_m,0,1.0/Tbeta_m]
    
    BB      = ba.calc_Bmag(TTbeam,TTmatch)
    return BB

def PlotPropagation(ebeam, vbeam, plasma):
    nstep     = len(ebeam)
    s         = np.zeros(nstep)
    beta      = np.zeros(nstep)
    v_beta    = np.zeros(nstep)
    rms_x_eps = np.zeros(nstep)
    J_kurt    = np.zeros(nstep)
    frac      = 1.00
    
    for i in range(0,nstep):
        ebeam_rms = ba.calc_ebeam_rms(ebeam,i,frac)
        s[i]     = ebeam[i]["s"]
        beta[i] = ebeam_rms["x_beta"]/(1e-2)
        v_beta[i] = vbeam[i]["beta"]/(1e-2)
        rms_x_eps[i] = ebeam_rms["x_eps"]/(1e-6)
        
        [u,v] = ba.real2norm_coords(ebeam[i]["x"],ebeam[i]["xp"],\
                                ebeam_rms["x_beta"],ebeam_rms["x_alpha"])
        J = (u**2+v**2)/2
        J_kurt[i] = stats.kurtosis(J,0,False,True)
    
    figA, (ax1, ax3) = plt.subplots(2, sharex=True, sharey=False)
    
    ax1.plot(s,v_beta/10,color='b',linestyle='dashed')
    ax1.plot(s,beta,color='b',linestyle='solid')
    ax1.set_xlim([0.5,3.0])
    ax1.set_ylim([0,4.0])
    ax1.set_ylabel(r'$\beta$ [cm]',color='b')
    ax1.tick_params('y',colors='b')
    
    npl = plasma["npl"]/plasma["bulk"]["npl0"]
    
    ax2  = ax1.twinx()
    ax2.plot(s,npl,color='g',linestyle='solid')
    ax2.set_ylabel(r'$n_p/n_{p,0}$',color='g')
    ax2.tick_params('y',colors='g')
    ax2.set_ylim([0,1.4])
    ax2.text(0.50, 0.80, r'$n_{p,0} = %2.1e$'%plasma["bulk"]["npl0"],
            verticalalignment='center', horizontalalignment='center',
            transform=ax2.transAxes,
            color='green', fontsize=12)
    
    BB = CalcBmag(ebeam, plasma)
    print('Bmag: ', BB)
    ax3.plot(s,rms_x_eps/rms_x_eps[0],color='k',linestyle='-')
    ax3.plot(s,BB*np.ones(len(s)),color='k',linestyle='-.')
    ax3.set_ylabel(r'$\varepsilon_n/\varepsilon_{n,0}$',color='k')
    ax3.tick_params('y',colors='k')
    ax3.set_xlim([0.5,3.0])
    ax3.set_xlabel('z [m]')
    ax3.set_ylim([0.975,1.075])
    
    ax4 = ax3.twinx()
    ax4.plot(s,J_kurt/J_kurt[0],color='r',linestyle='--')
    ax4.set_ylabel(r'$K_J/K_{J,0}$',color='r')
    ax4.tick_params('y',colors='r')
    ax4.set_ylim([0.975,1.075])
    """
    #These act goofy
    xlabel_locs = [0.5,1.0,1.5,2.0,2.5,3.0]
    xlabels = [0,0.5,1.0,1.5,2.0,2.5]
    plt.xticks(xlabel_locs, xlabels)
    
    ax1.set_xlim([0.5,3.0])
    ax3.set_xlim([0.5,3.0])
    #end goofs
    """
    figA.tight_layout()
    figA.subplots_adjust(hspace=0)
    
    plt.show()

def Calc2DTolerance(arr, x, y, thresh):
    minlocx = np.argmin(np.min(arr,1))
    minlocy = np.argmin(np.min(arr,0))
    dx = x[1]-x[0]
    xarr = arr[:,minlocy]
    print(" in +/- x:",str(Calc1DTolerance(xarr, dx, minlocx, thresh)/2))
    
    dy = y[1]-y[0]
    yarr = arr[minlocx,:]
    print(" in +/- y:",str(Calc1DTolerance(yarr, dy, minlocy, thresh)/2))

def Calc1DToleranceRange(arr, dx, minloc, thresh):
    i = minloc
    flag = 2
    top = i; bot = i
    
    while i < (len(arr) - 1):
        top = i
        if arr[i] < thresh:
            i = i + 1
        else:
            flag = flag - 1
            break
    
    i = minloc
    while i > 0:
        bot = i
        if arr[i] < thresh:
            i = i - 1
        else:
            flag = flag - 1
            break
    
    top = top - 1 + (thresh-arr[top-1])/(arr[top]-arr[top-1])
    bot = bot + 1 - (thresh-arr[bot+1])/(arr[bot]-arr[bot+1])
    
    tol_range = (top-bot)*dx
    if flag > 0:
        print("Warning: domain edge reached for " + str(thresh))
    return [tol_range,bot,top]
    
def Calc1DTolerance(arr, dx, minloc, thresh):
    result = Calc1DToleranceRange(arr, dx, minloc, thresh)
    return result[0]
    
def PlotContour(contour, x_arr, y_arr, x_label, y_label, log = 0):
        # find location of min(B)
    i_Bmin_x = np.argmin(np.min(contour,1))
    i_Bmin_y = np.argmin(np.min(contour,0))
    Bmin_x = x_arr[i_Bmin_x]
    Bmin_y = y_arr[i_Bmin_y]
    Bmin   = np.min(np.min(contour))
    
    print('matching x value: ',Bmin_x)
    print('matching y value: ',Bmin_y)
    print('matched B-mag: ',Bmin)
    print('emittance growth (%): ',100*(Bmin-1)/Bmin)
    
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,np.log10(contour),100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$log_{10}(B_m)$')
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    if log == 1:
        plt.yscale('log'); plt.xscale('log')
#    plt.title(r'beam matching for %s ramp'%shape_up)

    # thin line contour map of B
    levels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
    labels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
    
    fig, axes = plt.subplots(1,1, sharey=True)
    CS = plt.contour(X,Y,contour,levels,cmap=plt.get_cmap('Vega20b'))
    plt.clabel(CS,labels,fontsize=9,inline=1,fmt='%1.2f')
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$B_m$')
    cbar.set_ticks(levels)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    if log == 1:
        plt.yscale('log'); plt.xscale('log')
    plt.show()
#    plt.title(r'beam matching for %s ramp'%shape_up)

    return [Bmin_x,Bmin_y]