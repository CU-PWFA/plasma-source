#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 20:29:47 2018

@author: litos
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import mike_math as mm
from   calc_M import calc_M
import scipy.stats as stats
import scipy.optimize as opt
from collections import defaultdict
import beam_ana as ba
import scipy.spatial as spatial
from matplotlib.ticker import FormatStrFormatter

# do analysis on beam/plasma
nstep = len(ebeam)
npart = ebeam[0]["npart"]
s     = np.zeros(nstep)
gbC   = np.zeros(nstep)
eps   = np.zeros(nstep)
beta  = np.zeros(nstep)
alpha = np.zeros(nstep)
gamma = np.zeros(nstep)
rms_x = np.zeros(nstep)
rms_x_eps = np.zeros(nstep)
cent_x  = np.zeros(nstep)
cent_xp = np.zeros(nstep)
un_kurt = np.zeros(nstep)
J_kurt  = np.zeros(nstep)

kb = np.zeros(nstep)
vbeta  = np.zeros(nstep)
valpha = np.zeros(nstep)
dndz = np.zeros(nstep)

ds = ebeam[1]["s"]-ebeam[0]["s"]

frac = 1.0

for i in range(0,nstep):
    s[i]     = ebeam[i]["s"] # m\
    gbC[i]   = ebeam[i]["gbC"]
#    eps[i]   = ebeam[i]["eps"] # mm-mrad
    beta[i]  = ebeam[i]["beta"] # m
    alpha[i] = ebeam[i]["alpha"]
    gamma[i] = ebeam[i]["gamma"] # 1/m
#    ebeam_rms = ba.calc_ebeam_rms(ebeam,i,frac)
#    ebeam_cent = ba.calc_ebeam_cent(ebeam,i,frac)
#    cent_x[i]  = ebeam_cent["x"]/(1e-6)
#    cent_xp[i] = ebeam_cent["xp"]/(1e-6)
#    rms_x[i]     = ebeam_rms["x"]/(1e-6)
#    rms_x_eps[i] = ebeam_rms["x_eps"]/(1e-6)
    
    npl = plasma["npl"][i]
    wp    = (5.64e4)*np.sqrt(npl) # rad/s, plasma ang. freq.
    kp    = wp/nc.c               # 1/m, plasma wave number
    kb[i] = kp/np.sqrt(2*gbC[i])     # 1/m, betatron wave number
    
    vbeta[i]  = vbeam[i]["beta"]
    valpha[i] = vbeam[i]["alpha"]
    
    if i==0 or npl==0:
        dndz[i] = 0
    else:
        dndz[i] = (plasma["npl"][i]-plasma["npl"][i-1])/ \
        ((ebeam[i]["s"]-ebeam[i-1]["s"]))
    
#    [u,v] = ba.real2norm_coords(ebeam[i]["x"],ebeam[i]["xp"],\
#                                ebeam_rms["x_beta"],ebeam_rms["x_alpha"])
#    J = (u**2+v**2)/2
#    phi = np.arctan2(v,u)
#    un = np.sqrt(2*J*ebeam[i]["gbC"])*np.cos(phi)/np.sqrt(rms_x_eps[i])
#    un_kurt[i] = stats.kurtosis(un,0,True)
#    J_kurt[i] = stats.kurtosis(J,0,False,True)
    
dndz[0] = dndz[1]    

# normalized distance
skb = (s[2]-s[1])*np.abs(kb)

#i_flat_start = np.argwhere(plasma["s"]>=plasma["up_ramp"]["top_loc"])[0][0]
#
#Tbeta   = ebeam[i_flat_start]["beta"]
#Talpha  = ebeam[i_flat_start]["alpha"]
#Tgamma  = ebeam[i_flat_start]["gamma"]
#TgbC    = ebeam[i_flat_start]["gbC"]
##TgbC    = TgbC-0.734*plasma["bulk"]["dgds0"]*plasma["up_ramp"]["hw"]
#Twp0    = (5.64e4)*np.sqrt(plasma["bulk"]["npl0"]) # rad/s, plasma ang. freq.
#Tkp0    = Twp0/nc.c # m^-1, plasma wave number
#Tkb     = Tkp0/np.sqrt(2*TgbC)
#Tbeta_m = 1.0/Tkb
#TTbeam  = [Tbeta,Talpha,Tgamma]
#TTmatch = [Tbeta_m,0,1.0/Tbeta_m]
#
#BB      = ba.calc_Bmag(TTbeam,TTmatch)
#print('Bmag: ',BB)
    
ent_z_v = s[0]+alpha[0]/gamma[0]
#ent_beta_v = (1/gamma[0])/(1e-2)
#exit_z_v = s[nstep-1]+alpha[nstep-1]/gamma[nstep-1]
#exit_beta_v = (1/gamma[nstep-1])/(1e-2)

#%% beta and rms_x evolution through plasma

# which subfigure?
#subfig = 'a'
subfig = 'b'

figA, (ax1, ax3) = plt.subplots(2, sharex=True, sharey=False)

ax1.plot(skb,beta,color='b',linestyle='solid')
ax1.plot(skb,vbeta,color='b',linestyle='-.')
ax1.set_ylabel(r'$\beta_x$',color='b',fontsize=16)
ax1.tick_params('y',colors='b')
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_ylim([0,0.15])
#ax1.set_ylim([0.1,2])
#ax1.set_xlim([0.5,2.0])

#ax1.plot([ent_z_v,ent_z_v],[0,1],color='b',linestyle='dashed')
#ax1.plot([exit_z_v,exit_z_v],[0,1.0],color='b',linestyle='dashed')

#ax1.yaxis.set_ticks(np.arange(0.0, 3.5, 0.5))

## subfigure label
#if subfig=='a':
#    ax1.text(0.20, 0.85, r'(a)',
#        verticalalignment='center', horizontalalignment='center',
#        transform=ax1.transAxes,
#        color='black', fontsize=16)
#elif subfig=='b':
#    ax1.text(0.20, 0.85, r'(b)',
#        verticalalignment='center', horizontalalignment='center',
#        transform=ax1.transAxes,
#        color='black', fontsize=16)

npl = plasma["npl"]/plasma["bulk"]["npl0"]
ax2  = ax1.twinx()
ax2.plot(s,npl,color='g',linestyle='solid')
#ax2.plot(skb,dndz/kb,color='g',linestyle='-.')
ax2.set_ylabel(r'$n_p/n_{p,0}$',color='g',fontsize=16)
ax2.tick_params('y',colors='g')
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax2.set_ylim([0,1.2])
#ax2.set_ylim([1e-6,20])
#ax2.set_xlim([0.5,2.0])

#ax2.yaxis.set_ticks(np.arange(0.0, 1.4, 0.2))

## plasma density text
#if subfig=='a':
#    ax2.text(0.50, 0.85, r'$n_{p,0} = 5\times 10^{16} \,{\rm cm^{-3}}$',
#        verticalalignment='center', horizontalalignment='center',
#        transform=ax2.transAxes,
#        color='green', fontsize=16)
#elif subfig=='b':
#    ax2.text(0.50, 0.85, r'$n_{p,0} = 1\times 10^{18} \,{\rm cm^{-3}}$',
#        verticalalignment='center', horizontalalignment='center',
#        transform=ax2.transAxes,
#        color='green', fontsize=16)
#else:
#    ax2.text(0.50, 0.85, \
#             r'$n_{p,0} = %2.1e \,{\rm cm^{-3}}$'%plasma["bulk"]["npl0"],
#        verticalalignment='center', horizontalalignment='center',
#        transform=ax2.transAxes,
#        color='green', fontsize=16)
    
ax3.plot(s,alpha,color='k',linestyle='-')
ax3.plot(s,valpha,color='k',linestyle='-.')
ax3.set_ylabel(r'$\alpha_x$',color='k',fontsize=16)
ax3.tick_params('y',colors='k')
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_ylim([-4.5,4.5])
#ax3.set_ylim([0.001,90])
#if subfig=='a':
#    ax3.set_ylim([0.9875,1.0125]) # matched limits
#    ax3.yaxis.set_ticks(np.arange(0.990, 1.010, 0.005))
#elif subfig=='b':
#    ax3.set_ylim([0.9,1.9]) # mismatched limits
#    ax3.yaxis.set_ticks(np.arange(1.0, 2.0, 0.2))
#else:
#    ax3.set_ylim([0.9,1.9])
#ax3.set_xlim([0.5,2.0])

#ax4 = ax3.twinx()
#ax4.plot(skb,1./kb,color='r',linestyle='--')
#ax4.set_ylabel(r'$k_{\beta}^{-1} (m)$',color='r',fontsize=16)
#ax4.tick_params('y',colors='r')
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax4.set_ylim([0,1.1])
#ax4.set_ylim([0.001,90])
#if subfig=='a':
#    ax4.set_ylim([0.9875,1.0125]) # matched limits
#    ax4.yaxis.set_ticks(np.arange(0.990, 1.010, 0.005))
#elif subfig=='b':
#    ax4.set_ylim([0.9,1.9]) # mismatched limits
#    ax4.yaxis.set_ticks(np.arange(1.0, 2.0, 0.2))
#else:
#    ax4.set_ylim([0.9,1.9])
#ax4.set_xlim([0.5,2.0])

ax1.tick_params(top=True,bottom=True,left=True,direction='in',length=4)
ax2.tick_params(direction='in',length=4)
ax3.tick_params(top=True,bottom=True,left=True,direction='in',length=4)
#ax4.tick_params(direction='in',length=4)

#xlabel_locs = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
#if subfig=='a':
#    xlabels = []
#elif subfig=='b':
#    xlabels = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
#else:
#    xlabels = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
#plt.xticks(xlabel_locs, xlabels)

if subfig=='a':
    ax3.set_xlabel('',fontsize=16)
elif subfig=='b':
    ax3.set_xlabel('z [m]',fontsize=16)
else:
    ax3.set_xlabel('z [m]',fontsize=16)

# copy and past these two lines into console to remove
# x-axis title and tick marks after initial plot is made
#ax3.set_xlabel('',fontsize=16) # remove x-axis label
#ax3.get_xaxis().set_ticks([]) # remove x-axis ticks

# I know it appears redundant, but keep these here!
#ax1.set_xlim([0.0,2.0])
#ax3.set_xlim([0.0,2.0])

ax1.set_position([0.125,0.550,0.75,0.425])
ax2.set_position([0.125,0.550,0.75,0.425])
ax3.set_position([0.125,0.125,0.75,0.425])
#ax4.set_position([0.125,0.125,0.75,0.425])

#figA.tight_layout()
#figA.subplots_adjust(hspace=0)

plt.show()