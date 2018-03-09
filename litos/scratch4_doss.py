#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 20:53:15 2018

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
from matplotlib.ticker import MaxNLocator # added
from matplotlib import gridspec


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

frac = 1.0

for i in range(0,nstep):
    s[i]     = ebeam[i]["s"] # m\
    gbC[i]   = ebeam[i]["gbC"]
    eps[i]   = ebeam[i]["eps"] # mm-mrad
    beta[i]  = ebeam[i]["beta"]/(1e-2) # m
    alpha[i] = ebeam[i]["alpha"]
    gamma[i] = ebeam[i]["gamma"] # 1/m
    ebeam_rms = ba.calc_ebeam_rms(ebeam,i,frac)
    ebeam_cent = ba.calc_ebeam_cent(ebeam,i,frac)
    cent_x[i]  = ebeam_cent["x"]/(1e-6)
    cent_xp[i] = ebeam_cent["xp"]/(1e-6)
    rms_x[i]     = ebeam_rms["x"]/(1e-6)
    rms_x_eps[i] = ebeam_rms["x_eps"]/(1e-6)
    
    [u,v] = ba.real2norm_coords(ebeam[i]["x"],ebeam[i]["xp"],\
                                ebeam_rms["x_beta"],ebeam_rms["x_alpha"])
    J = (u**2+v**2)/2
    phi = np.arctan2(v,u)
    un = np.sqrt(2*J*ebeam[i]["gbC"])*np.cos(phi)/np.sqrt(rms_x_eps[i])
    un_kurt[i] = stats.kurtosis(un,0,True)
    J_kurt[i] = stats.kurtosis(J,0,False,True)
    
    
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
print('Bmag: ',BB)
    
#%% beta and rms_x evolution through plasma

figA = plt.figure(figsize=(8,6))

ax3 = figA.add_subplot(2,1,2)
ax3.plot(s,rms_x_eps/rms_x_eps[0],color='k',linestyle='-')
ax3.plot(s,BB*np.ones(len(s)),color='k',linestyle='-.')
ax3.set_ylabel(r'$\varepsilon_n/\varepsilon_{n,0}$',color='k',fontsize=16)
ax3.tick_params('y',colors='k')
#ax3.set_ylim([0.9875,1.0125]) # matched limits
ax3.set_ylim([0.9,1.9]) # mismatched limits
ax3.set_xlim([0.5,2.0])

ax4 = ax3.twinx()
ax4.plot(s,J_kurt/J_kurt[0],color='r',linestyle='--')
ax4.set_ylabel(r'$\kappa/\kappa_{0}$',color='r',fontsize=16)
ax4.tick_params('y',colors='r')
#ax4.set_ylim([0.9875,1.0125]) # matched limits
ax4.set_ylim([0.9,1.9]) # mismatched limits
ax4.set_xlim([0.5,2.0])

xlabel_locs = [0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]
xlabels = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5]
plt.xticks(xlabel_locs, xlabels)

ax3.set_xlabel('z [m]',fontsize=16)
ax3.tick_params(direction='inout')

ax1 = figA.add_subplot(2,1,1)
ax1.plot(s,rms_x,color='b',linestyle='solid')
ax1.set_ylabel(r'$\sigma_x$ [$\mu m$]',color='b',fontsize=16)
ax1.tick_params('y',colors='b')
ax1.set_ylim([0,6.5])
ax1.set_xlim([0.5,2.0])

npl = plasma["npl"]/plasma["bulk"]["npl0"]
ax2  = ax1.twinx()
ax2.plot(s,npl,color='g',linestyle='solid')
ax2.set_ylabel(r'$n_p/n_{p,0}$',color='g',fontsize=16)
ax2.tick_params('y',colors='g')
ax2.set_ylim([0,1.3])
ax2.set_xlim([0.5,2.0])

# plasma density text
#ax2.text(0.50, 0.85, r'$n_{p,0} = %2.1e \,{\rm cm^{-3}}$'%plasma["bulk"]["npl0"],
#        verticalalignment='center', horizontalalignment='center',
#        transform=ax2.transAxes,
#        color='green', fontsize=16)
ax2.text(0.50, 0.85, r'$n_{p,0} = 5\times 10^{16} \,{\rm cm^{-3}}$',
        verticalalignment='center', horizontalalignment='center',
        transform=ax2.transAxes,
        color='green', fontsize=16)
ax2.text(0.08, 0.95, '(a)',
        verticalalignment='center', horizontalalignment='center',
        transform=ax2.transAxes,
        color='black', fontsize=16)
xlabel_locs = [0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]
xlabels = []#[0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5]
plt.xticks(xlabel_locs, xlabels)

# copy and past these two lines into console to remove
# x-axis title and tick marks after initial plot is made
#ax3.set_xlabel('',fontsize=16) # remove x-axis label
#ax3.get_xaxis().set_ticks([]) # remove x-axis ticks

# I know it appears redundant, but keep these here!
ax1.set_xlim([0.5,2.0])
ax3.set_xlim([0.5,2.0])

#figA.tight_layout()

figA.subplots_adjust(hspace=0)
#xticklabels= ax1.get_xticklabels() + ax3.get_xticklabels()
#plt.setp(xticklabels, visible=True)
plt.show()