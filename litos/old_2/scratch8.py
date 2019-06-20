#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 17:20:08 2018

@author: mike
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
import matplotlib.animation as animation
import time as time

# do analysis on beam/plasma
nstep = len(ebeam)
npart = ebeam[0]["npart"]
s     = np.zeros(nstep)
gbC   = np.zeros(nstep)
eps   = np.zeros(nstep)
beta  = np.zeros(nstep)
alpha = np.zeros(nstep)
gamma = np.zeros(nstep)
cent_x  = np.zeros(nstep)
cent_xp = np.zeros(nstep)
rms_x = np.zeros(nstep)
rms_x_eps = np.zeros(nstep)

frac = 1.0

for i in range(0,nstep):
    s[i]     = ebeam[i]["s"] # m\

start_step = 0
stop_step  = nstep

ft_start   = plasma["up_ramp"]["L"] # m
ft_start_step = int(np.argwhere(s>=ft_start)[0])
#start_step = int(np.argwhere(s>=0.00)[0])
#stop_step  = int(np.argwhere(s>=0.50)[0])
#nstep      = int(stop_step-start_step)

npart = ebeam[0]["npart"]
kb0   = np.zeros(npart)
phi0  = np.zeros(npart)
x0    = np.zeros(npart)
Kx    = np.zeros(npart)

wp0   = (5.64e4)*np.sqrt(plasma["npl"][ft_start_step]) # rad/s, plasma ang. freq.
kp0   = wp0/nc.c               # 1/m, plasma wave number
kb0   = kp0/np.sqrt(2*ebeam[ft_start_step]["gb"])
phi0  = np.arctan2(-ebeam[ft_start_step]["xp"],kb0*ebeam[ft_start_step]["x"])
x0    = np.abs(ebeam[nstep-1]["x"]/np.cos(phi0))
Kx    = ebeam[ft_start_step]["gb"]*kb0*x0

K_cent = mm.calc_mean(Kx)
K_rms  = mm.calc_rms(Kx)
    
for i in range(start_step,stop_step):
    s[i]     = ebeam[i]["s"] # m\
    gbC[i]   = ebeam[i]["gbC"]
    eps[i]   = ebeam[i]["eps"] # mm-mrad
    beta[i]  = ebeam[i]["beta"] # m
    alpha[i] = ebeam[i]["alpha"]
    gamma[i] = ebeam[i]["gamma"] # 1/m
    # x-dimension
    ebeam_rms = ba.calc_ebeam_rms(ebeam,i,frac,'x')
    ebeam_cent = ba.calc_ebeam_cent(ebeam,i,frac,'x')
    cent_x[i]  = ebeam_cent["x"]
    cent_xp[i] = ebeam_cent["xp"]
    rms_x[i]     = ebeam_rms["x"]
    rms_x_eps[i] = ebeam_rms["x_eps"]
    # y-dimension
#    ebeam_rms = ba.calc_ebeam_rms(ebeam,i,frac,'y')
#    ebeam_cent = ba.calc_ebeam_cent(ebeam,i,frac,'y')
#    cent_x[i]  = ebeam_cent["y"]
#    cent_xp[i] = ebeam_cent["yp"]
#    rms_x[i]     = ebeam_rms["y"]
#    rms_x_eps[i] = ebeam_rms["y_eps"]

#%% K distribution
figC, (ax7) = plt.subplots(1, sharex=True, sharey=False)

ax7.hist(Kx,range=[160,240])
ax7.set_ylabel(r'<K>',color='k',fontsize=16)
ax7.tick_params('y',colors='k')
#ax7.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax7.set_xlim([0,20])
#ax7.set_ylim([0.001,90])

plt.show()

print('<K>: ',K_cent)
print('rms(K): ',K_rms)


#%% evolution through plasma


figA, (ax1, ax3, ax5) = plt.subplots(3, sharex=True, sharey=False)

ax1.plot(s,cent_x/(1e-6),color='b',linestyle='solid')
ax1.set_ylabel(r'$<\!x\!>$ ($\mu$ m)',color='b',fontsize=16)
ax1.tick_params('y',colors='b')
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax1.set_ylim([0,0.15])
#ax1.set_ylim([0.1,2])
#ax1.set_xlim([0.5,2.0])

ax2  = ax1.twinx()
ax2.plot(s,cent_xp/(1e-3),color='g',linestyle='solid')
#ax2.plot(skb,dndz/kb,color='g',linestyle='-.')
ax2.set_ylabel(r'$<\!x^{\prime}\!>$ (mrad)',color='g',fontsize=16)
ax2.tick_params('y',colors='g')
#ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax2.set_ylim([0,1.2])
#ax2.set_ylim([1e-6,20])
#ax2.set_xlim([0.5,2.0])
    
ax3.plot(s,rms_x/(1e-6),color='k',linestyle='-')
ax3.set_ylabel(r'$\sigma_x$ ($\mu$m)',color='k',fontsize=16)
ax3.tick_params('y',colors='k')
#ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax3.set_ylim([-4.5,4.5])
#ax3.set_ylim([0.001,90])

ax4 = ax3.twinx()
ax4.plot(s,rms_x_eps/(1e-6),color='r',linestyle='--')
ax4.set_ylabel(r'$\varepsilon$ (mm-mrad)',color='r',fontsize=16)
ax4.tick_params('y',colors='r')
#ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax4.set_ylim([0,10])
#ax4.set_ylim([0.001,90])

ax5.plot(s,K_cent,color='k',linestyle='-')
ax5.set_ylabel(r'<K>',color='k',fontsize=16)
ax5.tick_params('y',colors='k')
#ax5.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax5.set_ylim([0,10])
#ax5.set_ylim([0.001,90])

ax6 = ax5.twinx()
ax6.plot(s,K_rms,color='r',linestyle='--')
ax6.set_ylabel(r'$\sigma_K$',color='r',fontsize=16)
ax6.tick_params('y',colors='r')
#ax6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax6.set_ylim([0,10])
#ax6.set_ylim([0.001,90])

ax1.tick_params(top=True,bottom=True,left=True,direction='in',length=4)
ax2.tick_params(direction='in',length=4)
ax3.tick_params(top=True,bottom=True,left=True,direction='in',length=4)
ax4.tick_params(direction='in',length=4)
ax5.tick_params(top=True,bottom=True,left=True,direction='in',length=4)
ax6.tick_params(direction='in',length=4)

#xlabel_locs = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
#if subfig=='a':
#    xlabels = []
#elif subfig=='b':
#    xlabels = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
#else:
#    xlabels = [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
#plt.xticks(xlabel_locs, xlabels)


ax5.set_xlabel('z [m]',fontsize=16)

# copy and past these two lines into console to remove
# x-axis title and tick marks after initial plot is made
#ax3.set_xlabel('',fontsize=16) # remove x-axis label
#ax3.get_xaxis().set_ticks([]) # remove x-axis ticks

# I know it appears redundant, but keep these here!
#ax1.set_xlim([0.0,2.0])
#ax3.set_xlim([0.0,2.0])

#ax1.set_position([0.125,0.550,0.75,0.425])
#ax2.set_position([0.125,0.550,0.75,0.425])
#ax3.set_position([0.125,0.125,0.75,0.425])
#ax4.set_position([0.125,0.125,0.75,0.425])
#ax5.set_position([0.125,0.125,0.75,0.425])
#ax6.set_position([0.125,0.125,0.75,0.425])

#figA.tight_layout()
#figA.subplots_adjust(hspace=0)

plt.show()

#%%

# Animated plot

def update_plot(step, ebeam, scatxy, scatxxp):
    step = int(np.round(step))
    npart = ebeam[step]['npart']
    x   = ebeam[step]['x']/(1e-6)
    xp  = ebeam[step]['xp']/(1e-3)
    y   = ebeam[step]['y']/(1e-6)
    yp  = ebeam[step]['yp']/(1e-3)
    z   = ebeam[step]['z']/(1e-6)
    gb  = ebeam[step]['gb']
    xy  = np.zeros((npart,2))
    xxp = np.zeros((npart,2))
    yyp = np.zeros((npart,2))
    for i in range(0,npart):
        xy[i][:]  = [x[i],y[i]]
        xxp[i][:] = [x[i],xp[i]]
        yyp[i][:] = [y[i],yp[i]]
    scatxy.set_offsets(xy)
    scatxy.set_array(gb)
    scatxxp.set_offsets(xxp)
    scatxxp.set_array(gb)
    time.sleep(0.02)
    return scatxy, scatxxp,



figB, (ax5,ax6) = plt.subplots(2, sharex=True, sharey=False)

#scat = plt.scatter(xf,xpf,c=gbf,marker='.')
scatxy= ax5.scatter([],[],marker='.')
plt.jet()
ax5.set_xlabel(r'$x$ ($\mu$m)')
ax5.set_ylabel(r'$y$ ($\mu$m)')
ax5.set_xlim([-10,10])
ax5.set_ylim([-10,10])
#plt.xlim((1.1*min(xf),1.1*max(xf)))
#plt.ylim((1.1*min(xpf),1.1*max(xpf)))

scatxxp = ax6.scatter([],[],marker='.')
plt.jet()
ax6.set_xlabel(r'$x$ ($\mu$m)')
ax6.set_ylabel(r'$x^{\prime}$ (mrad)')
ax6.set_xlim([-10,10])
ax6.set_ylim([-1,1])

ani = animation.FuncAnimation(figB, update_plot,\
                              frames=np.linspace(start_step,stop_step,nstep),\
                              fargs=(ebeam, scatxy, scatxxp),blit=True)
plt.show()
