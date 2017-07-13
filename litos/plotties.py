#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:13:23 2017

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps
import particle_beam_propagation as pbp
import mike_math as mm
from   calc_M import calc_M
import scipy.stats as stats

#def plotties(ebeam,plasma):

# do analysis on beam/plasma
nstep = len(ebeam)
s     = np.zeros(nstep)
beta  = np.zeros(nstep)
rms_x = np.zeros(nstep)
frac  = 0.95

for i in range(0,len(ebeam)):
    s[i]     = ebeam[i]["s"] # m
    beta[i]  = ebeam[i]["beta"] # m
    rms_x[i] = mm.calc_rms(ebeam[i]["x"],frac)/(1e-6) # um


xi  = ebeam[0]["x"] # m
xf  = ebeam[len(ebeam)-1]["x"] # m
xpi = ebeam[0]["xp"] # rad
xpf = ebeam[len(ebeam)-1]["xp"] # rad

frac = 1
rms_xi = mm.calc_rms(xi,frac)
rms_xf = mm.calc_rms(xf,frac)
rms_xpi = mm.calc_rms(xpi,frac)
rms_xpf = mm.calc_rms(xpf,frac)
rms_xxpi = mm.calc_rms(xi*xpi,frac)
rms_xxpf = mm.calc_rms(xf*xpf,frac)

gbi = ebeam[0]["gb"]
gbf = ebeam[len(ebeam)-1]["gb"]
avg_gbi = np.mean(gbi)
avg_gbf = np.mean(gbf)

avg_x2i  = np.mean(xi**2)
avg_xp2i = np.mean(xpi**2)
avg_xxpi = np.mean(xi*xpi)

avg_x2f  = np.mean(xf**2)
avg_xp2f = np.mean(xpf**2)
avg_xxpf = np.mean(xf*xpf)

rms_x_epsi   = avg_gbi*np.sqrt(avg_x2i*avg_xp2i-avg_xxpi**2)
rms_x_betai  = avg_gbi*(rms_xi**2)/rms_x_epsi
rms_x_gammai = avg_gbi*(rms_xpi**2)/rms_x_epsi
rms_x_alphai = -avg_gbi*rms_xxpi/rms_x_epsi

rms_x_epsf   = avg_gbf*np.sqrt(avg_x2f*avg_xp2f-avg_xxpf**2)
rms_x_betaf  = avg_gbf*(rms_xf**2)/rms_x_epsf
rms_x_gammaf = avg_gbf*(rms_xpf**2)/rms_x_epsf
rms_x_alphaf = -avg_gbf*rms_xxpf/rms_x_epsf

kurti = stats.kurtosis(xi,0,True)
kurtf = stats.kurtosis(xf,0,True)


ift  = int(np.argwhere(s>plasma["up_ramp"]["top_loc"])[0])
Tbft = [ebeam[ift]["beta"],\
        ebeam[ift]["alpha"],\
        ebeam[ift]["gamma"]]

wp0  = (5.64e4)*np.sqrt(plasma["npl"][ift]) # rad/s, plasma ang. freq.
kp0  = wp0/nc.c # m^-1, plasma wave number
kb   = kp0/np.sqrt(2*ebeam[ift]["gbC"])
Tmft = [1.0/kb,0,kb]

Mft  = calc_M(Tbft,Tmft)





# make plots


# beta and rms_x evolution through plasma
figA = plt.figure()
ax1  = figA.add_subplot(111)
ax1.plot(s,plasma["npl"]*beta[0]/max(plasma["npl"]),color='g')
ax1.plot(s,beta,color='b')
ax1.set_ylim([0,1.1*beta[0]])
ax1.set_xlabel('s [m]')
ax1.set_ylabel(r'$\beta$ [m]',color='b')
ax1.tick_params('y',colors='b')

ax2  = ax1.twinx()
ax2.plot(s,rms_x,color='r')
ax2.set_ylabel(r'rms($x$) [$\mu$m]',color='r')
ax2.tick_params('y',colors='r')
ax2.set_ylim([0,1.1*rms_x[0]])

plt.title(r'ramp width = 5.0 cm')

figA.tight_layout()
plt.show()


# initial and final beam distribution
edge = (int(np.max([np.abs(np.min(xf/(1e-6))),np.max(xf/(1e-6))])/10)+1)*10
xbins = np.arange(-edge,edge,2*edge/51)

figB, ax3 = plt.subplots(1,1,sharey=True)
ax3.hist(xf/(1e-6),xbins,ls='solid',fc='none',edgecolor='b',\
                label='Final Beam Dist.')
ax3.hist(xi/(1e-6),xbins,ls='dashed',fc='none',edgecolor='r',\
                label='Initial Beam Dist.')
ax3.set_xlabel(r'x [$\mu$m]')
ax3.set_ylabel(r'macro particles')

plt.legend()

plt.title(r'lens-matched beam distribution')

figB.tight_layout()
plt.show()


#    return 0