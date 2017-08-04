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
npart = ebeam[0]["npart"]
s     = np.zeros(nstep)
eps   = np.zeros(nstep)
beta  = np.zeros(nstep)
alpha = np.zeros(nstep)
gamma = np.zeros(nstep)
rms_x = np.zeros(nstep)
rms_x_eps = np.zeros(nstep)
x_kurt  = np.zeros(nstep)
xp_kurt = np.zeros(nstep)
rn_kurt = np.zeros(nstep)


v_beta = np.zeros(nstep)
v_rms_x = np.zeros(nstep)

Psi_x = np.zeros([nstep,npart])
K_x   = np.zeros([nstep,npart])

PsiC_x = np.zeros(nstep)
KC_x   = np.zeros(nstep)
Krms_x = np.zeros(nstep)

K_twiss = np.zeros(nstep)

frac  = 1.00

for i in range(0,nstep):
    s[i]     = ebeam[i]["s"] # m
    eps[i]   = ebeam[i]["eps"] # mm-mrad
    beta[i]  = ebeam[i]["beta"] # m
    alpha[i] = ebeam[i]["alpha"]
    gamma[i] = ebeam[i]["gamma"] # 1/m
    ebeam_rms = pb.calc_ebeam_rms(ebeam,i,frac)
    rms_x[i]     = ebeam_rms["x"]
    rms_x_eps[i] = ebeam_rms["x_eps"]
    ebeam_kurt = pb.calc_ebeam_kurt(ebeam,plasma,i,frac)
    rn_kurt[i] = ebeam_kurt["rn"]
#    x_kurt[i]    = stats.kurtosis(ebeam[i]["x"],0,True)
    
    v_beta[i] = vbeam[i]["beta"]
    vbeam_rms = pb.calc_ebeam_rms(vbeam,i,frac)
    v_rms_x[i]   = vbeam_rms["x"]

    wp  = (5.64e4)*np.sqrt(plasma["npl"][i]) # rad/s, plasma ang. freq.
    kp  = wp/nc.c # m^-1, plasma wave number
    
    for j in range(0,npart):
        x      = ebeam[i]["x"][j]
        xp     = ebeam[i]["xp"][j]
        gb     = ebeam[i]["gb"][j]
        kb     = kp/np.sqrt(2*gb)
        beta_m = 1.0/kb
        R_x    = np.sqrt(x**2 + (beta_m*xp)**2)
        Psi_x[i,j] = np.arctan2(beta_m*xp,x)
        K_x[i,j]   = gb*kb*R_x
    
    # make Psi go from 0 to 2*pi
    Psi_x[(Psi_x<0)] += 2*np.pi

    frac = 1.0
    PsiC_x[i] = mm.calc_mean(Psi_x[i,:],frac)
    KC_x[i]   = mm.calc_mean(K_x[i,:],frac)
    Krms_x[i] = mm.calc_rms(K_x[i,:],frac)

    gbC    = ebeam[i]["gbC"]
    kb     = kp/np.sqrt(2*gbC)
    beta_m = 1.0/kb
    R_x    = np.sqrt(eps[i]*beta[i]/gbC+(beta_m**2)*eps[i]*gamma[i]/gbC)
    K_twiss[i] = gbC*kb*R_x
    

xi  = ebeam[0]["x"] # m
xf  = ebeam[len(ebeam)-1]["x"] # m
xpi = ebeam[0]["xp"] # rad
xpf = ebeam[len(ebeam)-1]["xp"] # rad

frac = 1.00
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
fft  = int(np.argwhere(s>=plasma["up_ramp"]["top_loc"]+plasma["bulk"]["L"])[0])

Tbft = [ebeam[ift]["beta"],\
        ebeam[ift]["alpha"],\
        ebeam[ift]["gamma"]]

wp0  = (5.64e4)*np.sqrt(plasma["npl"][ift]) # rad/s, plasma ang. freq.
kp0  = wp0/nc.c # m^-1, plasma wave number
kb   = kp0/np.sqrt(2*ebeam[ift]["gbC"])
Tmft = [1.0/kb,0,kb]

Mft  = calc_M(Tbft,Tmft)

sigr = np.sqrt(eps[0]/(gbC*kb))
lam = np.pi/(kb*(gbC**2))
print('sigma: ',sigr)
print('lambda: ',lam)
print('Lr: ',2*np.pi*(sigr**2)/lam)
print('Lg: ',(2*np.pi/kb)/(4*np.pi*np.sqrt(3)*(1e-3)))
print('eps_lim: ',gbC*lam/(4*np.pi))
print('lambda_beta: ',2*np.pi/kb)

# make plots


# beta and rms_x evolution through plasma
figA = plt.figure()
ax1  = figA.add_subplot(111)
ax1.plot(s,plasma["npl"]*beta[0]/max(plasma["npl"]),color='g')
ax1.plot(s,v_beta,color='b',linestyle='dashed')
ax1.plot(s,beta,color='b',linestyle='solid')
ax1.set_ylim([0,1.1*beta[0]])
ax1.set_xlabel('s [m]')
ax1.set_ylabel(r'$\beta$ [m]',color='b')
ax1.tick_params('y',colors='b')

ax2  = ax1.twinx()
ax2.plot(s,v_rms_x,color='r',linestyle='dashed')
ax2.plot(s,rms_x,color='r',linestyle='solid')
ax2.set_ylabel(r'$\sigma_r$ [$\mu$m]',color='r')
ax2.tick_params('y',colors='r')
ax2.set_ylim([0,1.1*max(rms_x)])

ax3  = ax1.twinx()
ax3.plot(s,rms_x_eps,color='k',linestyle='solid')
ax3.set_ylabel(r'$\varepsilon$ [mm-mrad]',color='k')
ax3.tick_params('y',colors='k')
ax3.set_ylim([0,1.1*max(rms_x_eps)])

ax4  = ax1.twinx()
ax4.plot(s,rn_kurt,color='k',linestyle='dashed')
ax4.set_ylabel(r'kurtosis',color='k')
ax4.tick_params('y',colors='k')
ax4.set_ylim([0,1.1*max(rn_kurt)])


plt.title(r'$\beta$ and beam size')

figA.tight_layout()
plt.show()


## initial and final beam distribution
#edge = (int(np.max([np.abs(np.min(xf/(1e-6))),np.max(xf/(1e-6))])/10)+1)*10
#xbins = np.arange(-edge,edge,2*edge/51)
#
#figB, ax3 = plt.subplots(1,1,sharey=True)
#ax3.hist(xf/(1e-6),xbins,ls='solid',fc='none',edgecolor='b',\
#                label='Final Beam Dist.')
#ax3.hist(xi/(1e-6),xbins,ls='dashed',fc='none',edgecolor='r',\
#                label='Initial Beam Dist.')
#ax3.set_xlabel(r'x [$\mu$m]')
#ax3.set_ylabel(r'macro particles')
#
#plt.legend()
#
#plt.title(r'lens-matched beam distribution')
#
#figB.tight_layout()
#plt.show()




## K evolution
#
#figC = plt.figure()
#ax4  = figC.add_subplot(111)
##ax4.plot(s,plasma["npl"]*max(test_Psi)/max(plasma["npl"]),color='g')
#
##ax4.plot(s,test_Psi,color='b')
##ax4.plot(s,test2_Psi,color='r')
#
##ax4.scatter(test_Psi,test_xp)
#
##ax4.plot(s,test_x,color='b')
##ax4.plot(s,test2_x,color='r')
#
#ax4.plot(s,KC_x,color='b')
#ax4.plot(s,KC_x+Krms_x,color='r')
#ax4.plot(s,KC_x-Krms_x,color='r')
#ax4.plot(s,K_twiss,color='g')
#
##ax4.hist(K_x[ift,:],21,ls='solid',fc='none',edgecolor='b',\
##                label='K Dist.')
##ax4.hist(K_x[fft,:],21,ls='solid',fc='none',edgecolor='r',\
##                label='K Dist.')
#
##ax4.scatter(ebeam[ift]["x"],K_x[ift,:])
#
##ax4.set_ylim([0,+1.1*2*np.pi])
#
#ax4.set_xlabel('s [m]')
##ax4.set_ylabel(r'$\Psi$ [rad]',color='b')
##ax4.tick_params('y',colors='b')
#
#ax4.set_ylabel(r'K',color='k')
#ax4.tick_params('y',colors='k')
#
#
#plt.title(r'K and/or $\Psi$')
#
#figC.tight_layout()
#plt.show()


#    return 0