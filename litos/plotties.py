#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:13:23 2017

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



def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)



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
    beta[i]  = ebeam[i]["beta"]/(1e-2) # m
    alpha[i] = ebeam[i]["alpha"]
    gamma[i] = ebeam[i]["gamma"] # 1/m
    ebeam_rms = pb.calc_ebeam_rms(ebeam,i,frac)
    rms_x[i]     = ebeam_rms["x"]/(1e-6)
    rms_x_eps[i] = ebeam_rms["x_eps"]/(1e-6)
    ebeam_kurt = pb.calc_ebeam_kurt(ebeam,plasma,ebeam_rms,i,frac)
    rn_kurt[i] = ebeam_kurt["rn"]
#    x_kurt[i]    = stats.kurtosis(ebeam[i]["x"],0,True)
    
    v_beta[i] = vbeam[i]["beta"]/(1e-2)
    vbeam_rms = pb.calc_ebeam_rms(vbeam,i,frac)
    v_rms_x[i]   = vbeam_rms["x"]/(1e-6)

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
    

i_flat_start = np.argwhere(plasma["s"]>=plasma["up_ramp"]["top_loc"])[0][0]
i_flat_stop  = np.argwhere(plasma["s"]<=\
                           plasma["up_ramp"]["top_loc"]+\
                           plasma["bulk"]["L"])[-1][0]
i_beta_min = np.argmin(beta)
i_beta_end = np.argwhere(beta<=beta[0])[-1][0]

i_rms_x_end = np.argwhere(rms_x<=rms_x[0])[-1][0]

xi  = ebeam[0]["x"]/(1e-6) # um
xf  = ebeam[i_beta_end]["x"]/(1e-6) # um
xpi = ebeam[0]["xp"] # rad
xpf = ebeam[i_rms_x_end]["xp"] # rad
#
#frac = 1.00
#rms_xi = mm.calc_rms(xi,frac)
#rms_xf = mm.calc_rms(xf,frac)
#rms_xpi = mm.calc_rms(xpi,frac)
#rms_xpf = mm.calc_rms(xpf,frac)
#rms_xxpi = mm.calc_rms(xi*xpi,frac)
#rms_xxpf = mm.calc_rms(xf*xpf,frac)
#
#gbi = ebeam[0]["gb"]
#gbf = ebeam[len(ebeam)-1]["gb"]
#avg_gbi = np.mean(gbi)
#avg_gbf = np.mean(gbf)
#
#avg_x2i  = np.mean(xi**2)
#avg_xp2i = np.mean(xpi**2)
#avg_xxpi = np.mean(xi*xpi)
#
#avg_x2f  = np.mean(xf**2)
#avg_xp2f = np.mean(xpf**2)
#avg_xxpf = np.mean(xf*xpf)
#
#rms_x_epsi   = avg_gbi*np.sqrt(avg_x2i*avg_xp2i-avg_xxpi**2)
#rms_x_betai  = avg_gbi*(rms_xi**2)/rms_x_epsi
#rms_x_gammai = avg_gbi*(rms_xpi**2)/rms_x_epsi
#rms_x_alphai = -avg_gbi*rms_xxpi/rms_x_epsi
#
#rms_x_epsf   = avg_gbf*np.sqrt(avg_x2f*avg_xp2f-avg_xxpf**2)
#rms_x_betaf  = avg_gbf*(rms_xf**2)/rms_x_epsf
#rms_x_gammaf = avg_gbf*(rms_xpf**2)/rms_x_epsf
#rms_x_alphaf = -avg_gbf*rms_xxpf/rms_x_epsf
#
#kurti = stats.kurtosis(xi,0,True)
#kurtf = stats.kurtosis(xf,0,True)
#
#
#ift  = int(np.argwhere(s>plasma["up_ramp"]["top_loc"])[0])
#fft  = int(np.argwhere(s>=plasma["up_ramp"]["top_loc"]+plasma["bulk"]["L"])[0])
#
#Tbft = [ebeam[ift]["beta"],\
#        ebeam[ift]["alpha"],\
#        ebeam[ift]["gamma"]]
#
#wp0  = (5.64e4)*np.sqrt(plasma["npl"][ift]) # rad/s, plasma ang. freq.
#kp0  = wp0/nc.c # m^-1, plasma wave number
#kb   = kp0/np.sqrt(2*ebeam[ift]["gbC"])
#Tmft = [1.0/kb,0,kb]
#
#Mft  = calc_M(Tbft,Tmft)
#
#sigr = np.sqrt(eps[0]/(gbC*kb))
#lam = np.pi/(kb*(gbC**2))
#print('sigma: ',sigr)
#print('lambda: ',lam)
#print('Lr: ',2*np.pi*(sigr**2)/lam)
#print('Lg: ',(2*np.pi/kb)/(4*np.pi*np.sqrt(3)*(1e-3)))
#print('eps_lim: ',gbC*lam/(4*np.pi))
#print('lambda_beta: ',2*np.pi/kb)





# make plots


# beta and rms_x evolution through plasma
#figA = plt.figure()
figA, ax1 = plt.subplots(figsize=(8,4))
#figA.subplots_adjust(right=0.5)
#figA.subplots_adjust(left=-0.5)

norm1 = min(v_beta)
#ax1  = figA.add_subplot(111)
#ax1.plot(s,plasma["npl"]*1.5*norm1/max(plasma["npl"]),color='g')
ax1.plot(s,v_beta,color='b',linestyle='dashed')
ax1.plot(s,beta,color='b',linestyle='solid')
#ax1.set_ylim([0,2.0*norm1])
ax1.set_ylim([0,20])
ax1.set_xlabel('s [m]')
ax1.set_ylabel(r'$\beta$ [cm]',color='b')
ax1.tick_params('y',colors='b')

norm2 = min(v_rms_x)
ax2  = ax1.twinx()
ax2.plot(s,v_rms_x,color='r',linestyle='dashed')
ax2.plot(s,rms_x,color='r',linestyle='solid')
ax2.set_ylabel(r'$\sigma_{r,{\rm rms}}$ [$\mu$m]',color='r')
ax2.tick_params('y',colors='r')
#ax2.set_ylim([0,1.5*norm2])
ax2.set_ylim([0,8])

ax3  = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.15))
make_patch_spines_invisible(ax3)
ax3.spines["right"].set_visible(True)
ax3.plot(s,rms_x_eps,color='k',linestyle='solid')
ax3.set_ylabel(r'$\varepsilon_{\rm rms}$ [mm-mrad]',color='k')
ax3.tick_params('y',colors='k')
#ax3.set_ylim([0.9*min(rms_x_eps),1.1*max(rms_x_eps)])
ax3.set_ylim([4.0,8.0])

norm4 = plasma["bulk"]["npl0"]
npl = plasma["npl"]/norm4
ax4  = ax1.twinx()
ax4.spines["left"].set_position(("axes", -0.20))
make_patch_spines_invisible(ax4)
ax4.spines["left"].set_visible(True)
ax4.yaxis.set_label_position('left')
ax4.spines["left"].set_visible(True)
ax4.yaxis.set_label_position('left')
ax4.yaxis.set_ticks_position('left')
ax4.plot(s,npl,color='g',linestyle='solid')
ax4.set_ylabel(r'$n_p$ [$10^{17}\,{\rm cm}^{-3}$]',color='g')
ax4.tick_params('y',colors='g')
#ax4.set_ylim([0,1.5*max(npl)])
ax4.set_ylim([0,1.4])

ax4.text(0.50, 0.75, r'$n_p \times 10$',
        verticalalignment='center', horizontalalignment='center',
        transform=ax4.transAxes,
        color='green', fontsize=12)

#ax4  = ax1.twinx()
#ax4.plot(s,rn_kurt,color='k',linestyle='dashed')
#ax4.set_ylabel(r'kurtosis',color='k')
#ax4.tick_params('y',colors='k')
#ax4.set_ylim([1.1*min(rn_kurt),1.1*max(rn_kurt)])

npl0 = plasma["bulk"]["npl0"]
#plt.title(r'beam matching into %.1e $cm^{-3}$ plasma source'%npl0)

figA.tight_layout()
plt.subplots_adjust(left=0.225)
plt.subplots_adjust(right=0.800)
plt.show()






# initial and final beam distribution

#edge = (int(np.max([np.abs(np.min(xf/(1e-6))),np.max(xf/(1e-6))])/10)+1)*10
edge = 80
nbins = 25
xbins = np.arange(-edge,edge,2*edge/nbins)

ki = stats.kurtosis(xi,0,True)
kf = stats.kurtosis(xf,0,True)

figB, ax3 = plt.subplots(1,1,sharey=True)
hi = ax3.hist(xf,xbins,ls='-.',fc='none',edgecolor='b',\
              normed=True,label='Final Beam Dist.\n kurtosis = %.1f'%ki)
hf = ax3.hist(xi,xbins,ls='dashed',fc='none',edgecolor='r',\
              normed=True,label='Initial Beam Dist.\n kurtosis = %.1f'%kf)




# best fit of data
(mui, sigmai) = stats.norm.fit(xi)
# add a 'best fit' line
xx = np.linspace(-edge,+edge,200)
yy = mlab.normpdf( xx, mui, sigmai)
li = plt.plot(xx, yy, 'r--', linewidth=1)

# best fit of data
(muf, sigmaf) = stats.norm.fit(xf)
# add a 'best fit' line
xx = np.linspace(-edge,+edge,200)
yy = mlab.normpdf( xx, muf, sigmaf)
lf = plt.plot(xx, yy, 'b-.', linewidth=1)


ax3.set_xlim([-edge,+edge])
ax3.set_xlabel(r'r [$\mu$m]')
ax3.set_ylabel(r'fraction of macro particles / %.1f $\mu$m'%(2*edge/nbins))

plt.legend()

#plt.title(r'lens-matched beam distribution')

figB.tight_layout()
plt.show()




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