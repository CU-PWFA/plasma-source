#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:39:58 2018

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

# do analysis on beam/plasma
nstep = len(ebeam)
npart = ebeam[0]["npart"]
s     = np.zeros(nstep)
gbC   = np.zeros(nstep)
eps   = np.zeros(nstep)
beta  = np.zeros(nstep)
alpha = np.zeros(nstep)
gamma = np.zeros(nstep)


bm = np.zeros(nstep)

beta_w = 0.20

gbC  = ebeam[1]["gbC"]

npl0 = plasma["bulk"]["npl0"]
wp0    = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq. (flat-top)
kp0    = wp0/nc.c               # 1/m, plasma wave number (flat-top)
kb0    = kp0/np.sqrt(2*gbC)     # 1/m, betatron wave number (flat-top)

hw  = plasma["up_ramp"]["hw"]
sig = hw/(np.sqrt(2*np.log(2)))

s0 = plasma["up_ramp"]["top_loc"]

def beta_func(s,a):
    return (1./kb0) + (((s-s0)**2)/beta_w + beta_w - 1./kb0)*np.arctan( (np.abs(s-s0)/(a*sig)) )*2/np.pi \
            - (beta_w-1./kb0)*(s-s0)/(5*np.pi*sig)

for i in range(0,nstep):
    s[i]     = ebeam[i]["s"] # m
    beta[i]  = ebeam[i]["beta"] # m
#    bm[i]    = 1./kb0 + ((s[i]-s0)**2)/beta_w + beta_w*(1.-np.exp( -((s[i]-s0)**2)/(4*sig**2) ))
#    bm[i]    = 1./kb0 + (((s[i]-s0)**2)/beta_w + beta_w)*(1.-np.exp( -((s[i]-s0)**2)/(150*sig**2) ))
#    bm[i]    = 1./kb0 + (((s[i]-s0)**2)/beta_w + beta_w)*(1.-np.exp( -(np.abs(s[i]-s0)/(25*sig)) ) )
#    bm[i]    = 1./kb0 + (((s[i]-s0)**2)/beta_w + beta_w)*(2/(1+np.exp( -(np.abs(s[i]-s0)/(25*sig)) ) )-1)
#    bm[i]    = 1./kb0 + (((s[i]-s0)**2)/beta_w + beta_w)*np.tanh( (np.abs(s[i]-s0)/(25*sig)) )
#    bm[i]    = 1./kb0 + (((s[i]-s0)**2)/beta_w + beta_w-1./kb0)*np.arctan( (np.abs(s[i]-s0)/(10*sig)) )*2/np.pi
#    bm[i]    = (1./kb0)*(1-np.arctan( (np.abs(s[i]-s0)/(10*sig)) )*2/np.pi) \
#               + (((s[i]-s0)**2)/beta_w + beta_w)*np.arctan( (np.abs(s[i]-s0)/(10*sig)) )*2/np.pi \
#               - (beta_w-1./kb0)*(s[i]-s0)/(5*np.pi*sig)
              
    bm[i]    = beta_func(s[i],10)

 #   bm[i]    = 1./kb0 + ((s[i]-s0)**2)/beta_w

popt, pcov = opt.curve_fit(beta_func,s,beta)

print(popt)
print(pcov)


figA, axA = plt.subplots(1,1,sharey=True)
plt_beta = axA.plot(s, beta, 'b--', linewidth=1)
plt_bm = axA.plot(s, bm, 'r', linewidth=1)
plt.xlabel(r'z (m)')
plt.ylabel(r'$\beta$ (m)')
#plt_beta.set_label('$\beta$ from sim')
axA.legend([r'$\beta$ from sim.','f(z)'])

