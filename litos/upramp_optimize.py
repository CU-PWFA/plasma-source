#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 17:00:28 2017

@author: litos
"""

import numpy as np
from scipy.optimize import minimize
from fitupramp import fitupramp

# constants
c  = 3e8 # m/s
me = 0.511e-3 # GeV

# define initial beam
gb0    = 20000
eps0   = 5e-6 # m-rad
beta0  = 0.10 # m
alpha0 = 0.00
gamma0 = (1+alpha0**2)/beta0 # 1/m
dE0    = 0.00
npart  = 0
beam0  = [gb0,eps0,beta0,alpha0,gamma0,dE0,npart]

# define plasma up-ramp
npl0     = 1e18 # cm^-3
shape    = 'gauss'

## define target values
wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
kb0 = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number
targ_beta  = 1.0/kb0
targ_alpha = 0.0
targ_gamma = (1.0+targ_alpha**2)/targ_beta

## initial guess
# x = [waist location, ramp ~half-length, ramp exponent]
waist0 = -0.38 # m, w.r.t. Lp
lp0    = 0.12 # m
P0     = 2.00
x0 = [waist0,lp0,P0]

# bounds on variables
bnds = [(-1.0,+0.0),(0.03,0.50),(2.00,2.00)]

# constraints on variables
#cons = ({'type': 'ineq', 'fun': lambda x:  x[1] - 5*x[2]})#,\
#        {'type': 'ineq', 'fun': lambda x:  x[1] -   x[0]})

# Lp : full length of plasma ramp
Lp0    = 1.00 # m
s = np.linspace(0,Lp0,1001)

## constant value arguments
args = [gb0,eps0,beta0,alpha0,gamma0,\
        npl0,shape,s,Lp0,\
        targ_beta,targ_alpha,targ_gamma]

## perform fit
fitres = minimize(fitupramp,x0,(args,),\
                  bounds=bnds,\
                  method='L-BFGS-B',\
                  options={'disp':True})

xfit = fitres.x
print("Fit results: %s" % str(xfit))

beam0 = [gb0,eps0,beta0,alpha0,gamma0,0,0]
[waist,lp,P] = xfit
s_w   = Lp0 + waist # m
beam0 = propbeamlast(beam0,[0,-s_w])
s     = np.linspace(0,Lp0,round(Lp0/0.001))
pargs = [Lp0,lp,2.0,1e-3]
plas  = plasma_ramp(npl0,shape,s[0:len(s)-1],pargs)
print(plas[len(plas)-1])
beam  = propbeamlast(beam0,s,plas)
[gb,eps,beta,alpha,gamma,dE,npart] = beam
Tb = [beta,alpha,gamma]
Tm = [targ_beta,targ_alpha,targ_gamma]
M = calc_M(Tb,Tm)
print(Tb)
print(Tm)
print(M)
