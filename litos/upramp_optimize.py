#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 17:00:28 2017

@author: litos
"""

import numpy as np
from scipy.optimize import minimize
from match_upramp import match_upramp
from mc_beam import mc_beam

# constants
c  = 3e8 # m/s
me = 0.511e-3 # GeV

# define plasma bulk (flat-top) properties
npl0   = 1e17 # cm^-3
L_ft   = 0.20 # m
dEds0  = 6.00 # GeV/m
dgds0  = dEds0/me

# define plasma up-ramp
shape_up = 'gauss'
L_up     = 1.00 # m, full length of propagation

# definte longitudinal propagation steps
s = np.linspace(0,L_up+L_ft,round((L_up+L_ft)/(1e-3))+1) # m, propagation steps

# define initial beam
gb0    = 20000
eps0   = 5e-6 # m-rad
beta0  = 0.10 # m
alpha0 = 0.00
gamma0 = (1+alpha0**2)/beta0 # 1/m
dE0    = 0.01
npart  = 100
beam0  = [gb0,eps0,beta0,alpha0,gamma0,dE0,npart]

# generate particles
[x,xp,ux,vx,gb] = mc_beam(beam0)
[y,yp,uy,vy,gb] = mc_beam(beam0)
z = np.zeros(npart)
parts0 = np.zeros((npart,6))
for i in range(0,npart):
    parts0[i][:] = [x[i],xp[i],y[i],yp[i],z[i],gb[i]]

## minimization variables, initial guess
waist0 = -0.30 # m, waist location w.r.t. flat top start
hw_up0 = 0.15 # m, half-width of up-ramp = sqrt(2*log(2))*sigma_gauss
x0 = [waist0,hw_up0] # initial guess
bnds = [(-1.0,+0.0),(0.03,0.50)] # bounds on variables

## constant value arguments
args = [parts0,npl0,dgds0,shape_up,s,L_up]

## perform fit
fitres = minimize(match_upramp,x0,(args,),\
                  bounds=bnds,\
                  method='L-BFGS-B',\
                  options={'disp':True,'eps': 1e-3})

xfit = fitres.x
print("Fit results: %s" % str(xfit))



"""

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

"""