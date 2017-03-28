#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:59:44 2017

@author: litos
"""

## common libraries
#import sys
from numpy import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt

## custom functions
from fitramp import fitramp
from propbeam import propbeam
from plasma_ramp import plasma_ramp
from draw_ellipse import draw_ellipse
from calc_Bmag import calc_Bmag

## define constants
c = 3e8 # m/s, speed of light

## set beam parameters
# use values for nominal vacuum waist
gb0   = 10000 # Lorentz factor
eps0   = 5e-6 # m-rad
beta0  = 0.10 # m
alpha0 = 0.0
gamma0 = (1+(alpha0**2))/beta0 # 1/m
T0     = [beta0,alpha0,gamma0] # Twiss vector
dgb    = 0.01 # energy spread

## set plasma parameters
# use values for main accelerator plasma
np0 = 5e16 # cm^-3, plasma number density
wp0 = (5.64e4)*sqrt(np0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave nutarget_betaer
kb0 = kp0/sqrt(2*gb0) # m^-1, betatron wave number
shape = 'gauss'
#shape = 'gen_gauss'
#shape = 'trap'
#shape = 'sigmoid'
#shape = 'gomp'
#shape = 'xu3'
#shape = 'xu4'
#shape = 'xu5'

#-------------------------
## optimize

## parameters to optimize:
# waist : location of vacuum waist in z
# lp : characteristic ramp length

## define target values
targ_beta  = 1/kb0
targ_alpha = 0
targ_gamma = kb0

## initial guess
# x = [waist location, ramp ~half-length, ramp exponent]
waist0 = -0.40 # m, w.r.t. Lp
lp0    = 0.15 # m
P0     = 2.00
x0 = [waist0,lp0,P0]

# bounds on variables
bnds = [(-1.0,+1.0),(0.03,0.50),(2.00,2.00)]

# constraints on variables
#cons = ({'type': 'ineq', 'fun': lambda x:  x[1] - 5*x[2]})#,\
#        {'type': 'ineq', 'fun': lambda x:  x[1] -   x[0]})

# Lp : full length of plasma ramp
Lp0    = 1.00 # m
z = linspace(0,2.0,2001)

## constant value arguments
args = [gb0,eps0,beta0,alpha0,gamma0,\
        np0,shape,z,Lp0,\
        targ_beta,targ_alpha,targ_gamma]

## perform fit
fitres = minimize(fitramp,x0,(args,),\
                  bounds=bnds,\
                  method='L-BFGS-B',\
                  options={'disp':False})
#                  constraints=cons,\

xfit = fitres.x

print("Fit results: %s" % str(xfit))

waistfit = xfit[0]
abs_waistfit = Lp0 + waistfit
beam0 = [gb0,eps0,beta0,alpha0,gamma0]
beam0 = propbeam(beam0,[0,-1*abs_waistfit])
beam0 = beam0[len(beam0)-1][:]

lpfit = xfit[1]
Pfit  = xfit[2]
pargs = [Lp0,lpfit,Pfit,1e-3]
#npl   = makeplasma(np0,shape,z[0:len(z)-1],pargs)
npl   = plasma_ramp(np0,shape,z[0:len(z)-1],pargs)

beam  = propbeam(beam0,z,npl)


# make energy slice ellipses
ngb = 11
igb = gb0*linspace((1-dgb),(1+dgb),ngb)

ibeam0 = zeros((ngb,5))
ibeam  = zeros((ngb,5))
iTbeam = zeros((ngb,3))
ixb    = zeros((ngb,2*101))
ixpb   = zeros((ngb,2*101))
ixpnb  = zeros((ngb,2*101))
for i in range(0,len(igb)):
    ib0 = [igb[i],eps0,beta0,alpha0,gamma0]
    ib0 = propbeam(ib0,[0,-1*abs_waistfit])
    ibeam0[i][:] = ib0[len(ib0)-1][:]
    ib  = propbeam(ibeam0[i][:],z,npl)
    ibeam[i][:] = ib[len(ib)-1][:]
    iTbeam[i][:] = [ibeam[i,2],ibeam[i,3],ibeam[i,4]]
    [ixb[i][:],ixpb[i][:],ixpnb[i][:]] = \
     draw_ellipse(eps,iTbeam[i][:],Ttarg,False,101)



## plot beam evolution
ymax_plt = 10*beam[len(beam)-1,2]
norm = ymax_plt/max(npl)
plt.plot([z],[beam[:,2]],'b.',[z],[npl*norm],'g.')
plt.xlabel(r'z (m)')
plt.ylabel(r'$\beta$ (m)')
plt.text(0.5*max(z),0.5*ymax_plt,\
         "lp = %.2f, waist = %.2f" % (lpfit,waistfit),\
         verticalalignment='center',\
         horizontalalignment='right')
plt.ylim((0,ymax_plt))
plt.axhline(targ_beta,color='r')
plt.show()

## plot beam ellipse
ifin  = len(beam)-1
eps   = beam[ifin,1]
Tbeam = [beam[ifin,2],beam[ifin,3],beam[ifin,4]]
Ttarg = [targ_beta,targ_alpha,targ_gamma]
Bmag = calc_Bmag(Tbeam,Ttarg)
[xb,xpb,xpnb] = draw_ellipse(eps,Tbeam,Ttarg,False)
[xt,xpt,xpnt] = draw_ellipse(eps,Ttarg,Ttarg,False)
# plot energy slices
for i in range(0,len(ixb)):
    plt.plot(ixb[i][:],ixpb[i][:],'.')
# plot matched ellipse
plt.plot(xt,xpt,'k.')
plt.title('Beam phase space ellipse')
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
#plt.axes().set_aspect('equal','datalim')
plt.text(0,0,"Bmag = %.3f" % Bmag,\
         verticalalignment='center',\
         horizontalalignment='center')
plt.show()

## plot normalized beam ellipse
# plot energy slices
for i in range(0,len(ixb)):
    plt.plot(ixb[i][:],ixpnb[i][:],'.')
# plot matched ellipse
plt.plot(xt,xpnt,'k.')
plt.title('Normalized beam phase space ellipse')
plt.xlabel(r'$\bar{x}$ (m)')
plt.ylabel(r'$\bar{x^{\prime}}$ (m)')
plt.axes().set_aspect('equal','datalim')
plt.text(0,0,"Bmag = %.3f" % Bmag,\
         verticalalignment='center',\
         horizontalalignment='center')
plt.show()
