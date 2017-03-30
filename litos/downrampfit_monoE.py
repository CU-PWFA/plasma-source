#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 19:36:37 2017

@author: litos
"""

## common libraries
#import sys
from numpy import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt

## custom functions
from fitdownramp import fitdownramp
from propbeam import propbeam
from plasma_ramp import plasma_ramp
from draw_ellipse import draw_ellipse
from calc_Bmag import calc_Bmag
from calc_M import calc_M

## define constants
c = 3e8 # m/s, speed of light

## set plasma parameters
# use values for main accelerator plasma
np0 = 5e16 # cm^-3, plasma number density
wp0 = (5.64e4)*sqrt(np0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
shape = 'gauss'
#shape = 'gen_gauss'
#shape = 'trap'
#shape = 'sigmoid'
#shape = 'gomp'
#shape = 'xu3'
#shape = 'xu4'
#shape = 'xu5'

## set beam parameters
# use values for nominal vacuum waist
gb0    = 24219 # Lorentz factor
eps0   = 1.1*(5e-6) # m-rad
kb0    = kp0/sqrt(2*gb0) # m^-1, betatron wave number
beta0  = (1/kb0) # m
alpha0 = 0.0
gamma0 = (1+(alpha0**2))/beta0 # 1/m
T0     = [beta0,alpha0,gamma0] # Twiss vector

#-------------------------
## optimize

## parameters to optimize:
# waist : location of vacuum waist in z
# lp : characteristic ramp length

## define target values
targ_beta  = 0.15 # m
targ_alpha = 0
targ_gamma = (1+(targ_alpha**2))/targ_beta # 1/m

## initial guess
# x = [waist location, ramp ~half-length, ramp exponent]
waist0 = 0.50 # m, w.r.t. Lp
lp0    = 0.15 # m
P0     = 2.00
x0 = [waist0,lp0,P0]

# bounds on variables
bnds = [(0.50,0.50),(0.10,0.20),(2.00,2.00)]

# constraints on variables
#cons = ({'type': 'ineq', 'fun': lambda x:  x[1] - 5*x[2]})#,\
#        {'type': 'ineq', 'fun': lambda x:  x[1] -   x[0]})

# Lp : full length of plasma ramp
Lp0    = 0.00 # m
z = linspace(0,2.0,2001)

## constant value arguments
args = [gb0,eps0,beta0,alpha0,gamma0,\
        np0,shape,z,Lp0,\
        targ_beta,targ_alpha,targ_gamma]

## perform fit
fitres = minimize(fitdownramp,x0,(args,),\
                  bounds=bnds,\
                  method='L-BFGS-B',\
                  options={'disp':False})
#                  constraints=cons,\

xfit = fitres.x

print("Fit results: %s" % str(xfit))


lpfit = xfit[1]
Pfit  = xfit[2]
pargs = [Lp0,lpfit,Pfit,1e-3]
npl   = plasma_ramp(np0,shape,z[0:len(z)-1],pargs,'down')

waistfit = xfit[0]
abs_waistfit = Lp0 + waistfit

iz_w = np.where(z<=Lp0+3*lp0*2)[0]
z_w = z[iz_w]
beam0 = [gb0,eps0,beta0,alpha0,gamma0]

beam = propbeam(beam0,z_w,npl[0:len(z_w)-1])
[gb,eps,beta,alpha,gamma] = beam[len(beam)-1][:]

print(beta)
print(alpha)
print(gamma)

dgb0 = 1.01*gb0
dbeam0 = [dgb0,eps0,beta0,alpha0,gamma0]
dbeam = propbeam(dbeam0,z,npl[0:len(z)-1])
[dgb,eps,dbeta,dalpha,dgamma] = dbeam[len(dbeam)-1][:]

dT = [dbeta,dalpha,dgamma]
Tnom = [beta,alpha,gamma]
dM = calc_M(dT,Tnom)

print(Tnom)
print(dT)
print(dM)


## plot beam evolution
ymax_plt = 1*beam[len(beam)-1,2]
norm = ymax_plt/max(npl)
plt.figure()
plt.plot([z],[beam[:,2]],'b.',[z],[npl*norm],'g.')
plt.xlabel(r'z (m)')
plt.ylabel(r'$\beta$ (m)')
plt.title(r"Optimized ramp with %s profile" % shape)
plt.text(0.5*max(z),0.5*ymax_plt,\
         r"$\sigma$ = %.2f m, waist = %.2f m" % (lpfit,waistfit),\
         verticalalignment='center',\
         horizontalalignment='center')
plt.ylim((0,ymax_plt))
plt.axhline(targ_beta,color='r')
plt.show()




#
#print(beam[iz_w[len(iz_w)-1]][:])
#
#T_w = [beam[iz_w[len(iz_w)-1]][2],\
#         beam[iz_w[len(iz_w)-1]][3],\
#         beam[iz_w[len(iz_w)-1]][4]]
#
#print(T_w)
#Ttarg = [targ_beta,targ_alpha,targ_gamma]
#M_w = calc_M(T_w,Ttarg)
#
#print(M_w)
#
#
#
## scan ramp lengths, for each ramp length look at beta
## at location of minimum |alpha|
#sigs = linspace(0.10,0.30,21)
#zmin = zeros(len(sigs))
#Mmin = zeros(len(sigs))
#for i in range(0,len(sigs)):
#    isig  = sigs[i]
#    iP    = 2
#    ipargs = [0,isig,iP,1e-3]
#    npl   = plasma_ramp(np0,shape,z[0:len(z)-1],ipargs,'down')
#    beam = propbeam(beam0,z,npl[0:len(z)-1])
#    jM = zeros(len(beam))
#    jgamma = zeros(len(beam))
#    for j in range(0,len(beam)):
#        Ts[j][:] = [beam[j][2],beam[j][3],beam[j][4]]
#        jT = [beam[j][2],beam[j][3],beam[j][4]]
#        jM[j] = calc_M(jT,Ttarg)
#        jgamma[j] = np.abs(beam[j][4])
#    imin_z = np.argmin(jgamma[1:len(jgamma)-1])
#    zmin[i] = z[imin_z]
#    Mmin[i] = jM[min_z]
##    print(zmin[i])
##    print(Mmin[i])
#    
#plt.figure()
#plt.plot(zmin,Mmin,'.')


"""

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


## plot beam evolution
ymax_plt = 10*beam[len(beam)-1,2]
norm = ymax_plt/max(npl)
plt.figure()
plt.plot([z],[beam[:,2]],'b.',[z],[npl*norm],'g.')
plt.xlabel(r'z (m)')
plt.ylabel(r'$\beta$ (m)')
plt.title(r"Matched ramp with %s profile" % shape)
plt.text(0.5*max(z),0.5*ymax_plt,\
         r"$\sigma$ = %.2f m, waist = %.2f m" % (lpfit,waistfit),\
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

plt.figure()
plt.plot(xb,xpb,'b.',xt,xpt,'r.')
plt.title('Beam phase space ellipse')
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.text(0,0,"Bmag = %.3f" % Bmag,\
         verticalalignment='center',\
         horizontalalignment='center')
plt.show()

plt.figure()
plt.plot(xb,xpnb,'b.',xt,xpnt,'r.')
plt.title('Normalized beam phase space ellipse')
plt.xlabel(r'$\bar{x}$ (m)')
plt.ylabel(r'$\bar{x^{\prime}}$ (m)')
plt.text(0,0,"Bmag = %.3f" % Bmag,\
         verticalalignment='center',\
         horizontalalignment='center')
plt.show()

"""