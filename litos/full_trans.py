#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 17:20:47 2017

@author: litos
"""

import numpy as np
from mc_beam import mc_beam
from propbeam import propbeamlast
from propbeam import propbeam
from propparts import propparts
from plasma_ramp import plasma_ramp
import matplotlib.pyplot as plt
import time as time
import matplotlib.animation as animation
from calc_M import calc_M
from calc_rms_parts import calc_rms_parts
from calc_twiss import calc_twiss

# constants
c  = 3e8 # m/s
me = 0.511e-3 # GeV

# define plasma bulk (flat-top) properties
npl0 = 1e17 # cm^-3
dEds   = 4.00 # GeV/m
dgds   = dEds/me
L_ft   = 0.50 # m

# define plasma up-ramp
shape_up = 'gauss'
hw_up    = 1.0*0.11843264 # m
L_up     = 1.5 # m

# define plasma down-ramp
shape_dn = 'gauss'
hw_dn    = hw_up # m
L_dn     = L_up  # m

# define initial beam
gb0    = 20000 # relativistic lorentz factor
eps0   = 5e-6  # m-rad, normalized emittance

# use this for defining matched beam (if desired)
wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
kb0 = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number

#beta0  = 1.0/kb0 # m
beta0  = 0.10 # m
alpha0 = 0.00
gamma0 = (1.0+alpha0**2)/beta0 # 1/m
dE0    = 0.001
npart  = 500
beam0  = [gb0,eps0,beta0,alpha0,gamma0,dE0,npart]


#sig0 = np.sqrt(eps0*beta0/gb0)
#print(sig0)


# define imaging spectrometer (not used yet)
#E_img  = gamma0*me+dEds*L_ft # GeV
#D      = 62e-3 # m


# set beam waist
s_w   = L_up - 0.38390091 # m
beam0 = propbeamlast(beam0,[0,-s_w])
T0    = beam0[0:5]
#print(T0)

# generate particles
[x,xp,ux,vx,gb] = mc_beam(T0,npart,dE0)
[y,yp,uy,vy,gb] = mc_beam(T0,npart,dE0)
z = np.zeros(npart)
parts0 = np.zeros((npart,6))
for i in range(0,npart):
    parts0[i][:] = [x[i],xp[i],y[i],yp[i],z[i],gb[i]]


frac = 0.9
[rms_parts,rms_X,rms_Y] = calc_twiss(parts0,frac)
print(rms_X[0])
print(rms_Y[0])


# simulate plasma up-ramp
s_up     = np.linspace(0,L_up,round(L_up/0.001))
pargs_up = [L_up,hw_up,2.0,1e-3]
plas_up  = plasma_ramp(npl0,shape_up,s_up[0:len(s_up)-1],pargs_up)
#beam_up  = propbeam(beam0,s_up,plas_up)
parts_up = propparts(parts0,s_up,plas_up)

# simulate plasma flat-top
s_ft     = np.linspace(0,L_ft,round(L_ft/0.001))
plas_ft  = npl0*np.ones(len(s_ft)-1)
#beam_ft  = propbeam(beam_up[len(beam_up)-1][:],s_ft,plas_ft,dgds)
parts_ft = propparts(parts_up[len(parts_up)-1],s_ft,plas_ft,dgds)

# simulate plasma down-ramp
s_dn     = np.linspace(0,L_dn,round(L_dn/0.001))
pargs_dn = [L_dn,hw_dn,2.0,1e-3]
plas_dn  = plasma_ramp(npl0,shape_dn,s_dn[0:len(s_dn)-1],pargs_dn,'down')
#beam_dn  = propbeam(beam_ft[len(beam_ft)-1][:],s_dn,plas_dn)
parts_dn = propparts(parts_ft[len(parts_ft)-1],s_dn,plas_dn)

# simulate wakeless plasma undulator
L_un = 0.66
s_un     = np.linspace(0,L_un,round(L_un/0.001))
plas_un  = (5e16)*np.ones(len(s_un)-1)
#beam_ft  = propbeam(beam_ft[len(beam_ft)-1][:],s_un,plas_un)
#parts_un = propparts(parts_ft[len(parts_ft)-1],s_un,plas_un)
#parts_un = propparts(parts0,s_un,plas_un)

# combine results
#beam = np.vstack((beam_up,beam_ft,beam_dn))
parts = np.vstack((parts_up,parts_ft,parts_dn))
#parts = np.vstack((parts_ft,parts_dn))
#parts = parts_un


partsMid = parts_ft[len(parts_ft)-1]
[rms_parts,rms_X,rms_Y] = calc_twiss(partsMid,frac)
print(rms_X[0])
print(rms_Y[0])


partsEnd = parts_dn[len(parts_dn)-1]
[rms_parts,rms_X,rms_Y] = calc_twiss(partsEnd,frac)
print(rms_X[0])
print(rms_Y[0])


"""

# plot results
def update_plot(j, parts, scat, words):
    j = int(np.round(j))
    npart = parts.shape[1]
    xf  = np.zeros(npart)   
    xpf = np.zeros(npart)
    yf  = np.zeros(npart)
    ypf = np.zeros(npart)
    zf  = np.zeros(npart)
    gbf = np.zeros(npart)
    xxpf= np.zeros((npart,2))
    xyf = np.zeros((npart,2))
    for i in range(0,npart):
        [xf[i],xpf[i],yf[i],ypf[i],zf[i],gbf[i]] =\
            parts[j][i][:]
        xxpf[i][:] = [xf[i],xpf[i]]
        xyf[i][:] = [xf[i],yf[i]]
    scat.set_offsets(xyf)
    scat.set_array(gbf)
#    plt.xlim((-1.1*max(abs(xf)),1.1*max(abs(xf))))
##    plt.ylim((-1.1*max(abs(xpf)),1.1*max(abs(xpf))))
#    plt.ylim((-1.1*max(abs(yf)),1.1*max(abs(yf))))
    plt.clim((min(gbf),max(gbf)))
    words.set_text("%d" % j)
#    time.sleep(0.02)
    return scat, words

nstep = parts.shape[0]

fig  = plt.figure()
scat = plt.scatter([],[],marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
#plt.ylabel(r'$x^{\prime}$')
plt.ylabel(r'$y$ (m)')
words = plt.text(0,0,'')

wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
kb0 = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number
sigm = np.sqrt(eps0/(gb0*kb0))

#plt.xlim((-10*sigm,10*sigm))
#plt.ylim((-10*sigm,10*sigm))
plt.xlim((-30e-6,30e-6))
plt.ylim((-30e-6,30e-6))

ani = animation.FuncAnimation(fig, update_plot,\
                              frames=np.linspace(0,nstep-1,round(nstep/1)),\
                              fargs=(parts, scat, words),blit=True)
plt.show()

"""