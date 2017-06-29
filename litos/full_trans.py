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
from propparts import proppartslast
from propparts import propparts
from plasma_ramp import plasma_ramp
import matplotlib.pyplot as plt
import time as time
import matplotlib.animation as animation
from calc_M import calc_M
from calc_rms import calc_rms
from calc_twiss import calc_twiss

# constants
c  = 3e8 # m/s
me = 0.511e-3 # GeV

# define plasma bulk (flat-top) properties
npl0   = 1e17 # cm^-3
dEds0  = -20.00 # GeV/m
dgds0  = dEds0/me
L_ft   = 0.43 # m

# define plasma up-ramp
shape_up = 'gauss'
hw_up    = 0.030 #0.137 #0.311 # m
L_up     = 0.30 # m

# define plasma down-ramp
shape_dn = 'gauss'
hw_dn    = hw_up # m
L_dn     = L_up  # m

# define plasma lens
npl0_lens   = 1e17 # cm^-3
s0_lens    = L_up - 0.044 # m
L_lens     = 100e-6 # m
dgds0_lens = dgds0*np.sqrt(npl0_lens/npl0)*(2*np.sqrt(npl0_lens/npl0)-1)

# define initial beam
gb0    = 20000 # relativistic lorentz factor
eps0   = 5e-6  # m-rad, normalized emittance

# use this for defining matched beam (if desired)
wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
kb0 = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number

beta0  = 1.0/kb0 # m
#beta0  = 0.10 # m
alpha0 = 0.00
gamma0 = (1.0+alpha0**2)/beta0 # 1/m
dE0    = 0.01
npart  = 1000
beam0  = [gb0,eps0,beta0,alpha0,gamma0,dE0,npart]

# define imaging spectrometer (not used yet)
#E_img  = gamma0*me+dEds*L_ft # GeV
#D      = 62e-3 # m

# generate particles
[x,xp,ux,vx,gb] = mc_beam(beam0)
[y,yp,uy,vy,gb] = mc_beam(beam0)
z = np.zeros(npart)
parts0 = np.zeros((npart,6))
for i in range(0,npart):
    parts0[i][:] = [x[i],xp[i],y[i],yp[i],z[i],gb[i]]

## set beam waist
#s_w   = L_up -0.067 #-0.350 #-0.963 # m
#parts0 = proppartslast(parts0,[0,-s_w])
#beam0  = propbeamlast(beam0,[0,-s_w])
#
## simulate plasma up-ramp
#s_up     = np.linspace(0,L_up,round(L_up/0.0001))
#pargs_up = [L_up,hw_up]
#[plas_up, dgds_up] = plasma_ramp(npl0,shape_up,s_up[0:len(s_up)-1],pargs_up,'up',dgds0)
#beam_up  = propbeam(beam0,s_up,plas_up,dgds_up)
#parts_up = propparts(parts0,s_up,plas_up,dgds_up)

# simulate plasma flat-top
s_ft     = np.linspace(0,L_ft,round(L_ft/0.001))
plas_ft  = npl0*np.ones(len(s_ft)-1)
dgds_ft  = dgds0*np.ones(len(s_ft)-1)
#beam_ft  = propbeam(beam_up[len(beam_up)-1],s_ft,plas_ft,dgds_ft)
#parts_ft = propparts(parts_up[len(parts_up)-1],s_ft,plas_ft,dgds_ft)
beam_ft  = propbeam(beam0,s_ft,plas_ft,dgds_ft)
parts_ft = propparts(parts0,s_ft,plas_ft,dgds_ft)

## simulate plasma down-ramp
#s_dn     = np.linspace(0,L_dn,round(L_dn/0.001))
#pargs_dn = [L_dn,hw_dn]
#[plas_dn, dgds_dn] = plasma_ramp(npl0,shape_dn,s_dn[0:len(s_dn)-1],pargs_dn,'down',dgds0)
#beam_dn  = propbeam(beam_ft[len(beam_ft)-1],s_dn,plas_dn)
#parts_dn = propparts(parts_ft[len(parts_ft)-1],s_dn,plas_dn)



#beam = np.vstack((beam_up,beam_ft))
#parts = np.vstack((parts_up,parts_ft))
beam = beam_ft
parts = parts_ft

sig_p  = np.zeros(len(parts))
eps_p  = np.zeros(len(parts))
beta_p = np.zeros(len(parts))
for i in range(0,len(parts)):
    [rms_parts,rms_X,rms_Y] = calc_twiss(parts[i])
    sig_p[i]   = rms_parts[0]
    eps_p[i]   = rms_X[0]
    beta_p[i]  = rms_X[1]

gb_b    = np.zeros(len(beam))
eps_b   = np.zeros(len(beam))
beta_b  = np.zeros(len(beam))
alpha_b = np.zeros(len(beam))
gamma_b = np.zeros(len(beam))
M       = np.zeros(len(beam))
for i in range(0,len(beam)):
    [gb_b[i],eps_b[i],beta_b[i],alpha_b[i],gamma_b[i],dE,npart] = beam[i]
    Tbeam  = [beta_b[i],alpha_b[i],gamma_b[i]]
    kb     = kp0/np.sqrt(2*gb_b[i])
    Tmatch = [1.0/kb,0,kb]
#    print(Tbeam)
#    print(Tmatch)
    M[i]   = calc_M(Tbeam,Tmatch)
#    print(M[i])

#plot_x = np.hstack((s_up,s_up[len(s_up)-1]+s_ft))
#plot_x = np.hstack((plot_x,s_up[len(s_up)-1]+s_ft[len(s_ft)-1]+s_dn))
plot_x = s_ft
plot_y1 = eps_p/eps_p[0]
plot_y2 = M
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(plot_x, plot_y1, s=10, c='b', marker=".", label='first')
#ax1.scatter(plot_x, plot_y2, s=10, c='r', marker=".", label='second')
plt.legend(loc='upper left');
plt.show()
plt.xlabel(r'$s$ (m)')
plt.ylabel(r'$\epsilon/\epsilon_0$')
#plt.ylim([0.9*min(plot_y1),1.1*max(plot_y1)])


x      = np.zeros(len(parts[0]))
xp     = np.zeros(len(parts[0]))
y      = np.zeros(len(parts[0]))
yp     = np.zeros(len(parts[0]))
z      = np.zeros(len(parts[0]))
gb     = np.zeros(len(parts[0]))

for j in range(0,len(parts[i])):
        [x[j],xp[j],y[j],yp[j],z[j],gb[j]] = parts[len(parts)-1,j,:]
fig = plt.figure()
plt.hist(x,25)
fig = plt.figure()
plt.hist(y,25)

print(M[len(M)-1])
print(np.sqrt(M[len(M)-1]))
print(eps_p[len(eps_p)-1]/eps_p[0])

#plot_x = np.hstack((s_up,s_up[len(s_up)-1]+s_ft))
#plot_y = beam[:,2]
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#
#ax1.scatter(plot_x, plot_y, s=10, c='b', marker=".", label='first')
#plt.legend(loc='upper left');
#plt.show()

"""





# simulate plasma down-ramp
s_dn     = np.linspace(0,L_dn,round(L_dn/0.001))
pargs_dn = [L_dn,hw_dn]
plas_dn  = plasma_ramp(npl0,shape_dn,s_dn[0:len(s_dn)-1],pargs_dn,'down')
beam_dn  = propbeam(beam_ft[len(beam_ft)-1],s_dn,plas_dn)
parts_dn = propparts(parts_ft[len(parts_ft)-1],s_dn,plas_dn)

# simulate wakeless plasma undulator
#L_un = 0.66
#s_un     = np.linspace(0,L_un,round(L_un/0.001))
#plas_un  = (5e16)*np.ones(len(s_un)-1)
#beam_ft  = propbeam(beam_ft[len(beam_ft)-1][:],s_un,plas_un)
#parts_un = propparts(parts_ft[len(parts_ft)-1],s_un,plas_un)
#parts_un = propparts(parts0,s_un,plas_un)

# combine results
#beam = np.vstack((beam_up,beam_ft,beam_dn))
beam = np.vstack((beam_ft,beam_dn))
#parts = np.vstack((parts_up,parts_ft,parts_dn))
parts = np.vstack((parts_ft,parts_dn))
#parts = parts_un

print('start:')
partsStart = parts_up[len(parts_ft)-1]
[rms_parts,rms_X,rms_Y] = calc_twiss(partsStart,frac)
print(rms_X[0])
print(rms_Y[0])
print(rms_X)
print(np.mean(partsStart[:,5]))

beamStart = beam_up[len(beam_up)-1]
Tstart = beamStart[2:5]
beta_m = np.sqrt(2*Tstart[0])/kp0
Tflat  = [beta_m,0,1.0/beta_m]
Mstart = calc_M(Tstart,Tflat)
print(Tstart)
print(Tflat)
print(Mstart)
print(beamStart[0])

print('mid:')
partsMid = parts_ft[len(parts_ft)-1]
[rms_parts,rms_X,rms_Y] = calc_twiss(partsMid,frac)
print(rms_X[0])
print(rms_Y[0])
print(np.mean(partsMid[:,5]))

beamMid = beam_ft[len(beam_ft)-1]
Tmid = beamMid[2:5]
beta_m = np.sqrt(2*Tmid[0])/kp0
Tflat  = [beta_m,0,1.0/beta_m]
Mmid = calc_M(Tmid,Tflat)
print(Mmid)
print(beamMid[0])

print('end:')
partsEnd = parts_dn[len(parts_dn)-1]
[rms_parts,rms_X,rms_Y] = calc_twiss(partsEnd,frac)
print(rms_X[0])
print(rms_Y[0])
print(np.mean(partsEnd[:,5]))
print(calc_rms(partsEnd[:,5]))

"""


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