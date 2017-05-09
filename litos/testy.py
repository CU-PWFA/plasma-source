#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 17:07:46 2017

@author: litos
"""

import numpy as np
import matplotlib.pyplot as plt
import time as time
import matplotlib.animation as animation
from mc_beam import mc_beam
from propparts import propparts

gb0   = 20000
eps   = 5e-6
beta  = 0.10
alpha = +1.2
gamma = (1+alpha**2)/beta
dE    = 0.01

Tbeam = [gb0,eps,beta,alpha,gamma]

npart = 1000

np0 = 1e17


[x,xp,ux,vx,gb] = mc_beam(Tbeam,npart,dE)
[y,yp,uy,vy,gb] = mc_beam(Tbeam,npart,dE)

# plot beam particles
plt.figure()
plt.scatter(x,xp,c=gb,marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(x),1.1*max(x)))
plt.ylim((1.1*min(xp),1.1*max(xp)))
plt.show()

plt.figure()
plt.scatter(ux,vx,c=gb,marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(ux),1.1*max(ux)))
plt.ylim((1.1*min(vx),1.1*max(vx)))
plt.show()

z = np.zeros(npart)
parts0 = np.zeros((npart,6))
for i in range(0,npart):
    parts0[i][:] = [x[i],xp[i],y[i],yp[i],z[i],gb[i]]
    
s = np.linspace(0,0.2,1000)

npl = np0*np.ones(len(s)-1)

parts = propparts(parts0,s,npl)

xf  = np.zeros(npart)   
xpf = np.zeros(npart)
yf  = np.zeros(npart)
ypf = np.zeros(npart)
zf  = np.zeros(npart)
gbf = np.zeros(npart)
for i in range(0,npart):
    [xf[i],xpf[i],yf[i],ypf[i],zf[i],gbf[i]] =\
        parts[parts.shape[0]-1][i][:]
        
plt.figure()
plt.scatter(xf,xpf,c=gbf,marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(xf),1.1*max(xf)))
plt.ylim((1.1*min(xpf),1.1*max(xpf)))
plt.show()



xf  = np.zeros(npart)   
xpf = np.zeros(npart)
yf  = np.zeros(npart)
ypf = np.zeros(npart)
zf  = np.zeros(npart)
gbf = np.zeros(npart)
for i in range(0,npart):
    [xf[i],xpf[i],yf[i],ypf[i],zf[i],gbf[i]] =\
        parts[0][i][:]
        
plt.figure()
plt.scatter(xf,xpf,c=gbf,marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(xf),1.1*max(xf)))
plt.ylim((1.1*min(xpf),1.1*max(xpf)))
plt.show()




def update_plot(j, parts, scat):
    j = int(np.round(j))
    npart = parts.shape[1]
    xf  = np.zeros(npart)   
    xpf = np.zeros(npart)
    yf  = np.zeros(npart)
    ypf = np.zeros(npart)
    zf  = np.zeros(npart)
    gbf = np.zeros(npart)
    xxpf= np.zeros((npart,2))
    for i in range(0,npart):
        [xf[i],xpf[i],yf[i],ypf[i],zf[i],gbf[i]] =\
            parts[j][i][:]
        xxpf[i][:] = [xf[i],xpf[i]]
    scat.set_offsets(xxpf)
    scat.set_array(gbf)
    time.sleep(0.02)
    return scat,

nstep = parts.shape[0]

fig  = plt.figure()

#scat = plt.scatter(xf,xpf,c=gbf,marker='.')
scat = plt.scatter([],[],marker='.')
plt.jet()
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$x^{\prime}$')
plt.xlim((1.1*min(xf),1.1*max(xf)))
plt.ylim((1.1*min(xpf),1.1*max(xpf)))

ani = animation.FuncAnimation(fig, update_plot,\
                              frames=np.linspace(0,nstep-1,nstep),\
                              fargs=(parts, scat),blit=True)
plt.show()








#
#fig = plt.figure()
#ax  = fig.add_subplot(111)
#
#line1, = ax.scatter(x,xpf,c=gbf,marker='.') # Returns a tuple of line objects, thus the comma
#
#plt.axis([1.1*min(xf),1.1*max(xf),1.1*min(xpf),1.1*max(xpf)])
#plt.xlabel(r'$x$ (m)')
#plt.ylabel(r'$x^{\prime}$')
#plt.jet()
#plt.ion()
#nstep = parts.shape[0]
#for j in range(0,nstep):
#    xf  = np.zeros(npart)   
#    xpf = np.zeros(npart)
#    yf  = np.zeros(npart)
#    ypf = np.zeros(npart)
#    zf  = np.zeros(npart)
#    gbf = np.zeros(npart)
#    for i in range(0,npart):
#        [xf[i],xpf[i],yf[i],ypf[i],zf[i],gbf[i]] =\
#            parts[j][i][:]
#
#    line1.set_xdata(xf)
#    line1.set_ydata(xpf)
#    line1.set_cdata(gbf)
#    fig.canvas.draw()
#
#
##    plt.cla()
##    plt.scatter(xf,xpf,c=gbf,marker='.')
###    plt.jet()
###    plt.xlabel(r'$x$ (m)')
###    plt.ylabel(r'$x^{\prime}$')
###    plt.xlim((1.1*min(xf),1.1*max(xf)))
###    plt.ylim((1.1*min(xpf),1.1*max(xpf)))
##    plt.show()
#
#    plt.pause(1)