#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:34:02 2017

@author: litos
"""

import numpy as np
from propbeam import propbeamlast
from propbeam import propbeam
from plasma_ramp import plasma_ramp
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from calc_M import calc_M
from plasma_lens import plasma_lens

# constants
c  = 3e8 # m/s
me = 0.511e-3 # GeV

# define plasma bulk (flat-top) properties
npl0   = 1e17 # cm^-3
dEds0  = 6.00 # GeV/m
dgds0  = dEds0/me
L_ft   = 0.50 # m

# define plasma up-ramp
shape_up = 'gauss'
L_up     = 0.50 # m

# define plasma lens
npl0_lens   = 1e17 # cm^-3
#s0_lens    = L_up - 0.10 # m
L_lens     = 100e-6 # m
dgds0_lens = dgds0*np.sqrt(npl0_lens/npl0)*(2*np.sqrt(npl0_lens/npl0)-1)

# define initial beam
gb0    = 20000 # relativistic lorentz factor
eps0   = 5.0e-6  # m-rad, normalized emittance

# use this for defining matched beam (if desired)
wp0 = (5.64e4)*np.sqrt(npl0) # rad/s, plasma ang. freq.
kp0 = wp0/c # m^-1, plasma wave number
kb0 = kp0/np.sqrt(2*gb0) # m^-1, betatron wave number

#beta0  = 1.0/kb0 # m
beta0  = 0.10 # m
alpha0 = 0.00
gamma0 = (1.0+alpha0**2)/beta0 # 1/m
dE0    = 0.01
npart  = 100
beam0  = [gb0,eps0,beta0,alpha0,gamma0,dE0,npart]

# scan parameters
nwaist = 51
waist  = np.linspace(-0.12,0.02,nwaist) # m, waist location w.r.t. L_up
#waist  = np.linspace(-0.98,-0.94,nwaist) # m, waist location w.r.t. L_up
#waist  = np.linspace(-0.48,-0.38,nwaist) # m, waist location w.r.t. L_up
nhw_up = 100
hw_up  = np.linspace(0.001,0.050,nhw_up) # m, HWHM of up-ramp
#hw_up  = np.linspace(0.30,0.32,nhw_up) # m, HWHM of up-ramp
#hw_up  = np.linspace(0.135,0.160,nhw_up) # m, HWHM of up-ramp
nlens  = 51
s0_lens = np.linspace(0.042,0.048,nlens)

#M = np.zeros([nwaist,nhw_up])
M = np.zeros([nwaist,nlens])
#M1 = np.zeros(nwaist)
for i in range(0,nwaist):
#    for j in range(0,nhw_up):
    for j in range(0,nlens):
            
        print(i*nhw_up+j+1)
        
        iwaist = waist[i]
#        jhw_up = hw_up[j]
        jhw_up = 0.03
        js0_lens = s0_lens[j]        

        # set beam waist
        s_w   = L_up + iwaist # m
        ibeam0  = propbeamlast(beam0,[0,-s_w])
        
        # make plasma up-ramp
        s_up     = np.linspace(0,L_up,round(L_up/0.0001))
        pargs_up = [L_up,jhw_up]
        [plas_up, dgds_up] =\
            plasma_ramp(npl0,shape_up,s_up[0:-1],pargs_up,'up',dgds0)
            
        # add plasma lens
        lens_args = [npl0_lens,L_up-js0_lens,L_lens]
        [plas_up, dgds_up] =\
            plasma_lens(plas_up,dgds_up,s_up[0:-1],lens_args,dgds0_lens)
            
        # simulate plasma up-ramp
        beam_up  = propbeamlast(ibeam0,s_up,plas_up,dgds_up)
        
#        # simulate plasma flat-top
#        s_ft     = np.linspace(0,L_ft,round(L_ft/0.001))
#        plas_ft  = npl0*np.ones(len(s_ft)-1)
#        dgds_ft  = dgds0*np.ones(len(s_ft)-1)
#        beam_ft  = propbeamlast(beam_up,s_ft,plas_ft,dgds_ft)
       
        Tbeam = beam_up[2:5]
        kb     = kp0/np.sqrt(2*beam_up[0])
        Tmatch = [1.0/kb,0,kb]
#        print(Tbeam)
#        print(Tmatch)
        M[i,j]   = calc_M(Tbeam,Tmatch)
#        M1[i] = calc_M(Tbeam,Tmatch)


#X = np.tile(waist.reshape(-1,1),(1,nhw_up))
#Y = np.tile(hw_up.T,(nwaist,1))

X = np.tile(waist.reshape(-1,1),(1,nlens))
Y = np.tile(s0_lens.T,(nwaist,1))

fig, axes = plt.subplots(1,1, sharey=True)
#fig = plt.figure()
#ax1 = fig.add_subplot(111)


plt.contourf(X,Y,np.log10(M),100,\
            cmap=cm.Vega20c,\
            linewidth=2.0)
#plt.pcolor(X,Y,-M)
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$log_{10}$(M)')

min_y = np.zeros(M.shape[0])
min_x = np.zeros(M.shape[0])
for i in range(0,M.shape[0]):
    min_y[i] = Y[0,i]
    min_x[i] = X[np.argmin(M[:,i]),0]

plt.scatter(min_x,min_y, s=10, c='b', marker=".", label='first')

plt.xlabel(r'$waist$')
#plt.ylabel(r'ramp half-length (m)')
plt.ylabel(r'plasma lens position (m)')
plt.title("Ramp Type: %s" % shape_up)
plt.show()
        
        
#fig = plt.figure()
#plt.scatter(waist,M[:,M.shape[1]-1])
        
fig = plt.figure()
plt.scatter(s_up[0:-1],plas_up)
        
        
        