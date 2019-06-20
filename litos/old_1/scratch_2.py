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

# constants
c  = 3e8 # m/s
me = 0.511e-3 # GeV

# define plasma bulk (flat-top) properties
npl0   = 1e18 # cm^-3
dEds0  = 6.00 # GeV/m
dgds0  = dEds0/me
L_ft   = 0.50 # m

# define plasma up-ramp
shape_up = 'gauss'
L_up     = 1.00 # m

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
waist  = np.linspace(-1.2,0,nwaist) # m, waist location w.r.t. L_up
#waist  = np.linspace(-0.98,-0.94,nwaist) # m, waist location w.r.t. L_up
#waist  = np.linspace(-0.48,-0.38,nwaist) # m, waist location w.r.t. L_up
nhw_up = 51
hw_up  = np.linspace(0.01,0.45,nhw_up) # m, HWHM of up-ramp
#hw_up  = np.linspace(0.30,0.32,nhw_up) # m, HWHM of up-ramp
#hw_up  = np.linspace(0.135,0.160,nhw_up) # m, HWHM of up-ramp

M = np.zeros([nwaist,nhw_up])
#M1 = np.zeros(nwaist)
for i in range(0,nwaist):
    for j in range(0,nhw_up):
        
        print(i*nhw_up+j+1)
        
        iwaist = waist[i]
        jhw_up = hw_up[j]
#        jhw_up = 0.05
        
        # set beam waist
        s_w   = L_up + iwaist # m
        ibeam0  = propbeamlast(beam0,[0,-s_w])
        
        # simulate plasma up-ramp
        s_up     = np.linspace(0,L_up,round(L_up/0.001))
        pargs_up = [L_up,jhw_up]
        [plas_up, dgds_up] =\
            plasma_ramp(npl0,shape_up,s_up[0:len(s_up)-1],pargs_up,'up',dgds0)
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


X = np.tile(waist.reshape(-1,1),(1,nhw_up))
Y = np.tile(hw_up.T,(nwaist,1))
fig = plt.figure()
plt.contour(X,Y,np.log10(M),100,\
            cmap=cm.Vega20c,\
            linewidth=2.0)
#plt.pcolor(X,Y,-M)
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$log_{10}$(M)')
plt.xlabel(r'$waist$')
plt.ylabel(r'ramp half-length (m)')
plt.title("Ramp Type: %s" % shape_up)
plt.show()
        
        
#fig = plt.figure()
#plt.scatter(waist,M[:,M.shape[1]-1])
        
        
        
        
        