#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:06:07 2017

@author: chris
"""
import sys
sys.path.insert(0, "../python")
import numpy as np
from propagation import propagation
from propagation import plasma
import os
import matplotlib.pyplot as plt

E0 = 26.90994
wx = 5.578897e-5*1e6
wy = 3.571468e-4*1e6
Px =1.559870e-9*1e6*1e6
Py =7.644225e-7*1e6*1e6
phi= 0.639129
zi = 0.005*1e6
n0 = 1e18*1e-17
z_off=1e-3*1e6
Lr=500e-6*1e6
Lz=500e-6*1e6

test_density = 0
free_space = 0
uniform_den = 0
gauss_den = 1

def Efunc(x,y):
    r2 = x**2 + y**2
    E = np.zeros(np.shape(r2))
    t1=np.exp(-(x**2)/(wx**2)-(y**2)/(wy**2))
    t2=np.exp(-1j*((x**2)/Px+(y**2)/Py))
    t3=np.exp(1j*phi)
    E=t1*t2*t3
    return E

# Temporal pulse
def Tfunc(t):
    tau = 50
    return np.exp(-1.38629*t**2 / (tau**2))

def Gaussian_density(x,y,z):
    nr=np.exp(-(np.power(z-zi,2)+np.power(x,2))/(2*np.power(Lr,2)))
    return n0*nr

def Gas_jet_density(x,y,z):
    nr=Gaussian_density(x,y,z)
    nz=np.exp(-(y+z_off)/Lz)
    return nr*nz
#9,9,7
# Setup the parameters
params = {'Nx' : 2**7,
          'Ny' : 2**7,
          'Nz' : 2**7,
          'Nt' : 2**6,
          'X' : 10.5e2/4,
          'Y' : 10.5e2,
          'Z' : 1e4,
          'T' : 150,
          'n0' : 5,
          'alpha' : 0.787,
          'EI' : 15.426,
          'E0' : E0,
          'lam' : 0.7853,
          'n' : 1.0,
          'angle' : 0.75
          }

def SetupArray(denfunc):
    nx=params['Nx']
    ny=params['Ny']
    nz=params['Nz']
    xstep=params['X']/nx
    ystep=params['Y']/ny
    zstep=params['Z']/nz
    n = np.zeros((nx,ny,nz))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                n[i][j][k]=denfunc(i*xstep,j*ystep,k*zstep)
    return n

path = '~/Desktop/FourierPlots/Users/doss/Documents/Research/Data/FACET_Refraction/'

if test_density == 1:
    h=SetupArray(Gaussian_density)
    hh=h[:,round(params['Ny']/2),:]
    plt.imshow(hh, interpolation="none", origin="lower",
               extent=[0,params['Z'],-.5*params['X'],.5*params['X']], aspect='50')
    CB=plt.colorbar()
    CB.set_label('n')
    plt.ylabel('x (microns)')
    plt.xlabel('z (microns)')
    plt.title('Gas Density; y=0')
    plt.show()

if free_space == 1:
    directory = 'free_space_propagation'
    params['path'] = path + directory+'/'
    # Create the directory if it doesn't exist
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
        # Simulate free space propagation
    propagation.laser_prop(params, Efunc)
    propagation.laser_prop_plot(params['path'])
    
if uniform_den == 1:
    directory = 'uniform_den_propagation_5e17'
    params['n0'] = 5
    params['path'] = path + directory+'/'
    # Create the directory if it doesn't exist
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
    # Run the simulation      
    plasma.plasma_refraction(params, Efunc, Tfunc)
    # Create the summary
    plasma.summary_plot(params['path'])
    
if gauss_den == 1:
    directory = 'guassian_den_propagation_1e18_gridtest'
    density = SetupArray(Gaussian_density)
    params['path'] = path + directory+'/'
    # Create the directory if it doesn't exist
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
    # Run the simulation      
    plasma.plasma_refraction(params, Efunc, Tfunc, density)
    # Create the summary
    plasma.summary_plot(params['path'])