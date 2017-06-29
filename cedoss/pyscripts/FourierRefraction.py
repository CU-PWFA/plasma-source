#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:06:07 2017
Using Robert's code as an example, this runs realistic refractin
of a density distribution as the gaussian beam passes through.
Saves a density distribution in the given directory

@author: chris
"""
import sys
sys.path.insert(0, "../../python")
import numpy as np
from propagation import propagation
from propagation import plasma
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#Script Selection
test_density = 1
free_space = 0
uniform_den = 0
gauss_den = 0
gasjet_den = 0

#This is the directory for saving
path = '/home/chris/Desktop/FourierPlots/real_FACET_Refraction/'
directory = 'gasjet_den_propagation_1e17'

#Parameters obtained using GaussianBeam.Prop_EPhase
E0 = 38.056404817
wx = 5.578897e-5*1e6
wy = 3.571468e-4*1e6
Px =1.559870e-9*1e6*1e6
Py =7.644225e-7*1e6*1e6
phi= 0.639129
zi = 0.005*1e6

#Parameters for the gas jet
n0 = 1e17*1e-17
z_off=1e-3*1e6
Lr=500e-6*1e6
Lz=500e-6*1e6

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

#DensityDistributions
def Gaussian_density(x,y,z):
    nr=np.exp(-(np.power(z-zi,2)+np.power(x,2))/(2*np.power(Lr,2)))
    return n0*nr

def Gas_jet_density(x,y,z):
    nr=Gaussian_density(x,y,z)
    nz=np.exp(-(y+z_off)/Lz)
    return nr*nz

# Setup the parameters, just like in Robert's code
params = {'Nx' : 2**9,
          'Ny' : 2**9,
          'Nz' : 2**7,
          'Nt' : 2**6,
          'X' : 10.5e2/4,
          'Y' : 10.5e2,
          'Z' : 1e4,
          'T' : 150,
          'n0' : 0,
          'alpha' : 0.787,
          'EI' : 15.426,
          'E0' : E0,
          'lam' : 0.7853,
          'n' : 1.0,
          'angle' : 0.75,
          'nFinal' : True
          }

#Create the density distribution as n(x,y,z)
def SetupArray(denfunc):
    nx=params['Nx']
    ny=params['Ny']
    nz=params['Nz']
    xstep=params['X']/nx
    ystep=params['Y']/ny
    zstep=params['Z']/nz
    n = np.zeros((nx,ny,nz))
    xstart=params['X']/2
    ystart=params['Y']/2
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                n[i][j][k]=denfunc((i*xstep)-xstart,(j*ystep)-ystart,k*zstep)
    return n

#Plots density distribution
if test_density == 1:
    h=SetupArray(Gas_jet_density)
    hh=h[:,round(params['Ny']/2),:]
    
    gridSize = (2,5)
    plt.figure(figsize=(16,9))
    gridspec.GridSpec(gridSize[0],gridSize[1])
    plt.subplot2grid(gridSize, (0,0), colspan=4)
    plt.imshow(hh,
               aspect='auto',
               extent=[0,params['Z'],-.5*params['X'],.5*params['X']])
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
    directory = 'uniform_den_propagation_1e17_new'
    params['n0'] = 1
    params['path'] = path + directory+'/'
    # Create the directory if it doesn't exist
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
    # Run the simulation      
    plasma.plasma_refraction(params, Efunc, Tfunc)
    # Create the summary
    plasma.summary_plot(params['path'])
    
if gauss_den == 1:
    directory = 'guassian_den_propagation_1e20'
    density = SetupArray(Gaussian_density)
    params['path'] = path + directory+'/'
    # Create the directory if it doesn't exist
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
    # Run the simulation      
    plasma.plasma_refraction(params, Efunc, Tfunc, density)
    # Create the summary
    plasma.summary_plot(params['path'])
    
if gasjet_den == 1:
    density = SetupArray(Gas_jet_density)
    params['path'] = path + directory+'/'
    # Create the directory if it doesn't exist
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
    # Run the simulation      
    plasma.plasma_refraction(params, Efunc, Tfunc, density)
    # Create the summary
    plasma.summary_plot(params['path'])