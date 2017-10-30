#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:06:07 2017
Using Robert's code as an example, this runs realistic refractin
of a density distribution as the gaussian beam passes through.
Saves a density distribution in the given directory

A clean copy of FourierRefraction, using LoadDensity from interpolated
OpenFOAM results and Robert's new refraction code

@author: chris
"""
import sys
sys.path.insert(0, "../../python")
import numpy as np
from propagation import propagation
from propagation import plasma
import os
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

#Script Selection
test_density = 0
free_space = 0
#Set to some integer to lower grid resolution by a power
reducer = 0

#DONT FORGET TO CHECK THE DIRECTORY AND ITS STUFF!
# initDensity only knows its shape, not dimensions
path = '/home/chris/Desktop/FourierPlots/ArJets/'
directory = 'Ar1_Big'

#Parameters obtained using GaussianBeam.Prop_EPhase.  Consult OpticalSetup
#In order to work with Robert's code, there are a bunch of 1e6 factors to
# get everything into micrometers
pulseParams = np.load(path+directory+'/pulseParams.npy').item()
E0 = pulseParams['E0']
wx = pulseParams['wx']
wy = pulseParams['wy']
Px = pulseParams['Px']
Py = pulseParams['Py']
phi = pulseParams['phi']

#Our Electric field is a function of the stuff above
def Efunc(x,y,params):
    r2 = x**2 + y**2
    E = np.zeros(np.shape(r2))
    t1 = np.exp(-(x**2)/(wx**2) - (y**2)/(wy**2))
    t2 = np.exp(-1j * ((x**2)/Px + (y**2)/Py))
    t3 = np.exp(1j * phi)
    E = t1 * t2 * t3
    return E

# Temporal pulse
def Tfunc(t):
    tau = 50
    return np.exp(-1.38629*t**2 / (tau**2))

#Assuming everything in params matches the density array we are loading
#Loads output from CSVInterpolater from the same directory this script
# outputs to
def LoadDensity():
    return np.load(path+directory+'/initDensity.npy')

# Setup the parameters, just like in Robert's code for Ar
params = {'Nx' : 2**(9-reducer),
          'Ny' : 2**(9-reducer),
          'Nz' : 2**(8-reducer),
          'Nt' : 2**7,
          'X' : 3e2,
          'Y' : 15e2,
          'Z' : 1.5e4,
          'T' : 150,
          'n0' : 1,
          'alpha' : 1.63,
          'EI' : 15.8,
          'E0' : E0,
          'lam' : 0.7853,
          'n' : 1.0,
          'nFinal' : True
          }

X = params['X'];   Y = params['Y']
Nx = params['Nx']; Ny = params['Ny']
x = plasma.create_centered_grid(X, Nx)
y = plasma.create_centered_grid(Y, Ny)
E_initial = Efunc(np.reshape(x, (Nx, 1)), np.reshape(y, (1, Ny)), params)
E_initial *= params['E0']
I_initial=np.square(np.abs(E_initial))*.0013272
ThrDim.PlotInitialIntensity(I_initial,x,y)

diff = pulseParams['zi'] - params['Z']/2
if diff < 0:
    print("Near focus case: " + str(diff))
elif diff > 0:
    print("Far focus case: " + str(diff))
else:
    print("Focus aligned w/ gas center")

#If we are actually running a simulation
if test_density == 0:
    #Create the folder if it not already exists
    params['path'] = path + directory+'/'
    if not os.path.exists(params['path']):
        os.makedirs(params['path'])
    #Simulate free space
    if free_space == 1:
        propagation.laser_prop(params, Efunc)
        propagation.laser_prop_plot(params['path'])
    #Simulate a density distribution
    else:
        density = LoadDensity()
        plasma.plasma_refraction(params, Efunc, Tfunc, density)
        plasma.summary_plot(params['path'])
#Otherwise, we are making sure our density distribution works
else:
    h = LoadDensity()
    h = ThrDim.RobertRoll(h)
    X = params['X']; Nx = params['Nx']
    Y = params['Y']; Ny = params['Ny']
    Z = params['Z']; Nz = params['Nz']
    y = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(-Y/2, Y/2, Ny, False)
    x = np.linspace(-Z/2, Z/2, Nz, False)
    ThrDim.ImageCut(h,x,y,z,0,0,0,1e-3,'(mm)','Plasma Density','e17(cm^-3)',1)
    