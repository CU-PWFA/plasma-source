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
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

#Script Selection
test_density = 1
free_space = 0
#Set to some integer to lower grid resolution by a power
reducer = 0

#This is the directory for saving
#path = '/home/chris/Desktop/FourierPlots/CompactOptics_DoubleJet/'
path = '/home/chris/Desktop/FourierPlots/CompactOptics_JetsInGas/'
directory = 'Exact_p2g8_Ar'
#directory = 'testdir'

#Density, converted to e17cm^-3  With a gas jet this is density at nozzle
#n0 = 1e19 * 1e-17

#Density of background gas, if using gas jets in a gas.  Set to 0 otherwise
n_b = 0

#Parameters obtained using GaussianBeam.Prop_EPhase.  Consult OpticalSetup
#In order to work with Robert's code, there are a bunch of 1e6 factors to
# get everything into micrometers
E0 = 33.3376399047
wx = 7.64942547328e-05 * 1e6
wy = 0.000339431499556 * 1e6
Px = 1.3083925296801255e-09 * 1e6 * 1e6
Py = 5.7480222100767626e-08 * 1e6 * 1e6
phi = 4.09336832829
zi = 0.005*1e6

#Gas Jet Parameters for ideal jets, converted to micrometers
z_off=1e-3 * 1e6
Lr=500e-6 * 1e6
Lz=800e-6 * 1e6

"""#This is for p1g8
#Tanh Parameters (in microns)
a = 385.588869202
b = 89.8895905719
n_0 = 1.32664984545e+19 * 1e-17
p_tanh=[a,b,n_0]
#Gauss Parameters, though all we need is sigma
sig = 2099.27300216
"""
#This is for p2g8
#Tanh Parameters (in microns)
a = 362.080596888
b = 55.9371159276
n_0 = 1.31332220886e+18 * 1e-17 / 4.
p_tanh=[a,b,n_0]
#Gauss Parameters, though all we need is sigma
sig = 2066.08516439

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

#DensityDistributions
def Gaussian_density(x,y,z):
    nr=np.exp(-(np.power(z-zi,2)+np.power(x,2))/(2*np.power(Lr,2)))
    return n_0*nr

def Gas_jet_density(x,y,z):
    nr=Gaussian_density(x,y,z)
    nz=np.exp(-(y+z_off)/Lz)
    return nr*nz

def Double_Jet(x,y,z):
    j1=Gas_jet_density(x,y,z)
    j2=Gas_jet_density(x,-y,z)
    return j1+j2

def Double_Jet_In_Gas(x,y,z):
    params['n0'] = n_b
    return Double_Jet(x,y,z) + n_b

def Tanh_Jet(x,y,z):
    nr = np.exp(-(np.power(z-zi,2)+np.power(x,2))/(2*np.power(sig,2)))
    ny = ((.5 + .5*np.tanh((y+p_tanh[0])/p_tanh[1])) * 
            (.5 - .5*np.tanh((y-p_tanh[0])/p_tanh[1])))
    return n_0*nr*ny

#Assuming everything in params matches the density array we are loading
#Loads output from CSVInterpolater from the same directory this script
# outputs to
def LoadDensity():
    return np.load(path+directory+'/initDensity.npy')

#False for uniform density, otherwise pick your poison
distribution = LoadDensity

# Setup the parameters, just like in Robert's code for Ar
params = {'Nx' : 2**(9-reducer),
          'Ny' : 2**(9-reducer),
          'Nz' : 2**(8-reducer),
          'Nt' : 2**6,
          'X' : 10.5e2/4,
          'Y' : 10.5e2,
          'Z' : 1e4,
          'T' : 150,
          'n0' : 0.1,
          'alpha' : 1.63,
          'EI' : 15.8,
          'E0' : E0,
          'lam' : 0.7853,
          'n' : 1.0,
          'angle' : 0.75,
          'nFinal' : True
          }

#Create the density distribution as n(x,y,z)
def SetupArray(denfunc):
    if denfunc is False:
        params['n0'] = n_0
        return False
    if (denfunc.__name__ == 'LoadDensity'):
        return LoadDensity()
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
        density = SetupArray(distribution)
        plasma.plasma_refraction(params, Efunc, Tfunc, density)
        plasma.summary_plot(params['path'])
#Otherwise, we are making sure our density distribution works
else:
    h=SetupArray(distribution)
    """
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
    """
    h = ThrDim.RobertRoll(h)
    X = params['X']; Nx = params['Nx']
    Y = params['Y']; Ny = params['Ny']
    Z = params['Z']; Nz = params['Nz']
    y = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(-Y/2, Y/2, Ny, False)
    x = np.linspace(-Z/2, Z/2, Nz, False)
    ThrDim.ImageCut(h,x,y,z,0,0,0,1e-3,'(mm)','Plasma Density','e17(cm^-3)',1)
    