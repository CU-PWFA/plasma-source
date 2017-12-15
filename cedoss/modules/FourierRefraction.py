#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:06:07 2017

Functions for running Robert's code.  Works with FourierRefractionLooper.py
which provides a params dictionary (includes variables for asym Gaussian pulse)
and an optional 3D density array for nonuniform gas propagation.

Other functions include appending the dictionary with pulse parameters,
plotting the initial intensity, and plotting initial density

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

def GetDefaultParams(reducer = 0, double_x = 0):
    def_params = {'Nx' : 2**(9-reducer)  * (1+double_x),
               'Ny' : 2**(9-reducer),
               'Nz' : 2**(8-reducer),
               'Nt' : 2**7,#7
               'X' : 3e2 * (1+double_x),
               'Y' : 15e2,
               'Z' : 1.5e4,
               'T' : 160,
               'n0' : 0.1,
               'alpha' : 1.63,
               'EI' : 15.8,
               'E0' : 1,
               'lam' : 0.7853,
               'n' : 1.0,
               'nFinal' : True,
              }
    return def_params

def RunRefraction(density,params):
    #Sanity check so that we are using correct pulse
    diff = params['d_zi'] - params['Z']/2
    if diff < 0:
        print("Near focus case: " + str(diff))
    elif diff > 0:
        print("Far focus case: " + str(diff))
    else:
        print("Focus aligned w/ gas center")
        
    if not os.path.exists(params['path']):
        print("Creating new directory")
        os.makedirs(params['path'])
    else: print("Directory already exists")
    
    plasma.plasma_refraction(params, Efunc, Tfunc, density)
    plasma.summary_plot(params['path'])

#Our Electric field is a function of the stuff above
def Efunc(x,y,params):
    r2 = x**2 + y**2
    E = np.zeros(np.shape(r2))
    t1 = np.exp(-(x**2)/(params['d_wx']**2) - (y**2)/(params['d_wy']**2))
    t2 = np.exp(-1j * ((x**2)/params['d_Px'] + (y**2)/params['d_Py']))
    t3 = np.exp(1j * params['d_phi'])
    E = t1 * t2 * t3
    return E

# Temporal pulse
def Tfunc(t):
    tau = 35
    return np.exp(-1.38629*t**2 / (tau**2))

def PlotInitialIntensity(params):
    X = params['X'];   Y = params['Y']
    Nx = params['Nx']; Ny = params['Ny']
    x = plasma.create_centered_grid(X, Nx)
    y = plasma.create_centered_grid(Y, Ny)
    E_initial = Efunc(np.reshape(x, (Nx, 1)), np.reshape(y, (1, Ny)), params)
    E_initial *= params['E0']
    I_initial=np.square(np.abs(E_initial))*.0013272
    ThrDim.PlotInitialIntensity(I_initial,x,y)

def LoadPulseParams(pulsefn,params):
    pulseParams = np.load(pulsefn).item()
    params['E0'] = pulseParams['E0']
    params['d_wx'] = pulseParams['wx']
    params['d_wy'] = pulseParams['wy']
    params['d_Px'] = pulseParams['Px']
    params['d_Py'] = pulseParams['Py']
    params['d_phi'] = pulseParams['phi']
    params['d_zi'] = pulseParams['zi']
    return params

def RunFreeSpace(params):
    if not os.path.exists(params['path']):
        print("Creating new directory")
        os.makedirs(params['path'])
    else: print("Directory already exists")
    propagation.laser_prop(params, Efunc)
    propagation.laser_prop_plot(params['path'])

def TestDensity(density,params):
    if density is False:
        print("Constant density of " + str(params['n0']))
    else:
        density = ThrDim.RobertRoll(density)
        X = params['X']; Nx = params['Nx']
        Y = params['Y']; Ny = params['Ny']
        Z = params['Z']; Nz = params['Nz']
        y = np.linspace(-X/2, X/2, Nx, False)
        z = np.linspace(-Y/2, Y/2, Ny, False)
        x = np.linspace(-Z/2, Z/2, Nz, False)
        ThrDim.ImageCut(density,x,y,z,0,0,0,1e-3,'(mm)','Plasma Density','e17(cm^-3)',1)