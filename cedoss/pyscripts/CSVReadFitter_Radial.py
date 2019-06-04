#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:16:34 2017

Import data from a csv written by paraview for the purpose of fitting 
data (probably density) to a function (probably Gaussian)

From paraview, take a radial "plot over line" and then export to csv using
file -> save data.  Depending on the axis, may need to adjust the # below in
axis = 'Points:#' from 1 to 3.
    
If exporting a wedge simulation, to get a full 3D sim from a 2D wedge do
1. Transform:  rotate by 0,90,0
2. Slice along x normal
3. Rotational Extrusion

@author: chris
"""

import csv
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim
from modules import ProfileAnalyzer as Prof

directory = '/home/chris/Desktop/CSVFiles/'
filename = '3DLam_R3_DenRadial.csv'

path = directory + filename

variable = 'Den'
axis = 'Points:0'  #Typically 0 for a radial cut, 1 for an axial cut
simlabel = "Radial-Simulation"
offset = 0e-6 #For if you didnt set up your OpenFOAM to have nice axes
zoom = 1e6 #For when you have meters and work with micrometers

guess = [1.75e20,100,0] #Tweak this if python doesnt want to fit
guess_tanh = [10,100,1e19]
guess_super = [1.75e20,10,0,2]

arr = []
dist = []

with open(path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        arr.append(float(row[variable]))
        dist.append(float(row[axis])-offset)
        
    arr = np.array(arr)
    dist = np.array(dist)
    
    arr = arr[round(len(arr)*1/4):round(len(arr)*3/4)]
    dist = dist[round(len(dist)*1/4):round(len(dist)*3/4)]
    
    dist=dist*1e6
#    arr=arr-arr[0]
    
    xwindow = dist[-1]; xstep = (dist[2]-dist[1])/50 #microns
    dist_precise = np.arange(-xwindow, xwindow, xstep)
    
    gfit = ThrDim.FitDataGaussian(arr, dist, guess)
    #tfit = ThrDim.FitDataDoubleTanh(arr2, dist2, guess_tanh)
    #tfit = ThrDim.FitDataDoubleTanhAbs(arr, dist, guess_tanh)
    lfit = ThrDim.FitDataLorentz(arr, dist, guess)
    #sfit = ThrDim.FitDataSuperGaussian(arr2, dist2, guess_super)
    
    plt.plot(dist,arr,label=simlabel)
    plt.title("Lorentzian Fit")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.plot(dist, ThrDim.Lorentz(lfit,dist),label="Lorentzian")
    plt.grid(); plt.legend(); plt.show()
    """
    plt.plot(dist,arr,label=simlabel)
    plt.title("Gaussian Fit")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.plot(dist, ThrDim.Gaussian(gfit,dist),label="Gaussian")
    plt.grid(); plt.legend(); plt.show() 
    
    plt.plot(dist,arr,label=simlabel)
    plt.title("SuperGaussian Fit")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.plot(dist, ThrDim.SuperGaussian(sfit,dist),label="SuperGaussian: p="+str(sfit[3]))
    plt.grid(); plt.legend(loc=4); plt.show()

    plt.plot(dist,arr,label=simlabel)
    plt.title("Double Tanh Fit")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.plot(dist, ThrDim.DoubleTanh(tfit,dist),label="Double Tanh")
    plt.plot(dist, ((.5 + .5*np.tanh((dist+tfit[0])/tfit[1]))) * tfit[2], linestyle='dashed', label="Tanh 1")
    plt.plot(dist, ((.5 - .5*np.tanh((dist-tfit[0])/tfit[1]))) * tfit[2], linestyle='dashed', label="Tanh 2")
    plt.grid(); plt.legend(); plt.show()
    """    
    
    ###  Calculates FWHM and flattop length
    fwhm = Prof.CalcFWHM(dist_precise,ThrDim.Gaussian(gfit,dist_precise))
    flattop = Prof.CalcTopWidth(dist_precise,ThrDim.Gaussian(gfit,dist_precise),.99)
    plt.plot(dist,arr, label = 'Simulation Profile')
    plt.plot(dist_precise, ThrDim.Gaussian(gfit,dist_precise),label="Gaussian")
    plt.title("Radial Distribution")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.axvspan(-fwhm/2, -flattop/2, facecolor = 'b', alpha = 0.3, label='FWHM: '+str(fwhm))
    plt.axvspan(-flattop/2, flattop/2, facecolor = 'g', alpha = 0.3, label='99% Len: '+str(flattop))
    plt.axvspan(flattop/2, fwhm/2, facecolor = 'b', alpha = 0.3)
    plt.xlim([-5000,5000])
    plt.grid(); plt.legend(); plt.show()
    
    #plt.plot(dist, ThrDim.SuperGaussian(sfit,dist) - ThrDim.Gaussian(gfit,dist))
    print('sig',gfit[1],'um')