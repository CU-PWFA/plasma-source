#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:16:34 2017

Import data from a cvs written by paraview for the purpose of fitting 
data (probably density) to a function (probably Gaussian)

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
filename = '3DLES_Wide_Axial.csv'
path = directory + filename

variable = 'Density'
axis = 'Points:1'  #Typically 0 for a radial cut, 1 for an axial cut
simlabel = "Axial-Simulation"
offset = 0e-6 #For if you didnt set up your OpenFOAM to have nice axes
trim = 38 #For if you want only the center region or data is NaN at edges
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
    
    dist=dist*1e6
    
    arr2 = arr[trim:-trim]
    dist2 = dist[trim:-trim]
    
#   gfit = ThrDim.FitDataGaussian(arr2, dist2, guess)
    #tfit = ThrDim.FitDataDoubleTanh(arr2, dist2, guess_tanh)
    tfit = ThrDim.FitDataDoubleTanhAbs(arr2, dist2, guess_tanh)
#   lfit = ThrDim.FitDataLorentz(arr2, dist2, guess)
    #sfit = ThrDim.FitDataSuperGaussian(arr2, dist2, guess_super)
    """
    plt.plot(dist,arr,label=simlabel)
    plt.title("Lorentzian Fit")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.plot(dist, ThrDim.Lorentz(lfit,dist),label="Lorentzian")
    plt.grid(); plt.legend(); plt.show()

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
    fwhm = Prof.CalcFWHM(dist,ThrDim.DoubleTanh(tfit,dist))
    flattop = Prof.CalcTopWidth(dist,ThrDim.DoubleTanh(tfit,dist),.95)
    plt.plot(dist,arr, label = 'Simulation Profile')
    plt.plot(dist, ThrDim.DoubleTanh(tfit,dist),label="Double Tanh")
    plt.title("Distribution")
    plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    plt.axvspan(-fwhm/2, -flattop/2, facecolor = 'b', alpha = 0.3, label='FWHM: '+str(fwhm))
    plt.axvspan(-flattop/2, flattop/2, facecolor = 'g', alpha = 0.3, label='95% Len: '+str(flattop))
    plt.axvspan(flattop/2, fwhm/2, facecolor = 'b', alpha = 0.3)
    plt.grid(); plt.legend(); plt.show()
    
    #plt.plot(dist, ThrDim.SuperGaussian(sfit,dist) - ThrDim.Gaussian(gfit,dist))