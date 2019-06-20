#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 16:01:38 2019

Import data from a csv written by paraview for the purpose of fitting 
data (probably density) to a function (probably Exponential)

From paraview, take a axial "plot over line" and then export to csv using
file -> save data.  Depending on the axis, may need to adjust the # below in
axis = 'Points:#' from 0 to 2.
    
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
off = 0
fac = -1
filename = '3DLam_R3_DenAxial.csv'

path = directory + filename

variable = 'Den'
axis = 'Points:1'  #Typically 0 for a radial cut, 1 for an axial cut
simlabel = "Axial-Simulation"
offset = -25000e-6 #For if you didnt set up your OpenFOAM to have nice axes
zoom = 1e6 #For when you have meters and work with micrometers

if off == 1:
    filename = 'WedgeRAS_R2_LowDenAxial_off.csv'
    path = directory + filename
    fac = -1
    axis = 'Points:0'  #Typically 0 for a radial cut, 1 for an axial cut

#guess = [1e15,0,1e12]
guess = [1.91288587e+20, 3.51434237e+03, -4.53692351e+02, 3.79269734e+13]
arr = []
dist = []

with open(path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        arr.append(float(row[variable]))
        dist.append(fac*float(row[axis])-offset)
        
    arr = np.array(arr)
    dist = np.array(dist)
    
    #arr = arr[round(len(arr)*1/4):round(len(arr)*3/4)]
    #dist = dist[round(len(dist)*1/4):round(len(dist)*3/4)]
    
    dist=-1*dist*1e6
#    arr=arr-arr[0]
    efit=ThrDim.FitDataSomething(arr, dist, ThrDim.LorentzOffset, guess)
    
    #xwindow = dist[-1]; xstep = (dist[2]-dist[1])/50 #microns
    #dist_precise = np.arange(-xwindow, xwindow, xstep)
    
    #plt.plot(dist,arr,label=simlabel)
    #plt.title("Data")
    #plt.xlabel("Distance (microns)"); plt.ylabel("Density (cm^-3)")
    #plt.plot(dist, ThrDim.Lorentz(lfit,dist),label="Lorentzian")
    #plt.grid(); plt.legend(); plt.show()