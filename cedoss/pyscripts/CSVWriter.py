#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 11:49:14 2017

read a numpy file and export a csv

@author: chris
"""
import sys
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim
import csv
import numpy as np

directory = '/home/chris/Desktop/figtester/'
npfilename = 'finalDensity.npy'
csvfilename= 'finalDensity.csv'

arr = np.load(directory+npfilename)
params = np.load(directory+'params.npy').item()

X = params['X']; Nx = params['Nx']
Y = params['Y']; Ny = params['Ny']
Z = params['Z']; Nz = params['Nz']

arr = ThrDim.RobertRoll(arr)
y = np.linspace(-X/2, X/2, Nx, False)
z = np.linspace(-Y/2, Y/2, Ny, False)
x = np.linspace(-Z/2, Z/2, Nz, False)

with open(directory+csvfilename,'w') as wrt:#csvfile:
    #wrt = csv.writer(csvfile)
    wrt.write("Density,x,y,z\n")
    for i in range(len(x))[int(len(x)/3):int(len(x)*2/3)]:
        xi = x[i]
        for j in range(len(y)):
            yj = y[j]
            for k in range(len(z)):
                arri = arr[i,j,k]
                row = str(arri)+","+str(xi)+","+str(yj)+","+str(z[k])+"\n"
                wrt.write(row)