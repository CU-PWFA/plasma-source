#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 11:49:14 2017

Reads Robert's output in the form of a folder with a ton of .npy's

@author: chris
"""
import sys
sys.path.insert(0, "../")
import numpy as np

directory = '/home/chris/Desktop/element_Plasma/'
#npfilename = 'Plasma_plasmaDensity_60.npy'
csvfilename= 'finalDensity.csv'

#arr = np.load(directory+npfilename)
params = np.load(directory+'Plasma_params.npy').item()

X = params['X']; Nx_r = params['Nx']
Y = params['Y']; Ny_r = params['Ny']
Z = params['Z']; Nz_r = params['Nz']

rescale = 1/62.5 #2X/Z

#arr = ThrDim.RobertRoll(arr)
x = np.linspace(-Z/2*rescale, Z/2*rescale, Nz_r, False)
r = np.linspace(0, X/2, Nx_r, False)
#z = np.linspace(-Y/2, Y/2, Ny_r, False)
Nx = Nz_r; Nr = Nx_r; Np = 256
phi = np.linspace(0,2*np.pi,Np)

arr2d=np.zeros((Nx,Nr))

for i in range(1,Nx):
    arr2d[i]=np.load(directory+"Plasma_plasmaDensity_"+str(i)+".npy")
arr2d = arr2d[:,round(Nr/2):]

rmax = 0; Nr=round(Nr/2)
for i in range(Nx-1):
    for j in range(Nr):
        if (arr2d[i,j] < 1e-04):
            if (j > rmax):
                rmax = j
            break;
r=r[:rmax]
arr2d=arr2d[:,:rmax]

#MAKE SURE TO CHANGE THE X LOOP TO TRUNCATE WHAT IS NEEDED
with open(directory+csvfilename,'w') as wrt:#csvfile:
    #wrt = csv.writer(csvfile)
    wrt.write("Density,x,y,z\n")
    for i in range(len(x)):
        xi = x[i]
        for j in range(len(r)):
            for k in range(len(phi)):
                yj = r[j]*np.sin(phi[k])
                zk = r[j]*np.cos(phi[k])
                arri = arr2d[i,j]
                row = str(arri)+","+str(xi)+","+str(yj)+","+str(zk)+"\n"
                wrt.write(row)
