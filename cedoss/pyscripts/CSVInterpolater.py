#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 11:19:02 2017

play around with Paraview clip scaling to get the correct dimensions,
and you want a little more in each direction so that your cartesian
grid is fully enclosed by the Paraview clip

Most recently, .6 .1 .03 seems to work pretty well to get a center region
the same size as the refraction code

@author: chris
"""

import csv
import numpy as np
from scipy.interpolate import griddata
import sys
import os
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

directory = '/home/chris/Desktop/CSVFiles/'
#filename = '3DLES_Wide_Axial.csv'
#filename = 'p2g8_box0.csv'
filename = 'WedgeRAS_R2_3DLowDen0.csv'
path = directory + filename

save_data = 0

rx=60e2
ry=24e2
rz=8e4

min_x = -rx/2; min_y = -ry/2; min_z = -rz/2
max_x =  rx/2; max_y =  ry/2; max_z =  rz/2

reducer = 1
nx = 2**(9-reducer)
ny = 2**(9-reducer)
nz = 2**(8-reducer)

bgr = 1e13 / 1e17

grid_x, grid_y, grid_z = np.mgrid[min_x:max_x:nx*1j, min_y:max_y:ny*1j, min_z:max_z:nz*1j]

variable = 'Density'
xaxis = 'Points:0'
yaxis = 'Points:1'
zaxis = 'Points:2'

arr=[]
x=[]; y=[]; z=[]

with open(path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        arr.append(float(row[variable])/1e17)
        x.append(float(row[xaxis])*1e6)
        y.append(float(row[yaxis])*1e6)
        z.append(float(row[zaxis])*1e6)
    arr=np.array(arr)
    x=np.array(x); y=np.array(y); z=np.array(z)
    
points = np.array([x,y,z]).T
values = arr

#grid_z = griddata(points, values, (grid_x, grid_y, grid_z), method='nearest')
grid_z = griddata(points, values, (grid_x, grid_y, grid_z), method='linear', fill_value = bgr)

if save_data == 1:
    sfolder = '/home/chris/Desktop/FourierPlots/ArJets/'
    sdirectory = 'WedgeRAS_R2/'
    spath = sfolder + sdirectory
    if not os.path.exists(spath):
        os.makedirs(spath)
    np.save(spath+'initDensity.npy', grid_z)

h = ThrDim.RobertRoll(grid_z)
y = np.linspace(-rx/2, rx/2, nx, False)
z = np.linspace(-ry/2, ry/2, ny, False)
x = np.linspace(-rz/2, rz/2, nz, False)
ThrDim.ImageCut(h,x,y,z,0,0,0,1e-3,'(mm)','Gas Density','e17(cm^-3)',1)