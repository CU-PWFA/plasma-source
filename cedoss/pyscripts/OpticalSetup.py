#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:14:55 2017

Dedicated Script for testing optical setup and looking for best way to get
the desired spot sizes from cylindrical lens

@author: chris
"""

import sys
sys.path.insert(0, "../")
import os
import numpy as np
import matplotlib.pyplot as plt

from modules import GaussianBeam as GB
from modules import ThreeDimensionAnalysis as ThrDim
from modules import FourierRefraction as frefract

c=2.998e8
P=60e9
n=1
epsilon0 = 8.854e-12
#cm_m=100
wavelength = 785.3e-9
w0 = 5e-3
setupTitle = "old"

#I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)

l_step = 5e-5
zi = 7.5e-3
window = zi
zoom=int(round(window/l_step))

path = '/home/chris/Desktop/DataLoads/PulseFilesNp/'
filename = 'pulseParams_narrow25.npy'
save = 0

calcdensity = 1

choice=10
#Really the only choice is 4
if choice==4:
    setupTitle = "Spherical_2Cylindrical"
    f1 = 0.100 #m
    f2 = 0.0125
    L1 = 0.030
    L2 = 0.088
    P=60e9
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,.200)
    GB.Prop_CylindricalLens(q_x,q_y,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L1,l_step) 
    GB.Prop_CylindricalLens(q_y,q_x,-2*f1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2-L1,l_step) 
    GB.Prop_CylindricalLens(q_x,q_y,-2*f2)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.4,l_step)
"""
The 10+ series are alterations of #4 to get different spot sizes
"""
if choice==10:
    setupTitle = "Spherical_2Cylindrical_adjusted"
    f1 = -0.075 #m
    f2 = -0.0125
    L1 = 0.040
    L2 = 0.088
    P=50e9
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,.200)
    GB.Prop_CylindricalLens(q_x,q_y,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L1,l_step) 
    GB.Prop_CylindricalLens(q_y,q_x,2*f1)#narrow
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2-L1,l_step) 
    GB.Prop_CylindricalLens(q_x,q_y,2*f2)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.4,l_step)
    
if choice==11:
    setupTitle = "Spherical_2Cylindrical_adjusted"
    f1 = -0.075 #m
    f2 = -0.0125
    L1 = 0.045
    L2 = 0.088
    P=50e9
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,.200)
    GB.Prop_CylindricalLens(q_x,q_y,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L1,l_step) 
    GB.Prop_CylindricalLens(q_y,q_x,2*f1)#narrow
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2-L1,l_step) 
    GB.Prop_CylindricalLens(q_x,q_y,2*f2)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.4,l_step)

I0 = 2*P/np.pi/np.power(w0,2)
Ei = np.sqrt(2*I0/c/n/epsilon0)*1e-9

#Get the total domain of spot sizes
zrange_tot=GB.Prop_GetRange(q_x)
wtotx=GB.Prop_SpotList(q_x,wavelength)
wtoty=GB.Prop_SpotList(q_y,wavelength)

#Plots spot size
plt.plot(zrange_tot,wtotx,label="Narrow Waist Spot Size (x)")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Laser Axis (z)")
plt.plot(zrange_tot,wtoty,label="Wide Waist Spot Size (y)")
plt.legend()
plt.show()

#Zoom in to a scale of the rayleigh length of wy
waist=wtotx.index(min(wtotx))
wx=wtotx[(waist-zoom):(waist+zoom)]
wy=wtoty[(waist-zoom):(waist+zoom)]
xrange=zrange_tot[(waist-zoom):(waist+zoom)]

#Plots spot size zoomed in around waist
plt.plot(xrange,wx,label="Narrow Waist Spot Size (x)")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Laser Axis (z)")
plt.plot(xrange,wy,label="Wide Waist Spot Size (y)")
plt.legend()
plt.show()



#Print some diagnostics for minimum spot sizes
GB.Prop_SpotInfo(q_x,wavelength,'w0x','z(w0x)')
GB.Prop_SpotInfo(q_y,wavelength,'w0y','z(w0y)')
print(str(wy[int(len(xrange)/2)]*1e6) + " at center")
print()
phase = GB.Prop_EPhase(q_x,q_y,zoom,wavelength,Ei,w0)
#phase = [E0,wx,wy,Px,Py,phi,zi]

scl = 1e6 #Robert's code uses micrometers

pulseParams = {'Description' : setupTitle,
               'Location' : path+filename,
               'power' : P,
               'wavelength' : wavelength,
               'w0' : w0,
               'E0' : phase[0],
               'wx' : phase[1] * scl,
               'wy' : phase[2] * scl,
               'Px' : phase[3] * scl * scl,
               'Py' : phase[4] * scl * scl,
               'phi' : phase[5],
               'zi' : zi * scl
               }

if save == 1:
    #Create the folder if it not already exists
    if not os.path.exists(path):
        os.makedirs(path)
    np.save(pulseParams['Location'],pulseParams)
    
if calcdensity == 1:
    params = frefract.GetDefaultParams()
    params['Nz'] = len(xrange)
    #params['Nx'] = int(params['Nx']/4)
    #params['Ny'] = int(params['Ny']/4)
    
    X = params['X']; Nx = params['Nx']
    Y = params['Y']; Ny = params['Ny']
    Z = params['Z']; Nz = params['Nz']

    y = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(-Y/2, Y/2, Ny, False)
    x = np.linspace(-Z/2, Z/2, Nz, False)
    
    wz=np.array(wy)*1e6
    wy=np.array(wx)*1e6
    I = GB.IntensityFromSpotSizes(wy/1e6,wz/1e6,x/1e6,y/1e6,z/1e6,I0/1e4,w0)
    ThrDim.ImageCut(I,x,y,z,0,0,0,1e-3,'(mm)','Intensity','W/cm^2',0)
    H = ThrDim.IonFracFromIntensity(I,params['EI'],35e-15)
    H = H * 0.1
    H = ThrDim.RobertRoll(ThrDim.RobertRoll(H))
    frefract.TestDensity(H,params)
    
    if save == 1:
        savefolder = '/home/chris/Desktop/FourierPlots/ArVarySpot_Calc/'
        savefolder = savefolder + 'case_2/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        np.save(savefolder+'finalDensity.npy',H)
        np.save(savefolder+'params.npy',params)
    
    
    
    
    
    
    
    
    
