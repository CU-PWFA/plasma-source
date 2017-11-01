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

c=2.998e8
P=60e9
n=1
epsilon0 = 8.854e-12
#cm_m=100
wavelength = 785.3e-9
w0 = 5e-3
setupTitle = "old"

#I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)
I0 = 2*P/np.pi/np.power(w0,2)
Ei = np.sqrt(2*I0/c/n/epsilon0)*1e-9

l_step = 5e-5
zi = 10e-3
window = zi
zoom=int(round(window/l_step))

path = '/home/chris/Desktop/FourierPlots/ArJets/'
directory = 'Ar1_Near6/'
save = 0

choice=4
#Really the only choice is 4
if choice==4:
    setupTitle = "Spherical_2Cylindrical"
    f1 = 0.100 #m
    f2 = 0.0125
    L1 = 0.030
    L2 = 0.088
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

if choice==3:#Most simple, but optic is close to gas jet
    f = 0.030 #m
    #  30mm is good, but expensive
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.500,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,1.000)
    GB.Prop_CylindricalLens(q_y,q_z,1.000)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.5-f,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-2*f)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)
    q_x=q_y; q_y=q_z

if choice==2: #A bunch, but it works well
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.500,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,1.000)
    GB.Prop_CylindricalLens(q_y,q_z,1.000)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.150,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.030,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-.2)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.020,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-.2)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.21,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-.1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)
    q_x=q_y; q_y=q_z

if choice==1:
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.5,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,1)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.1,l_step)
    
    GB.Prop_CylindricalLens(q_y,q_z,.3)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.159,l_step)

    GB.Prop_CylindricalLens(q_y,q_z,.018)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)
    q_x=q_y; q_y=q_z

if choice==0: #1 meter, 3 lens
    #Focuses wz over a distance of 1m with 2 cyl. lens to .35mm
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_CylindricalLens(q_y, q_z, .8)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.35,l_step)
    GB.Prop_CylindricalLens(q_y, q_z, -.105)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.15,l_step)
    
    GB.Prop_CylindricalLens(q_z, q_y, 1)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,1,l_step)
    q_x=q_y; q_y=q_z

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
print()
phase = GB.Prop_EPhase(q_x,q_y,zoom,wavelength,Ei,w0)
#phase = [E0,wx,wy,Px,Py,phi,zi]

scl = 1e6 #Robert's code uses micrometers

pulseParams = {'Description' : setupTitle,
               'Location' : path+directory,
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
    if not os.path.exists(pulseParams['Location']):
        os.makedirs(pulseParams['Location'])
    np.save(pulseParams['Location']+'pulseParams',pulseParams)