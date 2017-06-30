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

import numpy as np
import matplotlib.pyplot as plt

from modules import GaussianBeam as GB
from modules import Doss_Ionization as adkD

c=2.998e8
P=60e9
cm_m=100
wavelength = 785.3e-9
w0 = 5e-3

I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)

l_step = 5e-5
window = 5e-3
zoom=int(round(window/l_step))

choice=4

if choice==4:
    f1 = 0.100 #m
    f2 = 0.0125
    L1 = 0.030
    L2 = 0.088
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.1,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,.200)
    GB.Prop_CylindricalLens(q_y,q_z,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,L1,l_step) 
    GB.Prop_CylindricalLens(q_z,q_y,-2*f1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,L2-L1,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-2*f2)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)

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

#Get the total domain of spot sizes
xrange_tot=GB.Prop_GetRange(q_y)
wtoty=GB.Prop_SpotList(q_y,wavelength)
wtotz=GB.Prop_SpotList(q_z,wavelength)

#Plots spot size
plt.plot(xrange_tot,wtoty,label="Spot size in y")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.plot(xrange_tot,wtotz,label="Spot size in z")
plt.legend()
plt.show()

#Zoom in to a scale of the rayleigh length of wy
waist=wtoty.index(min(wtoty))
wy=wtoty[(waist-zoom):(waist+zoom)]
wz=wtotz[(waist-zoom):(waist+zoom)]
xrange=xrange_tot[(waist-zoom):(waist+zoom)]

#Plots spot size zoomed in around waist
plt.plot(xrange,wy,label="Spot size in y")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.plot(xrange,wz,label="Spot size in z")
plt.legend()
plt.show()

#Print some diagnostics for minimum spot sizes
GB.Prop_SpotInfo(q_y,wavelength,'w0y','x(w0y)')
GB.Prop_SpotInfo(q_z,wavelength,'w0z','x(w0z)')
print()
E0=adkD.Elefield(I0)
GB.Prop_EPhase(q_y,q_z,zoom,wavelength,E0,w0)