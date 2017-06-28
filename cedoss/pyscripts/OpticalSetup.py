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

from modules import GaussianBeam
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

choice=0
if choice==0: #1 meter, 3 lens
    #Focuses wz over a distance of 1m with 2 cyl. lens to .35mm
    q_y = GaussianBeam.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GaussianBeam.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GaussianBeam.Prop_CylindricalLens(q_y, q_z, .8)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,.35,l_step)
    GaussianBeam.Prop_CylindricalLens(q_y, q_z, -.105)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,.15,l_step)
    
    GaussianBeam.Prop_CylindricalLens(q_z, q_y, 1)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,1,l_step)

#Get the total domain of spot sizes
xrange_tot=GaussianBeam.Prop_GetRange(q_y)
wtoty=GaussianBeam.Prop_SpotList(q_y,wavelength)
wtotz=GaussianBeam.Prop_SpotList(q_z,wavelength)

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
GaussianBeam.Prop_SpotInfo(q_y,wavelength,'w0y','x(w0y)')
GaussianBeam.Prop_SpotInfo(q_z,wavelength,'w0z','x(w0z)')
print()
E0=adkD.Elefield(I0)
GaussianBeam.Prop_EPhase(q_y,q_z,zoom,wavelength,E0,w0)