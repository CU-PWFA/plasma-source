# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 00:39:44 2017

Tests the propagation portion of GaussianBeam.py.  Several cases of
different combinations of optics (selected with the choice variable)

@author: chris
"""
import sys
sys.path.insert(0, "../")

from modules import GaussianBeam
import numpy as np
import matplotlib.pyplot as plt

#Choice for Testing Gaussian Beam Propagation, 1-5 use old parameters
choice = 8
print(choice)

#Define some values for wavelength, and beam waist
if choice <= 5:
    wavelength = 100e-9
    w0 = np.sqrt(0.159*wavelength/2)*0.4
else: #Using more realistic parameters
    wavelength = 500e-9
    w0 = 1e-3
    
#Initial q:  complex q parameter, distance, index of refraction
q=GaussianBeam.Prop_Init_q(wavelength, w0, 0, 1)

if choice == 1:
    #4 Thin Lens
    q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
    q=GaussianBeam.Prop_ThinLens(q,10)
    q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
    q=GaussianBeam.Prop_ThinLens(q,10)
    q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
    q=GaussianBeam.Prop_ThinLens(q,10)
    q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
    q=GaussianBeam.Prop_ThinLens(q,10)
    q=GaussianBeam.Prop_FreeSpace(q,5,0.1)

if choice == 2:
    #Flat Interface with n from 1 to 1.5
    q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)
    q=GaussianBeam.Prop_FlatInterface(q,1.5)
    q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)

if choice == 3:
    #Thick Lens - n = 1.8, R = +/- 2, d = 1
    q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)
    q=GaussianBeam.Prop_ThickLens(q,1.8,1,2,-2,0.1)
    q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)
    
if choice == 4:
    #Curved Interface - n = 1.8, R = +2
    q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)
    q=GaussianBeam.Prop_CurvedInterface(q,1.8,2)
    q=GaussianBeam.Prop_FreeSpace(q,3.6,0.1)
    
if choice == 5:
    #4 Thin Lens with Custom (Matches Choice 1)
    looper = range(50)
    for x in looper:
        q=GaussianBeam.Prop_Custom(q,1,0.1,0,1,0.1,-1)
    q=GaussianBeam.Prop_Custom(q,1,0,-1/10,1,0,-1)
    for x in looper:
        q=GaussianBeam.Prop_Custom(q,1,0.1,0,1,0.1,-1)
    q=GaussianBeam.Prop_Custom(q,1,0,-1/10,1,0,-1)
    for x in looper:
        q=GaussianBeam.Prop_Custom(q,1,0.1,0,1,0.1,-1)
    q=GaussianBeam.Prop_Custom(q,1,0,-1/10,1,0,-1)
    for x in looper:
        q=GaussianBeam.Prop_Custom(q,1,0.1,0,1,0.1,-1)
    q=GaussianBeam.Prop_Custom(q,1,0,-1/10,1,0,-1)
    for x in looper:
        q=GaussianBeam.Prop_Custom(q,1,0.1,0,1,0.1,-1)
    
if choice == 6:
    #Thick lens at 1m, focuses beam around 3.5m
    q=GaussianBeam.Prop_FreeSpace(q,1,1e-5)
    q=GaussianBeam.Prop_ThickLens(q,1.4,1e-3,8,-8,1e-5)
    q=GaussianBeam.Prop_FreeSpace(q,5,1e-5)
    
if choice == 7:
    #Thick lens at 1m, focuses beam around 3.5m
    q=GaussianBeam.Prop_Init_q(wavelength, w0, -2, 1)
    q=GaussianBeam.Prop_FreeSpace(q,2,1e-5)
    q=GaussianBeam.Prop_ThickLens(q,1.4,1e-3,8,-8,1e-5)
    q=GaussianBeam.Prop_FreeSpace(q,6,1e-5)
    
if choice == 8:
    #Using cylindrical lens with axis in z direction
    q_y = list(q)
    q_z = list(q)
    q_y=GaussianBeam.Prop_FreeSpace(q_y,1,1e-5)
    q_z=GaussianBeam.Prop_FreeSpace(q_z,1,1e-5)
    GaussianBeam.Prop_CylindricalLens(q_z, q_y, 8)
    q_y=GaussianBeam.Prop_FreeSpace(q_y,6,1e-5)
    q_z=GaussianBeam.Prop_FreeSpace(q_z,6,1e-5)
    
#Gets the range over which beam propagated, as well as 
# radius of curvature, R, and beam spot size, w
if choice < 8:
    xrange=GaussianBeam.Prop_GetRange(q)
    R=GaussianBeam.Prop_RadiusList(q)
    R[0]=0  #The radius of curvature initially --> infinity
    w=GaussianBeam.Prop_SpotList(q,wavelength)

    #Plots spot size
    plt.plot(xrange,w)
    plt.title("Spot size from q parameter")
    plt.ylabel("Spot Size (m)")
    plt.xlabel("Distance (m)")
    plt.show()

    #Plots radius of curvature (Value can be very 
                            # large for when beam focuses on waist)
    plt.plot(xrange,R)
    plt.title("Curvature from q parameter")
    plt.ylabel("Radius of Curvature (m)")
    plt.xlabel("Distance (m)")
    plt.show()

if choice >= 8:
    xrange=GaussianBeam.Prop_GetRange(q_y)
    wy=GaussianBeam.Prop_SpotList(q_y,wavelength)
    wz=GaussianBeam.Prop_SpotList(q_z,wavelength)

    #Plots spot size
    plt.plot(xrange,wy)
    plt.title("Spot size from q parameter:  blue transverse, orange parallel")
    plt.ylabel("Spot Size (m)")
    plt.xlabel("Distance (m)")
    plt.plot(xrange,wz)
    plt.show()