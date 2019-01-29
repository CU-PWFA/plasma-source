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
from modules import TPLFocalLength as Foc

c=2.998e8               #Speed of light constant
P=60e9                  #Laser power in W
n=1.                    #Refractive index
epsilon0 = 8.854e-12    #Dielectric constant
wavelength = 785.3e-9   #Laser wavelength in m
w0 = 5e-3               #Initial spot size in m
E_ion = 13.6            #Ionization Energy:  13.6 for H, 15.8 for Ar, 
setupTitle = ""         #Name the setup 
reverse = 0

#I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)

l_step = 5e-5           #dz in m
zi = 7.5e-3             #Offset from small waist to plane we want to save params at
zoom=int(round(zi/l_step))

path = '/home/chris/Desktop/DataLoads/PulseFilesNp/'
filename = 'pulseParams_rotate.npy'

save = 0                #Set to 1 to save anything
calcdensity = 1         #Set to 1 to calc resulting plasma density w/out refraction
calcfocal = 1

foc_dom_fac = 2
radscl = 1   #Set to larger to increase the beam axis domain

choice=37         #Set to one of the setups below

if choice==38:#Case 6 but opposite  Gives l = 775 um 
    setupTitle = "Spherical_2Cylindrical_reversed"
    reverse = 1
    zi = 5e-2             #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))
    radscl = 3
    
    setupTitle = "Spherical_2Cylindrical"
    f1 = 0.075 #m
    f2 = 0.0125
    L1 = 0.037
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
    GB.Prop_CylindricalLens(q_y,q_x,-2*f1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2-L1,l_step) 
    GB.Prop_CylindricalLens(q_x,q_y,-2*f2)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, 0.5,l_step)

if choice==37:#Gives l = 775 um 
    setupTitle = "Spherical_2Cylindrical_reversed"
    reverse = 1
    zi = 5e-2             #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))
    radscl = 3
    
    setupTitle = "Spherical_2Cylindrical"
    f1 = 0.075 #m
    f2 = 0.025
    L1 = 0.037
    L2 = 0.076
    P=83e9 #83 GW for 597 um, 218 GW for 996 um
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
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, 0.5,l_step)

if choice==35:#Case 5 but opposite  Gives l = 174 um
    setupTitle = "Spherical_2Cylindrical_reversed"
    choice = 5
    reverse = 1
    zi = 1.5e-2             #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))
    
if choice==36:#Case 6 but opposite  Gives l - 199 um 
    setupTitle = "Spherical_2Cylindrical_reversed"
    choice = 6
    reverse = 1
    zi = 1e-2             #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))

if choice==6:#Case 4 but wider ~40 um
    setupTitle = "Spherical_2Cylindrical"
    f1 = 0.075 #m
    f2 = 0.0125
    L1 = 0.037
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
    GB.Prop_CylindricalLens(q_y,q_x,-2*f1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2-L1,l_step) 
    GB.Prop_CylindricalLens(q_x,q_y,-2*f2)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.4,l_step)

if choice==5:#Case 4 but wider ~40 um
    setupTitle = "Spherical_2Cylindrical"
    f1 = 0.040 #m
    f2 = 0.0125
    L1 = 0.064 # 0.065
    L2 = 0.088
    P=66e9 # 80e9
    
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

#Was my favorite setup in the beginning
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
The 20+ series is a quest for 100 um lens.  The issue is that you need
to ionize over a range of 24cm, waaay too long for a TPL
"""
if choice==20:
    setupTitle = "Spherical_2Cylindrical_200"
    f2 = -0.0125
    P=100e9
    zi = 0.1             #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    #10 cm to narrow waist focuser
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, .1, l_step)
    GB.Prop_CylindricalLens(q_y,q_x,.200)
    
    #8.1 cm to narrow waist defocuser
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, .081, l_step)
    GB.Prop_CylindricalLens(q_y,q_x,-.04)
    
    #1.9 cm to wide waist focuser
    GB.Prop_Cylindrical_FreeSpace(q_y,q_x, .019, l_step)
    GB.Prop_CylindricalLens(q_x,q_y,.200)
    
    #9 cm to wide wast defocuser
    GB.Prop_Cylindrical_FreeSpace(q_y,q_x, .09, l_step)
    GB.Prop_CylindricalLens(q_x,q_y,-.02)
    
    #Ride it out
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.45,l_step)

if choice==21: #100 um waist, but spherical only
    setupTitle = "Spherical_2Cylindrical_200"
    f2 = -0.0125
    P=120e9
    zi = 0.1             #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, .1, l_step)
    GB.Prop_CylindricalLens(q_y,q_x,.200)
    GB.Prop_CylindricalLens(q_x,q_y,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, .081, l_step)
    GB.Prop_CylindricalLens(q_y,q_x, -.04)
    GB.Prop_CylindricalLens(q_x,q_y, -.04)
    
    #GB.Prop_Cylindrical_FreeSpace(q_x,q_y, 0.05, l_step) 
    #GB.Prop_CylindricalLens(q_x,q_y, 2 * -0.015)
    
    #GB.Prop_Cylindrical_FreeSpace(q_x,q_y, 0.02, l_step) 
    #GB.Prop_CylindricalLens(q_y,q_x,2*f2)#narrow
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.5,l_step)
    
if choice==22:#Fairly good ~70 um waist
    setupTitle = "Spherical_2Cylindrical_200"
    f2 = -0.0125
    P=50e9
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,.200)
    GB.Prop_CylindricalLens(q_x,q_y,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, .072, l_step)
    GB.Prop_CylindricalLens(q_y,q_x, -.06)
    GB.Prop_CylindricalLens(q_x,q_y, -.06)
    
    #GB.Prop_Cylindrical_FreeSpace(q_x,q_y, 0.05, l_step) 
    #GB.Prop_CylindricalLens(q_x,q_y, 2 * -0.015)
    
    #GB.Prop_Cylindrical_FreeSpace(q_x,q_y, 0.02, l_step) 
    #GB.Prop_CylindricalLens(q_y,q_x,2*f2)#narrow
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.5,l_step)
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
if reverse == 1:
    plt.plot(zrange_tot,wtoty,label="Wide Waist Spot Size")
plt.plot(zrange_tot,wtotx,label="Narrow Waist Spot Size")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Laser Axis (z)")
if reverse != 1:
    plt.plot(zrange_tot,wtoty,label="Wide Waist Spot Size")
plt.legend()
plt.grid(); plt.show()

#Zoom in to a scale of the rayleigh length of wy
waist=wtotx.index(min(wtotx))
wx=np.array(wtotx[(waist-zoom):(waist+zoom)])
wy=np.array(wtoty[(waist-zoom):(waist+zoom)])
xrange=zrange_tot[(waist-zoom):(waist+zoom)]

if reverse == 1:
    q_d = q_y
    q_y = q_x
    q_x = q_d
    wd = np.copy(wy)
    wy = np.copy(wx)
    wx = wd

#Plots spot size zoomed in around waist
plt.plot(xrange,wx*1e6,label="Beam Axis Waist (x)")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (um)")
plt.xlabel("Laser Axis (z)")
plt.plot(xrange,wy*1e6,label="Transverse Waist (y)")
plt.legend()
plt.grid(); plt.show()

#Print some diagnostics for minimum spot sizes
GB.Prop_SpotInfo(q_x,wavelength,'w0x','z(w0x)')
GB.Prop_SpotInfo(q_y,wavelength,'w0y','z(w0y)')
if reverse != 1:
    print(str(wy[int(len(xrange)/2)]*1e6) + " at center")
else:
    print(str(wx[int(len(xrange)/2)]*1e6) + " at center")
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

if calcfocal == 1:
    den = 1
    #den = 10
    params = frefract.GetDefaultParams()
    X = params['X']*foc_dom_fac*radscl; Nx = params['Nx']
    y = np.linspace(-X/2, X/2, Nx, False)/1e6
    
    wmult = wy * wx
    wmult_min = np.argmin(wmult)
    
    #wz_center = wy[int(len(wy)/2)]
    #wy_center = wx[int(len(wx)/2)]
    print("center change: ",int(len(wy)/2), " to " ,wmult_min)
    wz_center = wy[wmult_min]
    wy_center = wx[wmult_min]
    
    
    I = GB.IntensityFromSpotSizes_1D(wy_center,wz_center,y,I0/1e4,w0)
    plt.plot(y*1e6,I)
    plt.title("Intensity along beam axis")
    plt.xlabel(r'$\mathrm{Beam\,\,Axis\,\,[\mu m]}$')
    plt.ylabel(r'$\mathrm{Intensity \, [W/cm^2]}$')
    plt.grid(); plt.show()
    
    params['EI'] = E_ion
    H = ThrDim.IonFracFromIntensity_1D(I,params['EI'],35e-15)
    H = H * den
    plt.plot(y*1e6,H)
    plt.title("Density along beam axis")
    plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
    plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
    plt.grid(); plt.show()
    
    focal = Foc.Calc_Focus(H, y*1e6)
    thick = Foc.Calc_Square_Lens(den*1e17, focal)
    print(thick,"um equivalent thickness")
    
    guess = [thick, thick/6, den]
    fit = ThrDim.FitDataDoubleTanh(H, y*1e6, guess)

if calcdensity == 1:
    den = 0.5
    params = frefract.GetDefaultParams()
    params['Nz'] = len(xrange)
    params['Nx'] = int(params['Nx']/4)
    params['Ny'] = int(params['Ny']/4)
    params['Z'] = zi*1e6*2
    
    X = params['X']*radscl; Nx = params['Nx']
    Y = params['Y']; Ny = params['Ny']
    Z = params['Z']; Nz = params['Nz']

    y = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(-Y/2, Y/2, Ny, False)
    x = np.linspace(-Z/2, Z/2, Nz, False)
    
    wz=np.array(wy)*1e6
    wy=np.array(wx)*1e6
    I = GB.IntensityFromSpotSizes(wy/1e6,wz/1e6,x/1e6,y/1e6,z/1e6,I0/1e4,w0)
    ThrDim.ImageCut(I,x,y,z,0,0,0,1e-3,'(mm)','Intensity','W/cm^2',0)
    params['EI'] = E_ion
    H = ThrDim.IonFracFromIntensity(I,params['EI'],35e-15)
    ThrDim.ImageCut(H,x,y,z,0,0,0,1e-3,'(mm)','Ion. Frac.','%',0)
    H = H * den
    H = ThrDim.RobertRoll(ThrDim.RobertRoll(H))
    frefract.TestDensity(H,params)
    
    if save == 1:
        savefolder = '/home/chris/Desktop/FourierPlots/ArVarySpot_Calc/'
        savefolder = savefolder + 'case_2/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        np.save(savefolder+'finalDensity.npy',H)
        np.save(savefolder+'params.npy',params)