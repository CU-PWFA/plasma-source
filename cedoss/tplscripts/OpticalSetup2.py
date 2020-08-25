#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:07:31 2020

Dedicated Script for testing optical setup and looking for best way to getthe desired spot sizes from cylindrical lens

This is the sequel version, created after added M-squared of non-ideal
beams into the analysis

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
P=60e9                  #Default laser power in W
n=1.                    #Refractive index
epsilon0 = 8.854e-12    #Dielectric constant
wavelength = 800.0e-9 #785.3e-9   #Laser wavelength in m
w0 = 5e-3               #Initial spot size in m
E_ion = 24.6           #Ionization Energy:  13.6 for H, 15.8 for Ar, 24.6 for He 27.6 for Ar++
setupTitle = ""         #Name the setup 
reverse = 0
fs_duration = 40e-15

M2 = 2#np.sqrt(2)

#I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)

l_step = 5e-5           #dz in m
zi = 7.5e-3             #Offset from small waist to plane we want to save params at
zoom=int(round(zi/l_step))

path = '/home/chris/Desktop/DataLoads/PulseFilesNp/'
filename = 'pulseParams_testm22.npy'

save = 0             #Set to 1 to save anything
calcdensity = 1 #SLOW   #Set to 1 to calc resulting plasma density w/out refraction
calcfocal = 0

foc_dom_fac = 2
radscl = 1   #Set to larger to increase the beasm axis domain

choice=592         #Set to one of the setups below

if choice == 83:  #To be designed for the limits of a window at 38.67 cm or 21.27 cm for <1e12 W/cm2
    setupTitle = "Spherical_2Cylindrical_reversed"
    w0 = 8.66e-3 #calculated from 30mm diameter flattop estimation
    reverse = 0
    zi = 1.0e-2#2e-2#38.67e-2            #Offset from small waist to plane we want to save params at
    zoom=int(round(zi/l_step))
    radscl = 0.4
    M2 = 4.6
    setupTitle = "Spherical_2Cylindrical"

    f0 = 0.20

    f2 = -0.1000 #m
    L2 = 0.109

    f1 = 0.700 #m
    L1 = 0.420

    Lmax = 2.5

    P = 300e9

     #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)

    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    #GB.Prop_CylindricalLens(q_y,q_x,2*f0)
    GB.Prop_CylindricalLens(q_x,q_y,2*f0)

    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2,l_step) 
    GB.Prop_CylindricalLens(q_x,q_y,2*f2)

    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L1-L2,l_step) 
    GB.Prop_CylindricalLens(q_y,q_x,2*f1)

    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, Lmax,l_step)

if choice == 591:  #Single Spherical Lens
    setupTitle = "1Spherical"
    w0 = 14.3e-3 #calculated from 30mm diameter flattop estimation
    #w0 = 3.575e-3
    reverse = 0
    zi = 2e-2
    zoom = int(round(zi/l_step))
    radscl = 0.3
    
    f0 = 0.646
    Lmax = 1.0
    #P = 12.5e9
    P = 125e9
    
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    
    GB.Prop_CylindricalLens(q_y,q_x,2*f0)
    GB.Prop_CylindricalLens(q_x,q_y,2*f0)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,Lmax,l_step)
    
if choice == 592:  #Single Spherical Lens, same as above but zoomed in more
    setupTitle = "1Spherical"
    w0 = 14.3e-3 #calculated from 30mm diameter flattop estimation
    #w0 = 3.575e-3
    reverse = 0
    zi = 0.40e-2#0.75e-2#2e-2
    zoom = int(round(zi/l_step))
    radscl = 0.3
    
    f0 = 0.646
    Lmax = 1.0
    #P = 12.5e9
    P = 125e9
    
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    
    GB.Prop_CylindricalLens(q_y,q_x,2*f0)
    GB.Prop_CylindricalLens(q_x,q_y,2*f0)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,Lmax,l_step)
    
if choice == 593:  #Single Spherical Lens, large M2 and tons of power
    setupTitle = "1Spherical"
    w0 = 4.2e-3
    M2 = 16*np.sqrt(2)
    reverse = 0
    zi = 2.00e-2#0.75e-2#2e-2
    zoom = int(round(zi/l_step))
    radscl = 1.0
    
    f0 = 0.646
    Lmax = 1.0
    #P = 12.5e9
    P = 1250e9
    
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    
    GB.Prop_CylindricalLens(q_y,q_x,2*f0)
    GB.Prop_CylindricalLens(q_x,q_y,2*f0)
    
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,Lmax,l_step)
    
if choice == 600:  #Single Spherical Lens, same as above but zoomed in more
    setupTitle = "1Spherical"
    w0 = 4.2e-3
    reverse = 0
    zi = 2.00e-2#0.75e-2#2e-2
    zoom = int(round(zi/l_step))
    radscl = 0.3
    
    f1 = 1.0
    f2 = 0.3
    Lmax = 1.0
    #P = 12.5e9
    P = 250e9
    
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_x,q_y,2*f1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,0.7,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,2*f2)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,0.4,l_step)
    
if choice == 601:  #Single Spherical Lens, same as above but zoomed in more
    setupTitle = "1Spherical"
    w0 = 4.2e-3
    reverse = 0
    zi = 2.00e-2#0.75e-2#2e-2
    zoom = int(round(zi/l_step))
    radscl = 0.3
    
    f1 = 1.5
    f2 = 0.3
    Lmax = 1.0
    #P = 12.5e9
    P = 250e9
    
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_x,q_y,2*f1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,1.2,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,2*f2)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,0.4,l_step)
    
if choice == 602:  #Single Spherical Lens, same as above but zoomed in more
    setupTitle = "1Spherical"
    w0 = 16.8e-3#4.2e-3
    reverse = 0
    zi = 0.30e-2#0.75e-2#2e-2
    zoom = int(round(zi/l_step))
    radscl = 0.3
    
    f1 = 2.0
    f2 = 0.3
    Lmax = 1.0
    #P = 12.5e9
    P = 25e9
    
    q_x = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    q_y = GB.Prop_Init_q(wavelength, w0, -.1, 1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,.1,l_step)
    GB.Prop_CylindricalLens(q_x,q_y,2*f1)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,1.7,l_step)
    GB.Prop_CylindricalLens(q_y,q_x,2*f2)
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,0.4,l_step)

I0 = 2*P/np.pi/np.power(w0,2)
Ei = np.sqrt(2*I0/c/n/epsilon0)*1e-9

#Get the total domain of spot sizes
zrange_tot=GB.Prop_GetRange(q_x)
wtotx=GB.Prop_SpotList(q_x,wavelength)
wtoty=GB.Prop_SpotList(q_y,wavelength)

#Plots spot size
fig, ax1 = plt.subplots(figsize=(14,4))
#plt.rcParams.update({'font.size': 12})
plt.subplot(121)
if reverse == 1:
    plt.plot(zrange_tot,wtoty,c='orange',ls='solid',label="Wide Waist Spot Size")
    if M2 != 1:
        plt.plot(zrange_tot,np.array(wtoty)*np.sqrt(M2),c='orange',ls='--',label="Wide M2")
plt.plot(zrange_tot,wtotx,c='blue',ls='solid',label="Narrow Waist Spot Size")
if M2 != 1:
    plt.plot(zrange_tot,np.array(wtotx)*np.sqrt(M2),c='blue',ls='--',label="Narrow M2")
plt.title("Spot size from q parameter")
plt.ylabel(r'Spot Size [$m$]')
plt.xlabel(r'Laser Axis (x) [$m$]')
if reverse != 1:
    plt.plot(zrange_tot,wtoty,c='orange',ls='solid',label="Wide Waist Spot Size")
    if M2 != 1:
        plt.plot(zrange_tot,np.array(wtoty)*np.sqrt(M2),c='orange',ls='--',label="Wide M2")
plt.legend()
plt.grid()#; plt.show()

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
plt.subplot(122)
plt.plot(xrange,wx*1e6,c='blue',ls='solid',label="Beam Axis Waist (x)")
if M2 != 1:
    plt.plot(xrange,wx*np.sqrt(M2)*1e6,c='blue',ls='--',label="Beam Axis w/ M2 (x)")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (um)")
plt.xlabel("Laser Axis (x)")
plt.plot(xrange,wy*1e6,c='orange',ls='solid',label="Transverse Waist (y)")
if M2 != 1:
    plt.plot(xrange,wy*np.sqrt(M2)*1e6,c='orange',ls='--',label="Transverse w/ M2 (y)")
plt.legend()
plt.grid(); plt.show()

print(zoom)
print(zrange_tot[waist-zoom])
print(xrange[0])

#Print some diagnostics for minimum spot sizes
GB.Prop_SpotInfo(q_x,wavelength,'w0x','z(w0x)')
GB.Prop_SpotInfo(q_y,wavelength,'w0y','z(w0y)')
if reverse != 1:
    print(str(wy[int(len(xrange)/2)]*1e6) + " x " + str(wx[int(len(xrange)/2)]*1e6) + " at center")
    if M2 != 1:
        print(str(wy[int(len(xrange)/2)]*1e6*np.sqrt(M2)) + " x " + str(wx[int(len(xrange)/2)]*1e6*np.sqrt(M2)) + " with M2 factor")
else:
    print(str(wx[int(len(xrange)/2)]*1e6) + " x " + str(wy[int(len(xrange)/2)]*1e6) + " at center")
    if M2 != 1:
        print(str(wx[int(len(xrange)/2)]*1e6*np.sqrt(M2)) + " x " + str(wy[int(len(xrange)/2)]*1e6*np.sqrt(M2)) + " with M2 factor")
print()
#phase = GB.Prop_EPhase(q_x,q_y,zoom,wavelength,Ei,w0)
phase = GB.Prop_EPhase_M2(q_x,q_y,zoom,wavelength,Ei,w0,M2=M2)
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
#sys.exit()
I_start=GB.GaussianBeamIntensity_SpotArray_2D(I0/1e4,wx[0],wy[0],w0,0,0)
print()
print("Intensity at start: ",I_start/1e12,"x10^12 W/cm2")
print("Change zi to find I at first lens.")

if save == 1:
    #Create the folder if it not already exists
    if not os.path.exists(path):
        os.makedirs(path)
    np.save(pulseParams['Location'],pulseParams)

wy = wy*np.sqrt(M2)
wx = wx*np.sqrt(M2)

if calcfocal == 1:
    den = 1
    #den = 10
    params = frefract.GetDefaultParams()
    X = params['X']*foc_dom_fac*radscl; Nx = params['Nx']
    y = np.linspace(-X/2, X/2, Nx, False)/1e6
    
    centerchange = 0
    if centerchange == 1:
        wmult = wy * wx
        wmult_min = np.argmin(wmult)
    else:
        if min(wy) < min(wx):
            wmult_min = np.argmin(wy)
        else:
            wmult_min = np.argmin(wx)
    
    #wz_center = wy[int(len(wy)/2)]
    #wy_center = wx[int(len(wx)/2)]
    print("center change: ",int(len(wy)/2), " to " ,wmult_min)
    print(xrange[int(len(wy)/2)]," to ",xrange[wmult_min])
    wz_center = wy[wmult_min]
    wy_center = wx[wmult_min]
    
    
    I = GB.IntensityFromSpotSizes_1D(wy_center,wz_center,y,I0/1e4,w0)

    plt.plot(y*1e6,I)
    plt.title("Intensity along beam axis")
    plt.xlabel(r'$\mathrm{Beam\,\,Axis\,\,[\mu m]}$')
    plt.ylabel(r'$\mathrm{Intensity \, [W/cm^2]}$')
    plt.grid(); plt.show()
    
    params['EI'] = E_ion
    H = ThrDim.IonFracFromIntensity_1D(I,params['EI'],fs_duration)
    H = H * den
    plt.plot(y*1e6,H)
    plt.title("Ion.Frac. along beam axis")
    plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
    #plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
    plt.grid(); plt.show()
    
    focal = Foc.Calc_Focus(H, y*1e6)
    thick = Foc.Calc_Square_Lens(den*1e17, focal)
    print(thick,"um equivalent thickness")
    
    guess = [thick, thick/6, den]
    fit = ThrDim.FitDataDoubleTanh(H, y*1e6, guess)

if calcdensity == 1:
    den = 0.3
    params = frefract.GetDefaultParams()
    params['Nz'] = len(xrange)
    params['Nx'] = int(params['Nx']/4)
    params['Ny'] = int(params['Ny']/4)
    params['Z'] = zi*1e6*2
    
    X = params['X']*radscl; Nx = params['Nx']
    Y = params['Y']; Ny = params['Ny']
    Z = params['Z']; Nz = params['Nz']

    y = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(-Y/8, Y/8, Ny, False)
    x = np.linspace(-Z/2, Z/2, Nz, False)
    
    wz=np.array(wy)*1e6
    wy=np.array(wx)*1e6
    I = GB.IntensityFromSpotSizes(wy/1e6,wz/1e6,x/1e6,y/1e6,z/1e6,I0/1e4,w0)
    ThrDim.ImageCut(I,x,y,z,0,0,0,1e-3,'(mm)','Intensity','W/cm^2',0)
    params['EI'] = E_ion
    H = ThrDim.IonFracFromIntensity(I,params['EI'],fs_duration)
    ThrDim.ImageCut(H,x,y,z,0,0,0,1e-3,'(mm)','Ion. Frac.','%',0)
    H = H * den
    #H = ThrDim.RobertRoll(ThrDim.RobertRoll(H))
    ThrDim.ImageCut(H,x,y,z,0,0,0,1e-3,'(mm)','Plasma Density','e17(cm^-3)',1)
    #frefract.TestDensity(H,params)
    
    if save == 1:
        savefolder = '/home/chris/Desktop/FourierPlots/ArVarySpot_Calc/'
        savefolder = savefolder + 'case_2/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        np.save(savefolder+'finalDensity.npy',H)
        np.save(savefolder+'params.npy',params)