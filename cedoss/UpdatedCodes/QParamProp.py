#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:06:06 2023

Dedicated Script for testing optical setup and looking for best way to getthe desired spot sizes from cylindrical lens

This is a cleaned up version in June of 2023.  Includes M^2 functoinality

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

#Define some constants
c=2.998e8               #Speed of light constant
epsilon0 = 8.854e-12    #Dielectric constant

#For a simple analysis of the ionization process, define the ionization energy
#Note, this does not include refraction.  This is only for calculating ADK ionization
# rate from the laser intensity profile.
E_ion = 15.4    #Ionization Energy:  13.6 for H, 15.8 for Ar, 24.6 for He 27.6 for Ar++
den = 1         #Density of the uniform gas, in 10e17 cm-3.  

#Flags for the script.  True turns on, False turns off
save = False            #Save electric field profile at zi to the filename below.  Also, if calcdensity is True, saves the plasma density profile as a npy
calcdensity = False     #SLOW   Calculates the resulting 3D plasma density w/out refraction
calcfocal = True        #Takes the profile at the focus and calculates the 1D plasma lens thickness and focal length
offsetlooper = False    #Generates a quick scan of the lens thickness vs ebeam horizontal offset
reverse = False         #For cylindrical focus, flips which axis is the ebeam axis
single=False            #For plots that have cylindrical symmetery, only plot one transverse component

#If we are saving, then save to path+filename
path = '/home/chris/Desktop/DataLoads/PulseFilesNp/'
filename = 'pulseParams_example.npy'

#Effectively chooses which inputs we are using
#These blocks define both the laser parameters and the optical setup
choice=2        #Set to one of the setups below
    
if choice == 1:  #1 Spherical Lens, similar to old OAP at FACET-II
    setupTitle = "SingleSpherical"
    
    #Laser Parameters
    wavelength = 800e-9             #Laser Wavelength in m
    w0 = 13.6e-3                    #Initial Laser Spot Size in m.  13.6 mm roughly corresponds to 30 mm flattop
    M2 = 2.0                        #M^2 Value.  M2=1 corresponds to an ideal Gaussian Beam
    fs_duration = 70e-15            #Temporal pulse width in s
    laser_J = 10e-3                 #Energy of the laser in J
    laser_P = laser_J / fs_duration #Power of the laser in W
    
    #Parameters on the window size at the focus.  "zi" is especially important here
    zi = 1e-2       #Distance from focus as which to start the zoomed window.  ALSO, sets the location where the E-field is saved
    radscl = 3.0    #Relative size of the transverse window size compared to 30 um (yeah, I know.  kinda arbitrary and bad)
    l_step = 5e-5   #Longitudinal step size in m, make sure all lengths are divisible by l_step
    zoom = int(round(zi/l_step))
    
    #Define the optic and subsequent free space propagation
    f0 = 0.646#0.646
    Lmax = 1.000#1.0
    
    #Initialize the q parameters of each 
    initial_position = -0.1 #Position with respect to the waist to initialize laser at, in m
    initial_index = 1.0 #Refractive index at the start of the propagation
    q_x = GB.Prop_Init_q(wavelength, w0, initial_position, initial_index)
    q_y = GB.Prop_Init_q(wavelength, w0, initial_position, initial_index)
    
    #Propagate laser from the start to the optic at z=0.
    #This function propagates both transvese components
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,-1*initial_position,l_step)
    
    #Propagate both transverse components through the spherical lens.
    #NOTE: transfer matrices use radius of curvature, so I multiply focal length by 2 here
    GB.Prop_SphericalLens(q_y,q_x,2*f0)

    #Lastly, propagate the beam after the lens through free space
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,Lmax,l_step)

if choice == 2:  #Three crossed cylidnrical lenses for an unexpanded beam
    setupTitle = "CrossedCylindrical"
    
    #Laser Parameters
    wavelength = 800e-9                 #Laser Wavelength in m
    w0 = 4.77e-3                        # #Initial Laser Spot Size in m.  4.77 mm roughly corresponds to 10 mm flattop
    M2 = 1.0                            #M^2 Value.  M2=1 corresponds to an ideal Gaussian Beam
    fs_duration = 70e-15                #Temporal pulse width in s
    laser_J = 10e-3                     #Energy of the laser in J
    laser_P = laser_J / fs_duration     #Power of the laser in W
    
    #Parameters on the window size at the focus.  "zi" is especially important here
    zi = 6.0e-2     #Distance from focus as which to start the zoomed window.  ALSO, sets the location where the E-field is saved
    radscl = 2.0    #Relative size of the transverse window size compared to 30 um (yeah, I know.  kinda arbitrary and bad)
    l_step = 5e-5   #Longitudinal step size in m, make sure all lengths are divisible by l_step
    zoom = int(round(zi/l_step))
    
    #Lens 1 for horizontal, at z=0
    f0 = 0.20

    #Lens 2 for horizontal, at z=109 mm
    f2 = -0.1000 #m
    L2 = 0.109

    #Lens 3 for vertical, at z=417.25 mm
    f1 = 0.700 #m
    L1 = 0.41725

    #Free space distance after Lens 3
    Lmax = 1.7

    #Initialize the q parameters of each 
    initial_position = -0.1 #Position with respect to the waist to initialize laser at, in m
    initial_index = 1.0 #Refractive index at the start of the propagation
    q_x = GB.Prop_Init_q(wavelength, w0, initial_position, initial_index)
    q_y = GB.Prop_Init_q(wavelength, w0, initial_position, initial_index)
    
    #Free space to Lens 1
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,-1*initial_position,l_step)
    
    #Thin Lens 1 modifying q_x (with radius of curvature 2*f0)
    GB.Prop_CylindricalLens(q_x,q_y,2*f0)
    
    #Free space to Lens 2
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L2,l_step) 
    
    #Lens 2 modifying q_x
    GB.Prop_CylindricalLens(q_x,q_y,2*f2)

    #Free space to Lens 3
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y,L1-L2,l_step) 
    
    #Lens 3 modifying q_y this time
    GB.Prop_CylindricalLens(q_y,q_x,2*f1)

    #Free space through the focus and beyond
    GB.Prop_Cylindrical_FreeSpace(q_x,q_y, Lmax,l_step)

#Calculate the initial power and electric field strength at the start
# using the laser parameters
I0 = 2*laser_P/np.pi/np.power(w0,2)
Ei = np.sqrt(2*I0/c/initial_index/epsilon0)*1e-9

#Get the total domain of spot sizes
zrange_tot=GB.Prop_GetRange(q_x)
wtotx=GB.Prop_SpotList(q_x,wavelength)
wtoty=GB.Prop_SpotList(q_y,wavelength)

#Plots spot size
fig, ax1 = plt.subplots(figsize=(14,4))
#plt.rcParams.update({'font.size': 12})
plt.subplot(121)
if reverse:
    if not single:
        plt.plot(zrange_tot,wtoty,c='orange',ls='solid',label="Transverse Waist Spot Size")
        if M2 != 1:
            plt.plot(zrange_tot,np.array(wtoty)*np.sqrt(M2),c='orange',ls='--',label="Transverse with M2")
plt.plot(zrange_tot,wtotx,c='blue',ls='solid',label="Longitudinal Waist Spot Size")#"Narrow Waist Spot Size")
if M2 != 1:
    plt.plot(zrange_tot,np.array(wtotx)*np.sqrt(M2),c='blue',ls='--',label="Longitudinal with M2")#"Narrow M2")
plt.title("Spot size from q parameter")
plt.ylabel(r'Spot Size [$mm$]')
plt.xlabel(r'Laser Axis [$m$]')
if not reverse:
    if not single:
        plt.plot(zrange_tot,wtoty,c='orange',ls='solid',label="Transverse Waist Spot Size")
        if M2 != 1:
            plt.plot(zrange_tot,np.array(wtoty)*np.sqrt(M2),c='orange',ls='--',label="Transverse with M2")
plt.legend()
plt.grid()#; plt.show()

#Zoom in to a scale of the rayleigh length of wy
waist=wtotx.index(min(wtotx))
wx=np.array(wtotx[(waist-zoom):(waist+zoom)])
wy=np.array(wtoty[(waist-zoom):(waist+zoom)])
xrange=zrange_tot[(waist-zoom):(waist+zoom)]

if reverse:
    q_d = q_y
    q_y = q_x
    q_x = q_d
    wd = np.copy(wy)
    wy = np.copy(wx)
    wx = wd

#Plots spot size zoomed in around waist
plt.subplot(122)
plt.plot(xrange,wx*1e6,c='blue',ls='solid',label="Longitudinal Waist Spot Size")#"Beam Axis Waist (x)")
if M2 != 1:
    plt.plot(xrange,wx*np.sqrt(M2)*1e6,c='blue',ls='--',label="Longitudinal with M2")#"Beam Axis w/ M2 (x)")
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size [$\mu m$]")
plt.xlabel("Laser Axis [$m$]")
if not single:
    plt.plot(xrange,wy*1e6,c='orange',ls='solid',label="Transverse Waist Spot Size")
    if M2 != 1:
        plt.plot(xrange,wy*np.sqrt(M2)*1e6,c='orange',ls='--',label="Transverse with M2")
plt.legend()
plt.grid(); plt.show()

print(zoom)
print(zrange_tot[waist-zoom])
print(xrange[0])

#Print some diagnostics for minimum spot sizes
GB.Prop_SpotInfo(q_x,wavelength,'w0x','z(w0x)')
GB.Prop_SpotInfo(q_y,wavelength,'w0y','z(w0y)')
if not reverse:
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
#For quick reference; phase = [E0,wx,wy,Px,Py,phi,zi]

scl = 1e6 #Robert's code uses micrometers
pulseParams = {'Description' : setupTitle,
               'Location' : path+filename,
               'power' : laser_P,
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

I_start=GB.GaussianBeamIntensity_SpotArray_2D(I0/1e4,wx[0],wy[0],w0,0,0)
print()
print("Intensity at start: ",I_start/1e12,"x10^12 W/cm2")
print("Change zi to find I at first lens.")

if save:
    #Create the folder if it not already exists
    if not os.path.exists(path):
        os.makedirs(path)
    np.save(pulseParams['Location'],pulseParams)

wy = wy*np.sqrt(M2)
wx = wx*np.sqrt(M2)

#Calculate the focal length by finding the 1D plasma density along e-beam axis
if calcfocal:
    params = frefract.GetDefaultParams()
    X = params['X']*radscl; Nx = params['Nx']
    y = np.linspace(-X/2, X/2, Nx, False)/1e6
    
    centerchange = False
    if centerchange:
        wmult = wy * wx
        wmult_min = np.argmin(wmult)
    else:
        if min(wy) < min(wx):
            wmult_min = np.argmin(wy)
        else:
            wmult_min = np.argmin(wx)
    
    if offsetlooper:
        l_arr = np.zeros(len(wy))
        o_arr = (np.array(range(len(l_arr)))-int(0.5*len(l_arr)))*(xrange[1]-xrange[0])*1e3
        for i in range(len(l_arr)):
            wz_current = wy[i]
            wy_current = wx[i]
            I = GB.IntensityFromSpotSizes_1D(wy_current,wz_current,y,I0/1e4,w0)            
            params['EI'] = E_ion
            H = ThrDim.IonFracFromIntensity_1D(I,params['EI'],fs_duration)*den
            focal = Foc.Calc_Focus(H, y*1e6)
            l_arr[i] = Foc.Calc_Square_Lens(den*1e17, focal)
        plt.title("Longitudinal Plasma Thickness vs Horizontal Offset")
        plt.xlabel("Horizontal Laser Focus Offset "+r'$(mm)$')
        plt.ylabel("Plasma Thickness in Electron Beam Axis "+r'$(\mu m)$')
        plt.plot(o_arr,l_arr)
        plt.grid(); plt.show()
        print("Maximum Thickness: ",max(l_arr))
    
    #wz_center = wy[int(len(wy)/2)]
    #wy_center = wx[int(len(wx)/2)]
    print("center change: ",int(len(wy)/2), " to " ,wmult_min)
    print(xrange[int(len(wy)/2)]," to ",xrange[wmult_min])
    wz_center = wy[wmult_min]
    wy_center = wx[wmult_min]
    
    
    I = GB.IntensityFromSpotSizes_1D(wy_center,wz_center,y,I0/1e4,w0)
    print(wy_center,wz_center,w0,I0/1e4)

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

#Calculate the 3D plasma density assuming no refraction
if calcdensity:
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
    ThrDim.ImageCut(H,x,y,z,0,0,0,1e-3,'(mm)','Ionization Fraction','',0)
    H = H * den
    #H = ThrDim.RobertRoll(ThrDim.RobertRoll(H))
    ThrDim.ImageCut(H,x,y,z,0,0,0,1e-3,'(mm)','Plasma Density','e17(cm^-3)',1)
    #frefract.TestDensity(H,params)
    
    if save:
        savefolder = path#'/home/chris/Desktop/FourierPlots/ArVarySpot_Calc/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        np.save(savefolder+'finalDensity.npy',H)
        np.save(savefolder+'params.npy',params)