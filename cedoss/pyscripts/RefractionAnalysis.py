#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:17:01 2017

Loads a 3D density profile from the output of FourierRefraction and analyzes
the 2D planes and variance cuts of the profile

@author: chris
"""

import sys
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim
import numpy as np

cuts=1
x_correct=1
z_correct=1

#size of window in micrometers
y_window = 100
z_window = 400

folder = '/home/chris/Desktop/FourierPlots/real_FACET_Refraction/'
directory = 'gasjet_den_propagation_1e20/'
path = folder+directory

nplot = np.load(path+'finalDensity.npy')
params = np.load(path+'params.npy').item()
X = params['X']
Nx = params['Nx']
Y = params['Y']
Ny = params['Ny']
Z = params['Z']
Nz = params['Nz']

den = ThrDim.RobertRoll(nplot)
y = np.linspace(-X/2, X/2, Nx, False)
z = np.linspace(-Y/2, Y/2, Ny, False)
x = np.linspace(-Z/2, Z/2, Nz, False)

ywidth=int(round(len(y)/(y[0]/(-.5*y_window))))
yi=int(round((len(y)-ywidth)/2))
yf=yi+ywidth

zwidth=int(round(len(z)/(z[0]/(-.5*z_window))))
zi=int(round((len(z)-zwidth)/2))
zf=zi+zwidth

den=den[:,yi:yf,zi:zf]
y=y[yi:yf]
z=z[zi:zf]

x_off=0
if x_correct == 1:
    denx = list(den[:,round(len(y)/2),round(len(z)/2)])
    loc = denx.index(max(denx))
    x_off = loc - round(len(denx)/2)
    print('Corrected x plane at '+str(x[loc])+' microns')

ThrDim.ImageCut(den,x,y,z,x_off,0,0,1e-3,
                '(mm)','Plasma Density','e17(cm^-3)',1)
if cuts == 1:
    zoom=1#1e-3
    x_axis=[i*zoom for i in x]
    y_axis=[i*zoom for i in y]
    z_axis=[i*zoom for i in z]
    x_step=x[1]-x[0]
    y_step=y[1]-y[0]
    z_step=z[1]-z[0]
    
    den_plane_zx=np.transpose(den[:,round(len(y)/2),:])
    label=['Density along gas jet (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in +/- x(microns)']
    ThrDim.VarianceCut(den_plane_zx,z_axis,x_off,3,1,x_step,label,True)
    
    z_off=0
    if x_correct == 1:
        denz = list(den[round(len(x)/2)+x_off,round(len(y)/2),:])
        zloc = denz.index(max(denz))
        z_off = zloc - round(len(denz)/2)
        print('Corrected z plane at '+str(z[zloc])+' microns')
    
    den_plane_yz=den[round(len(x)/2)+x_off,:,:]
    label=['Density along beam axis (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in y(microns)']
    ThrDim.VarianceCut(den_plane_yz,y_axis,z_off,5,10,z_step,label)
    ThrDim.VarianceCut(den_plane_yz,y_axis,z_off,5,-10,z_step,label)
    label=['Density along beam axis (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in +/- y(microns)']
    ThrDim.VarianceCut(den_plane_yz,y_axis,z_off,10,5,z_step,label,True)
    
    den_plane_yx=np.transpose(den[:,:,round(len(z)/2)+z_off])
    label=['Density along beam axis (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in x(microns)']
    ThrDim.VarianceCut(den_plane_yx,y_axis,x_off,5,2,x_step,label)
    ThrDim.VarianceCut(den_plane_yx,y_axis,x_off,5,-2,x_step,label)
    label=['Density along beam axis (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in +/- x(microns)']
    ThrDim.VarianceCut(den_plane_yx,y_axis,x_off,10,1,x_step,label,True)