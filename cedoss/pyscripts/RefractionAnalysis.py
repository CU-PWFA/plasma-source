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
max_corrector=0

#size of window in micrometers
y_window = 100
z_window = 400

folder = '/home/chris/Desktop/FourierPlots/CompactOptics_DoubleJet/'
#folder = '/home/chris/Desktop/FourierPlots/CompactOptics/'
#folder = '/home/chris/Desktop/FourierPlots/real_FACET_Refraction/'
directory = 'gasjet_den_propagation_1e18_Lz500/'
path = folder+directory

nplot = np.load(path+'finalDensity.npy')
params = np.load(path+'params.npy').item()

X = params['X']; Nx = params['Nx']
Y = params['Y']; Ny = params['Ny']
Z = params['Z']; Nz = params['Nz']

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

x_off=0; y_off=0; z_off=0
if max_corrector == 1:
    maximum = ThrDim.GetMaximumOffset(den)
    x_off=maximum[0]; y_off=maximum[1]; z_off=maximum[2]
    print('-Corrected beam position to maximum of '+str(maximum[3])+' e17cm^-3')
    print('-Corrected x plane at '+str(x[x_off + round(len(den[:,0,0])/2)])+' microns')
    print('-Corrected y plane at '+str(y[y_off + round(len(den[0,:,0])/2)])+' microns')
    print('-Corrected z plane at '+str(z[z_off + round(len(den[0,0,:])/2)])+' microns')

ThrDim.ImageCut(den,x,y,z,x_off,y_off,z_off,1e-3,
                '(mm)','Plasma Density','e17(cm^-3)',1)

if cuts == 1:
    x_step=x[1]-x[0]
    y_step=y[1]-y[0]
    z_step=z[1]-z[0]
    
    den_plane_zx=np.transpose(den[:,round(len(y)/2),:])
    label=['Density along gas jet (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in +/- x(microns)']
    ThrDim.VarianceCut(den_plane_zx,z,x_off,3,1,x_step,label,True)
    
    den_plane_yz=den[round(len(x)/2)+x_off,:,:]
    label=['Density along beam axis (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in z(microns)']
    ThrDim.VarianceCut(den_plane_yz,y,z_off,5,10,z_step,label)
    ThrDim.VarianceCut(den_plane_yz,y,z_off,5,-10,z_step,label)
    label=['Density along beam axis (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in +/- z(microns)']
    ThrDim.VarianceCut(den_plane_yz,y,z_off,10,5,z_step,label,True)
    
    den_plane_yx=np.transpose(den[:,:,round(len(z)/2)+z_off])
    label=['Density along beam axis (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in x(microns)']
    ThrDim.VarianceCut(den_plane_yx,y,x_off,5,2,x_step,label)
    ThrDim.VarianceCut(den_plane_yx,y,x_off,5,-2,x_step,label)
    label=['Density along beam axis (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni (e17 cm^-3)',
           'Offset in +/- x(microns)']
    ThrDim.VarianceCut(den_plane_yx,y,x_off,4,1,x_step,label,True)