# -*- coding: utf-8 -*-
"""

@author: chris
"""
import sys


import GaussianBeam
import numpy as np
import matplotlib.pyplot as plt
import ADK_Combined as adkD
import DensityDistributions as dendis
import ThreeDimensionAnalysis as ThrDim

#Flags for plotting different things
only_w = 0
cut_plots_H = 1
contour_plots = 1
den_plot = 1
cut_plots_den = 1

c=2.998e8
P=60e9
cm_m=100
wavelength = 785.3e-9
w0 = 5e-3

I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)
delt_t = 50e-15
chi = 15.426

L_z = 500e-6
L_r = 500e-6
nozzle_loc = 1000e-6
nozzle_den = 1e18

l_step = 5e-5 #1e-4
zoom=100

y_window = 100e-6
z_window = 400e-6
t_step = 1e-6

choice=1
if choice==0:
    #Focuses wz over a distance of 4m with 1 cyl. lens to .2mm
    q_y = GaussianBeam.Prop_Init_q(wavelength, w0, -3.5, 1)
    q_z = GaussianBeam.Prop_Init_q(wavelength, w0, -3.5, 1)
    GaussianBeam.Prop_CylindricalLens(q_y, q_z, 8)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,3.5,1e-4)
    GaussianBeam.Prop_CylindricalLens(q_z, q_y, 1)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,1,1e-4)

if choice==1:
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

if only_w == 1:
    #Print some diagnostics for minimum spot sizes
    GaussianBeam.Prop_SpotInfo(q_y,wavelength,'w0y','x(w0y)')
    GaussianBeam.Prop_SpotInfo(q_z,wavelength,'w0z','x(w0z)')
    
    E0=adkD.Elefield(I0)
    GaussianBeam.Prop_EPhase(q_y,q_z,zoom,wavelength,E0,w0)
    sys.exit()

#Create a 3D array for the intensity of this guassian beam as a function of
# both w(z) and r
yrange = ThrDim.BuildDomain(y_window,t_step)
zrange = ThrDim.BuildDomain(z_window,t_step)
I = ThrDim.IntensityFromSpotSizes(wy,wz,xrange,yrange,zrange,t_step,I0,w0)

#1-to-1 corresponance between intensity and ADK ionization
H = ThrDim.IonFracFromIntensity(I,chi,delt_t)

#Shift the origin to the waist of wy
waist_loc = xrange[round(len(xrange)/2)]
xrange=[i-waist_loc for i in xrange]

if cut_plots_H==1:
    #Make some plots for H vs y with small variances in z
    yrange_cut=yrange[round(len(yrange)*1/4):round(len(yrange)*3/4)]
    axis=[i*1e6 for i in yrange_cut]
    
    H_plane=H[round(len(xrange)/2),round(len(yrange)*1/4):round(len(yrange)*3/4),:]
    label=['Ion. Frac. through waist (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni/n0',
           'Offset in z(microns)']
    ThrDim.VarianceCut(H_plane,axis,0,5,20,t_step*1e6,label)
    
    H_plane=H[:,round(len(yrange)*1/4):round(len(yrange)*3/4),round(len(zrange)/2)]
    H_plane=np.transpose(H_plane)
    label=['Ion. Frac. through waist (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni/n0',
           'Offset in x(microns)']
    ThrDim.VarianceCut(H_plane,axis,0,3,2,l_step*1e6,label)

if contour_plots==1:
    #Countour plots for I and H on interesting planes for our 3D plasma
    ThrDim.ImageCut(I,xrange,yrange,zrange,0,0,0,1000,'(mm)','Intensity','(W/cm^2)')
    ThrDim.ImageCut(H,xrange,yrange,zrange,0,0,0,1000,'(mm)','Ionization_Frac','(ni/n0)',1)

if den_plot + cut_plots_den > 0:
    den=np.empty([len(xrange),len(yrange),len(zrange)])
    for x in np.arange(len(xrange)):
        for y in np.arange(len(yrange)):
            for z in np.arange(len(zrange)):
                den[x][y][z]=H[x][y][z]*dendis.SimpleGasJet(nozzle_den,
                   L_z,L_r,nozzle_loc,xrange[x],yrange[y],zrange[z])
if den_plot == 1:
    ThrDim.ImageCut(den,xrange,yrange,zrange,0,0,0,1000,'(mm)','Plasma Density','(cm^-3)',1)
    
if cut_plots_den == 1:
    y_axis=[i*1e6 for i in yrange]
    z_axis=[i*1e6 for i in zrange]
    
    den_plane_zx=np.transpose(den[:,round(len(yrange)/2),:])
    label=['Density along gas jet (Vary Laser Distance)',
           'Radius from axis (microns)',
           'ni',
           'Offset in x(microns)']
    ThrDim.VarianceCut(den_plane_zx,z_axis,0,2,2,l_step*1e6,label)
    
    den_plane_yz=den[round(len(xrange)/2),:,:]
    label=['Density along beam axis (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni',
           'Offset in x(microns)']
    ThrDim.VarianceCut(den_plane_yz,y_axis,0,5,20,t_step*1e6,label)
    ThrDim.VarianceCut(den_plane_yz,y_axis,0,5,-20,t_step*1e6,label)
    label=['Density along beam axis (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni',
           'Offset in +/- x(microns)']
    ThrDim.VarianceCut(den_plane_yz,y_axis,0,11,10,t_step*1e6,label,True)
    
    z_cut=list(den_plane_zx[:,round(len(xrange)/2)])
    offset=round(z_axis[z_cut.index(max(z_cut))])
    label=['Density along beam axis z='+str(offset)+' (Vary Jet Distance)',
           'Radius from axis (microns)',
           'ni',
           'Offset in +/- x(microns)']
    ThrDim.VarianceCut(den_plane_yz,y_axis,int(offset),11,10,t_step*1e6,label,True)