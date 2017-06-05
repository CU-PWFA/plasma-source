# -*- coding: utf-8 -*-
"""

@author: chris
"""

import GaussianBeam
import numpy as np
import matplotlib.pyplot as plt
import ADK_Combined as adk
import DensityDistributions as dendis
import sys

c=2.998e8
P=12e9
cm_m=100
wavelength = 785.3e-9
w0 = 5e-3

I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)
delt_t = 100e-15
chi = adk.Get_chi_H()
choice=0
if choice==0:
    q_y = GaussianBeam.Prop_Init_q(wavelength, w0, -3.5, 1)
    q_z = GaussianBeam.Prop_Init_q(wavelength, w0, -3.5, 1)
    GaussianBeam.Prop_CylindricalLens(q_y, q_z, 8)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,3.5,1e-4)
    GaussianBeam.Prop_CylindricalLens(q_z, q_y, 1)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,1,1e-4)

if choice==1:
    q_y = GaussianBeam.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GaussianBeam.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GaussianBeam.Prop_CylindricalLens(q_y, q_z, .8)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,.35,1e-4)
    GaussianBeam.Prop_CylindricalLens(q_y, q_z, -.105)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,.15,1e-4)
    
    GaussianBeam.Prop_CylindricalLens(q_z, q_y, 1)
    GaussianBeam.Prop_Cylindrical_FreeSpace(q_y,q_z,1,1e-4)

xrange_tot=GaussianBeam.Prop_GetRange(q_y)
wtoty=GaussianBeam.Prop_SpotList(q_y,wavelength)
wtotz=GaussianBeam.Prop_SpotList(q_z,wavelength)

print(min(wtoty))
print(min(wtotz))
print(xrange_tot[wtoty.index(min(wtoty))])
print(xrange_tot[wtotz.index(min(wtotz))])

#Plots spot size
plt.plot(xrange_tot,wtoty)
plt.title("Spot size from q parameter:  blue transverse, orange parallel")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.plot(xrange_tot,wtotz)
plt.show()

waist=wtoty.index(min(wtoty))
zoom=100
wy=wtoty[(waist-zoom):(waist+zoom)]
wz=wtotz[(waist-zoom):(waist+zoom)]

xrange=xrange_tot[(waist-zoom):(waist+zoom)]

#Plots spot size zoomed in around waist
plt.plot(xrange,wy)
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.plot(xrange,wz)
plt.show()


"""REMOVE ONCE CHOICE 1 IS DEEMED COOL ENOUGH TO ANALYZE"""
if choice==1:
    sys.exit()

#Create a 3D array for the intensity of this guassian beam as a function of
# both w(z) and r
rrange=np.arange(-200e-6,200e-6,.5e-5)
I=np.empty([len(xrange),len(rrange),len(rrange)])
i_counter=np.arange(len(xrange))
j_counter=np.arange(len(rrange))
k_counter=np.arange(len(rrange))
for i in i_counter:
    for j in j_counter:
        for k in k_counter:
            I[i][j][k]=GaussianBeam.GaussianBeamIntensity_SpotArray_2D(I0,wy[i],wz[i],w0,rrange[j],rrange[k])

H=np.empty([len(xrange),len(rrange),len(rrange)])
for i in i_counter:
    for j in j_counter:
        H[i][j] = adk.TotalProb(adk.Elefield(I[i][j]),delt_t,chi,1)

waist_loc = xrange[round(len(xrange)/2)]
xrange=[i-waist_loc for i in xrange]

cut_plots = 0
if cut_plots==1:
    H_cut=H[:,round(len(rrange)*1/4):round(len(rrange)*3/4),:]
    rrange_cut=rrange[round(len(rrange)*1/4):round(len(rrange)*3/4)]
    H_center    = H_cut[round(len(xrange)/2),:,round(len(rrange)/2)]
    H_center_dz = H_cut[round(len(xrange)/2),:,round(len(rrange)/2)+50]
    H_center_dz2 = H_cut[round(len(xrange)/2),:,round(len(rrange)/2)+100]
    rrange_zoomed=[i*1e6 for i in rrange_cut]
    plt.plot(rrange_zoomed,H_center)
    plt.title("Ion. Frac. through waist")
    plt.ylabel("ni/n0")
    plt.xlabel("Radius from axis (microns)")
    plt.plot(rrange_zoomed,H_center_dz) #Orange is 50microns in z direction
    plt.plot(rrange_zoomed,H_center_dz2) #Green is 100microns in z direction
    plt.show()

contour_plots = 1
if contour_plots==1:
    rrange_z=[i*1000 for i in rrange]
    #Plot intensity at x=waist
    I_atwaist=np.transpose(I[round(len(xrange)/2)])
    plt.imshow(I_atwaist, interpolation="none", origin="lower",
                 extent=[rrange_z[0],rrange_z[-1],rrange_z[0],rrange_z[-1]], aspect='1')
    CB=plt.colorbar()
    CB.set_label('I (I0 = 3e10')
    plt.ylabel('z (mm)')
    plt.xlabel('y (mm)')
    plt.title('Intensity; x=waist')
    plt.show()

    H_atwaist=np.transpose(H[round(len(xrange)/2)])
    plt.imshow(H_atwaist, interpolation="none", origin="lower",
                 extent=[rrange_z[0],rrange_z[-1],rrange_z[0],rrange_z[-1]], aspect='1')
    CB=plt.colorbar()
    CB.set_label('ni/n0')
    plt.ylabel('z (mm)')
    plt.xlabel('y (mm)')
    plt.title('ADK Ion.; x=waist')
    plt.show()

    #Plot intensity at y=0
    I_yzero=I[:,round(len(rrange)/2),:]
    plt.imshow(I_yzero, interpolation="none", origin="lower",
               extent=[rrange_z[0],rrange_z[-1],xrange[0],xrange[-1]],aspect='10')
    CB=plt.colorbar()
    CB.set_label('I (I0 = 3e10')
    plt.ylabel('x (m)')
    plt.xlabel('z (mm)')
    plt.title('Intensity; y=0')
    plt.show()
    
    H_yzero=H[:,round(len(rrange)/2),:]
    plt.imshow(H_yzero, interpolation="none", origin="lower",
               extent=[rrange_z[0],rrange_z[-1],xrange[0],xrange[-1]], aspect='10')
    CB=plt.colorbar()
    CB.set_label('ni/n0')
    plt.ylabel('x (m)')
    plt.xlabel('z (mm)')
    plt.title('ADK Ion.; y=0')
    plt.show()

    #Plot intensity at z=0
    I_zzero=I[:,:,round(len(rrange)/2)]
    plt.imshow(I_zzero, interpolation="none", origin="lower",
               extent=[rrange_z[0],rrange_z[-1],xrange[0],xrange[-1]],aspect='10')
    CB=plt.colorbar()
    CB.set_label('I (I0 = 3e10')
    plt.ylabel('x (m)')
    plt.xlabel('y (mm)')
    plt.title('Intensity; z=0')
    plt.show()

    H_zzero=H[:,:,round(len(rrange)/2)]
    plt.imshow(H_zzero, interpolation="none", origin="lower",
               extent=[rrange_z[0],rrange_z[-1],xrange[0],xrange[-1]], aspect='10')
    CB=plt.colorbar()
    CB.set_label('ni/n0')
    plt.ylabel('x (m)')
    plt.xlabel('y (mm)')
    plt.title('ADK Ion.; z=0')
    plt.show()
    
den_plot = 1
if den_plot == 1:
    den=np.empty([len(xrange),len(rrange),len(rrange)])
    for x in i_counter:
        for y in j_counter:
            for z in k_counter:
                den[x][y][z]=H[x][y][z]*dendis.SimpleGasJet(1e18,500e-6,500e-6,200e-6,xrange[x],rrange[y],rrange[z])

    rrange_z=[i*1000 for i in rrange]

    den_atwaist=np.transpose(den[round(len(xrange)/2)])
    plt.imshow(den_atwaist, interpolation="none", origin="lower",
                 extent=[rrange_z[0],rrange_z[-1],rrange_z[0],rrange_z[-1]], aspect='1')
    CB=plt.colorbar()
    CB.set_label('ni')
    plt.ylabel('z (mm)')
    plt.xlabel('y (mm)')
    plt.title('Plasma Density; x=waist')
    plt.show()

    den_yzero=den[:,round(len(rrange)/2),:]
    plt.imshow(den_yzero, interpolation="none", origin="lower",
               extent=[rrange_z[0],rrange_z[-1],xrange[0],xrange[-1]], aspect='10')
    CB=plt.colorbar()
    CB.set_label('ni')
    plt.ylabel('x (m)')
    plt.xlabel('z (mm)')
    plt.title('Plasma Density; y=0')
    plt.show()

    den_zzero=den[:,:,round(len(rrange)/2)]
    plt.imshow(den_zzero, interpolation="none", origin="lower",
               extent=[rrange_z[0],rrange_z[-1],xrange[0],xrange[-1]], aspect='10')
    CB=plt.colorbar()
    CB.set_label('ni')
    plt.ylabel('x (m)')
    plt.xlabel('y (mm)')
    plt.title('ADK Ion.; z=0')
    plt.show()