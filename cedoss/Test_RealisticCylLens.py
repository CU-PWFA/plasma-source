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

sys.path.insert(0, "../python")

from ionization import ionization
from ionization import adk

#Flags for plotting different things
cut_plots_H = 0
contour_plots = 0
den_plot = 0
cut_plots_den = 1

c=2.998e8
P=15e9
cm_m=100
wavelength = 785.3e-9
w0 = 5e-3

I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)
delt_t = 100e-15
chi = adkD.Get_chi_H()

L_z = 500e-6
L_r = 500e-6
nozzle_loc = 200e-6
nozzle_den = 1e18

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
zoom=50
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
yrange=np.arange(-50e-6,50e-6,1e-6)
zrange=np.arange(-200e-6,200e-6,1e-6)

I=np.empty([len(xrange),len(yrange),len(zrange)])
i_counter=np.arange(len(xrange))
j_counter=np.arange(len(yrange))
k_counter=np.arange(len(zrange))
for i in i_counter:
    for j in j_counter:
        for k in k_counter:
            I[i][j][k]=GaussianBeam.GaussianBeamIntensity_SpotArray_2D(I0,wy[i],wz[i],w0,yrange[j],zrange[k])

H=np.empty([len(xrange),len(yrange),len(zrange)])\

square_wave=0
if square_wave == 1:
    for i in i_counter:
        for j in j_counter:
            H[i][j] = adkD.TotalProb(adkD.Elefield(I[i][j]),delt_t,chi,1)
else:
    for i in i_counter:
        for j in j_counter:
            E = ionization.field_from_intensity(I[i][j]/1e14)
            H[i][j] = adk.gaussian_frac(chi,E,delt_t/1e-15,1)

waist_loc = xrange[round(len(xrange)/2)]
xrange=[i-waist_loc for i in xrange]

if cut_plots_H==1:
    H_cut=H[:,round(len(yrange)*1/4):round(len(yrange)*3/4),:]
    yrange_cut=yrange[round(len(yrange)*1/4):round(len(yrange)*3/4)]
    H_center    = H_cut[round(len(xrange)/2),:,round(len(zrange)/2)]
    H_center_dz = H_cut[round(len(xrange)/2),:,round(len(zrange)/2)+20]
    H_center_dz2 = H_cut[round(len(xrange)/2),:,round(len(zrange)/2)+40]
    H_center_dz3 = H_cut[round(len(xrange)/2),:,round(len(zrange)/2)+60]
    rrange_zoomed=[i*1e6 for i in yrange_cut]
    plt.plot(rrange_zoomed,H_center)
    plt.title("Ion. Frac. through waist")
    plt.ylabel("ni/n0")
    plt.xlabel("Radius from axis (microns)")
    plt.plot(rrange_zoomed,H_center_dz) #Orange is 20 microns in z direction
    plt.plot(rrange_zoomed,H_center_dz2) #Green is 40 microns in z direction
    plt.plot(rrange_zoomed,H_center_dz3) #Red is 60 microns
    plt.show()

if contour_plots==1:
    xrange_z=[i*1000 for i in xrange]
    yrange_z=[i*1000 for i in yrange]
    zrange_z=[i*1000 for i in zrange]
    #Plot intensity at x=waist
    I_atwaist=np.transpose(I[round(len(xrange)/2)])
    plt.imshow(I_atwaist, interpolation="none", origin="lower",
                 extent=[yrange_z[0],yrange_z[-1],zrange_z[0],zrange_z[-1]], aspect='.5')
    CB=plt.colorbar()
    CB.set_label('I (I0 = 3.8e10')
    plt.ylabel('z (mm)')
    plt.xlabel('y (mm)')
    plt.title('Intensity; x=waist')
    plt.show()

    H_atwaist=np.transpose(H[round(len(xrange)/2)])
    plt.imshow(H_atwaist, interpolation="none", origin="lower",
                 extent=[yrange_z[0],yrange_z[-1],zrange_z[0],zrange_z[-1]], aspect='.5')
    CB=plt.colorbar()
    CB.set_label('ni/n0')
    plt.ylabel('z (mm)')
    plt.xlabel('y (mm)')
    plt.title('ADK Ion.; x=waist')
    plt.show()

    #Plot intensity at y=0
    I_yzero=I[:,round(len(yrange)/2),:]
    plt.imshow(I_yzero, interpolation="none", origin="lower",
               extent=[zrange_z[0],zrange_z[-1],xrange_z[0],xrange_z[-1]],aspect='.02')
    CB=plt.colorbar()
    CB.set_label('I (I0 = 3.8e10')
    plt.ylabel('x (mm)')
    plt.xlabel('z (mm)')
    plt.title('Intensity; y=0')
    plt.show()
    
    H_yzero=H[:,round(len(yrange)/2),:]
    plt.imshow(H_yzero, interpolation="none", origin="lower",
               extent=[zrange_z[0],zrange_z[-1],xrange_z[0],xrange_z[-1]], aspect='.02')
    CB=plt.colorbar()
    CB.set_label('ni/n0')
    plt.ylabel('x (mm)')
    plt.xlabel('z (mm)')
    plt.title('ADK Ion.; y=0')
    plt.show()

    #Plot intensity at z=0
    I_zzero=I[:,:,round(len(zrange)/2)]
    plt.imshow(I_zzero, interpolation="none", origin="lower",
               extent=[yrange_z[0],yrange_z[-1],xrange_z[0],xrange_z[-1]],aspect='.01')
    CB=plt.colorbar()
    CB.set_label('I (I0 = 3.8e10')
    plt.ylabel('x (mm)')
    plt.xlabel('y (mm)')
    plt.title('Intensity; z=0')
    plt.show()

    H_zzero=H[:,:,round(len(zrange)/2)]
    plt.imshow(H_zzero, interpolation="none", origin="lower",
               extent=[yrange_z[0],yrange_z[-1],xrange_z[0],xrange_z[-1]], aspect='.01')
    CB=plt.colorbar()
    CB.set_label('ni/n0')
    plt.ylabel('x (mm)')
    plt.xlabel('y (mm)')
    plt.title('ADK Ion.; z=0')
    plt.show()
    
if den_plot + cut_plots_den > 0:
    den=np.empty([len(xrange),len(yrange),len(zrange)])
    for x in i_counter:
        for y in j_counter:
            for z in k_counter:
                den[x][y][z]=H[x][y][z]*dendis.SimpleGasJet(nozzle_den,L_z,L_r,nozzle_loc,xrange[x],yrange[y],zrange[z])
if den_plot == 1:
    
    xrange_z=[i*1000 for i in xrange]
    yrange_z=[i*1000 for i in yrange]
    zrange_z=[i*1000 for i in zrange]

    den_atwaist=np.transpose(den[round(len(xrange)/2)])
    plt.imshow(den_atwaist, interpolation="none", origin="lower",
                 extent=[yrange_z[0],yrange_z[-1],zrange_z[0],zrange_z[-1]], aspect='.5')
    CB=plt.colorbar()
    CB.set_label('n')
    plt.ylabel('z (mm)')
    plt.xlabel('y (mm)')
    plt.title('Plasma Density; x=waist')
    plt.show()

    den_yzero=den[:,round(len(yrange)/2),:]
    plt.imshow(den_yzero, interpolation="none", origin="lower",
               extent=[zrange_z[0],zrange_z[-1],xrange_z[0],xrange_z[-1]], aspect='.02')
    CB=plt.colorbar()
    CB.set_label('n')
    plt.ylabel('x (mm)')
    plt.xlabel('z (mm)')
    plt.title('Plasma Density; y=0')
    plt.show()

    den_zzero=den[:,:,round(len(zrange)/2)]
    plt.imshow(den_zzero, interpolation="none", origin="lower",
               extent=[yrange_z[0],yrange_z[-1],xrange_z[0],xrange_z[-1]], aspect='.01')
    CB=plt.colorbar()
    CB.set_label('n')
    plt.ylabel('x (mm)')
    plt.xlabel('y (mm)')
    plt.title('Plasma Density; z=0')
    plt.show()
    
if cut_plots_den == 1:
    xrange_z=[i*1000 for i in xrange]
    yrange_z=[i*1000 for i in yrange]
    zrange_z=[i*1000 for i in zrange]
    
    n_vs_z   = den[round(len(xrange_z)/2),round(len(yrange_z)/2),:]
    n_vs_y   = den[round(len(xrange_z)/2),:,round(len(zrange_z)/2)]
    n_vs_y1  = den[round(len(xrange_z)/2),:,round(len(zrange_z)/2)+20]
    n_vs_y2  = den[round(len(xrange_z)/2),:,round(len(zrange_z)/2)+40]
    n_vs_y3  = den[round(len(xrange_z)/2),:,round(len(zrange_z)/2)+60]
    n_vs_y2d = den[round(len(xrange_z)/2),:,round(len(zrange_z)/2)-40]
    
    plt.plot(zrange_z,n_vs_z)
    plt.title("Plasma density along gas jet")
    plt.ylabel("n")
    plt.xlabel("Radius from laser axis (mm)")
    plt.show()
    
    plt.plot(yrange_z,n_vs_y)
    plt.title("Plasma density along beam axis")
    plt.ylabel("n")
    plt.xlabel("Radius from beam axis (mm)")
    plt.show()
    
    plt.plot(yrange_z,n_vs_y)
    plt.title("Plasma density along beam axis")
    plt.ylabel("n")
    plt.xlabel("Radius from beam axis (mm)")
    plt.plot(yrange_z,n_vs_y1)
    plt.plot(yrange_z,n_vs_y2)
    plt.plot(yrange_z,n_vs_y3)
    plt.show()
    
    plt.plot(yrange_z,n_vs_y2)
    plt.title("Plasma density along beam axis")
    plt.ylabel("n")
    plt.xlabel("Radius from beam axis (mm)")
    plt.plot(yrange_z,n_vs_y2d)
    plt.show()
    
    compare=0
    if compare ==1:
        den2=np.empty([len(xrange),len(yrange),len(zrange)])
        for x in i_counter:
            for y in j_counter:
                for z in k_counter:
                    den2[x][y][z]=H[x][y][z]*dendis.SimpleGasJet(nozzle_den,L_z,L_r,2*nozzle_loc,xrange[x],yrange[y],zrange[z])
                    n_vs_z2   = den2[round(len(xrange_z)/2),round(len(yrange_z)/2),:]
                    plt.plot(zrange_z,n_vs_z)
                    plt.plot(zrange_z,n_vs_z2)
                    plt.title("Plasma density along gas jet")
                    plt.ylabel("n")
                    plt.xlabel("Radius from laser axis (mm)")
                    plt.show()
                
                