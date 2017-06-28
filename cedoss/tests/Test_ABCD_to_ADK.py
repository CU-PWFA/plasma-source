# -*- coding: utf-8 -*-
"""
Created on Thu May 18 15:59:33 2017

Takes the ABCD method of propagating a gaussian beam through an optical
system.  Then the intensity of the beam is found, relative to the intensity
of the beam at the original beam waist.  Then the intensity is related
to the ADK model for ionization of a gas

@author: chris
"""
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "../")

from modules import GaussianBeam
from modules import Doss_Ionization as adk

wavelength = 500e-9
k=2*np.pi/wavelength
w0 = 1e-3
I0 = 0.6e14
delt_t = 100e-15
chi = adk.Get_chi_H()

q=GaussianBeam.Prop_Init_q(wavelength, w0, -6, 1)
q=GaussianBeam.Prop_FreeSpace(q,6,1e-3)
#q=GaussianBeam.Prop_ThickLens(q,1.4,1e-3,8,-8,1e-4)
q=GaussianBeam.Prop_ThinLens(q,3)
q=GaussianBeam.Prop_FreeSpace(q,8,1e-3)

xrange=GaussianBeam.Prop_GetRange(q)
w=GaussianBeam.Prop_SpotList(q,wavelength)

#Plots spot size
plt.plot(xrange,w)
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.show()

#Create a 2D array for the intensity of this guassian beam as a function of
# both w(z) and r
rrange=np.arange(-0.0015,0.0015,0.00001)
I=np.empty([len(xrange),len(rrange)])
i_counter=np.arange(len(xrange))
j_counter=np.arange(len(rrange))
for i in i_counter:
    for j in j_counter:
        I[i][j]=GaussianBeam.GaussianBeamIntensity_SpotArray(I0,w[i],w0,rrange[j])

H=np.empty([len(xrange),len(rrange)])
for i in i_counter:
    H[i] = adk.TotalProb(adk.Elefield(I[i]),delt_t,chi,1)

#Get the Intensity and Ion. Frac. through the axis and across the waist,
# and plot them
center_I = [middle[round(len(rrange)/2)] for middle in I]
waist_I = I[w.index(min(w))]

center_H = [middle[round(len(rrange)/2)] for middle in H]
waist_H = H[w.index(min(w))]

plt.plot(xrange,center_I)
plt.title("Intensity along propagation axis")
plt.ylabel("Intensity (W/cm^2)")
plt.xlabel("Distance (m)")
plt.show()

plt.plot(rrange,waist_I)
plt.title("Intensity through waist")
plt.ylabel("Intensity (W/cm^2)")
plt.xlabel("Radius from axis (m)")
plt.show()

plt.plot(xrange,center_H)
plt.title("Ion. Frac. along propagation axis")
plt.ylabel("ni/n0")
plt.xlabel("Distance (m)")
plt.show()

plt.plot(rrange,waist_H)
plt.title("Ion. Frac. through waist")
plt.ylabel("ni/n0")
plt.xlabel("Radius from axis (m)")
plt.show()

plt.plot(rrange[100:130],waist_H[100:130])
plt.title("Ion. Frac. through waist (zoomed)")
plt.ylabel("ni/n0")
plt.xlabel("Radius from axis (m)")
plt.show()

#Plot intensity
plt.imshow(I, interpolation="none", origin="lower",
           extent=[rrange[0],rrange[-1],xrange[0],xrange[-1]], aspect='.0002')
CB=plt.colorbar()
CB.set_label('I (I0 = 0.6e14')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('Intensity')
plt.show()

plt.imshow(H, interpolation="none", origin="lower",
           extent=[rrange[0],rrange[-1],xrange[0],xrange[-1]], aspect='.0002')
CB=plt.colorbar()
CB.set_label('ni/n0')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('ADK Ion.')
plt.show()