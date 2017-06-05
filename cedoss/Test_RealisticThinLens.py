# -*- coding: utf-8 -*-
"""
Created on Fri May 26 09:55:45 2017

@author: chris
"""

import GaussianBeam
import numpy as np
import matplotlib.pyplot as plt
import ADK_Combined as adk

c=2.998e8
P=4e9
cm_m=100
wavelength = 785.3e-9
w0 = 5e-3

I0 = 2*P/(np.pi*np.power(w0,2))*np.power(1/cm_m,2)
delt_t = 100e-15
chi = adk.Get_chi_H()

q=GaussianBeam.Prop_Init_q(wavelength, w0, -0.5, 1)
q=GaussianBeam.Prop_FreeSpace(q,0.5,1e-4)
#q=GaussianBeam.Prop_ThickLens(q,1.4,1e-3,8,-8,1e-4)
q=GaussianBeam.Prop_ThinLens(q,.5)
q=GaussianBeam.Prop_FreeSpace(q,1,1e-4)

xrange_tot=GaussianBeam.Prop_GetRange(q)
w_tot=GaussianBeam.Prop_SpotList(q,wavelength)
"""
print(min(w_tot))
print(GaussianBeam.Gauss_zR(min(w_tot),wavelength))

#Plots spot size
plt.plot(xrange_tot,w_tot)
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.show()
"""
waist=w_tot.index(min(w_tot))
w=w_tot[(waist-100):(waist+100)]
xrange=xrange_tot[(waist-100):(waist+100)]
"""
#Plots spot size zoomed in around waist
plt.plot(xrange,w)
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.show()
"""
#Create a 2D array for the intensity of this guassian beam as a function of
# both w(z) and r
rrange=np.arange(-0.0002,0.0002,0.000001)
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

waist_off1_H = H[w.index(min(w))+10]
waist_off2_H = H[w.index(min(w))+20]
waist_off3_H = H[w.index(min(w))+30]
waist_off4_H = H[w.index(min(w))+40]
"""
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

plt.plot([i*100 for i in xrange],center_H)
plt.title("Ion. Frac. along propagation axis")
plt.ylabel("ni/n0")
plt.xlabel("Distance (cm)")
plt.show()

plt.plot(rrange,waist_H)
plt.title("Ion. Frac. through waist")
plt.ylabel("ni/n0")
plt.xlabel("Radius from axis (m)")
plt.show()
"""
middle=round(len(rrange)/2)
width=50
rrange_middle=rrange[middle-width:middle+width]*1000
waist_H_middle=waist_H[middle-width:middle+width]
waist_off1_H_middle=waist_off1_H[middle-width:middle+width]
waist_off2_H_middle=waist_off2_H[middle-width:middle+width]
waist_off3_H_middle=waist_off3_H[middle-width:middle+width]
waist_off4_H_middle=waist_off4_H[middle-width:middle+width]

#top=np.where(waist_H_middle > 0.95)
#print(rrange_middle[top[0][0]])

plt.plot(rrange_middle,waist_H_middle)
plt.title("Ion. Frac. through waist (zoomed)")
plt.ylabel("ni/n0")
plt.xlabel("Radius from axis (mm)")
plt.plot(rrange_middle,waist_off1_H_middle)
plt.plot(rrange_middle,waist_off2_H_middle)
plt.plot(rrange_middle,waist_off3_H_middle)
plt.plot(rrange_middle,waist_off4_H_middle)
plt.show()
"""
#Plot intensity
plt.imshow(I, interpolation="none", origin="lower",
           extent=[rrange[0],rrange[-1],xrange[0],xrange[-1]], aspect='.05')
CB=plt.colorbar()
CB.set_label('I (I0 = 1e10)')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('Intensity')
plt.show()

plt.imshow(H, interpolation="none", origin="lower",
           extent=[rrange[0],rrange[-1],xrange[0],xrange[-1]], aspect='.05')
CB=plt.colorbar()
CB.set_label('ni/n0')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('ADK Ion.')
plt.show()
"""