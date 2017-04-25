# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 00:39:44 2017

@author: chris
"""

import GaussianBeam
import numpy as np
import matplotlib.pyplot as plt

wavelength = 100e-9
k = np.pi*2/(wavelength)
w0 = np.sqrt(0.159*wavelength/2)*0.4
zR = np.power(w0,2)*np.pi/wavelength

q_init = complex(1e-25,zR)
q=[q_init]

#4 Thin Lens
"""
q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
q=GaussianBeam.Prop_ThinLens(q,10)
q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
q=GaussianBeam.Prop_ThinLens(q,10)
q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
q=GaussianBeam.Prop_ThinLens(q,10)
q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
q=GaussianBeam.Prop_ThinLens(q,10)
q=GaussianBeam.Prop_FreeSpace(q,5,0.1)
xrange = np.arange(0,25,.1)
"""

q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)
q=GaussianBeam.Prop_FlatInterface(q,1,1.5)
q=GaussianBeam.Prop_FreeSpace(q,2.6,0.1)
xrange = np.arange(0,5,.1)

R=GaussianBeam.Prop_RadiusList(q)
w1=GaussianBeam.Prop_SpotList(q[:25],wavelength,1)
w2=GaussianBeam.Prop_SpotList(q[26:],wavelength,1.5)
w=list(w1)
w.extend(w2)

plt.plot(xrange,w)
plt.title("Spot size from q parameter")
plt.ylabel("Spot Size (m)")
plt.xlabel("Distance (m)")
plt.show()

plt.plot(xrange[1:],R[1:-1],)
plt.title("Curvature from q parameter")
plt.ylabel("Radius of Curvature (m)")
plt.xlabel("Distance (m)")
plt.show()