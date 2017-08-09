#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:49:29 2017

@author: mike
"""

import numpy as np
import beam_ana as ba

#step = 2244
step = 0
#step = 709

lips = ba.calc_95ellipse(vbeam,step)

print(lips)

[u,v] = ba.real2norm_coords(ebeam[step]["x"],ebeam[step]["xp"],\
                            ebeam[step]["beta"],ebeam[step]["alpha"])

figX, axX1 = plt.subplots(1,1,sharey=True)
Xs = axX1.scatter(u,v,c='r',s=0.5)

cent = lips["center"]
rot  = lips["rot"]

pt1 = cent+[lips["radii"][0],0]
r1 = np.dot(pt1,rot)
r1mag = np.sqrt(r1[0]**2+r1[1]**2)
print(r1)
Xr1 = axX1.scatter(r1[0],r1[1],c='k',s=2)

pt2 = cent+[0,lips["radii"][1]]
r2 = np.dot(pt2,rot)
r2mag = np.sqrt(r2[0]**2+r2[1]**2)
print(r2)
Xr2 = axX1.scatter(r2[0],r2[1],c='k',s=2)

focmag = np.sqrt((r1mag/2)**2+(r2mag/2)**2)
if np.abs(r2mag>r1mag):
    theta = np.arctan2(r2[1],r2[0])
else:
    theta = np.arctan2(r1[1],r1[0])
foc1 = [focmag*np.cos(theta),focmag*np.sin(theta)]
foc2 = [focmag*np.cos(theta+np.pi),focmag*np.sin(theta+np.pi)]
Xf1 = axX1.scatter(foc1[0],foc1[1],c='b',s=2)
Xf2 = axX1.scatter(foc2[0],foc2[1],c='b',s=2)

def calc_s(u,v,foc1,foc2):
    s = np.sqrt((u-foc1[0])**2+(v-foc1[0])**2) +\
        np.sqrt((u-foc2[0])**2+(v-foc2[0])**2)
    return s

s = np.zeros(len(u))
for i in range(len(u)):
    s[i] = calc_s(u[i],v[i],foc1,foc2)

i_s = np.argsort(s)
i_s = i_s[:-round(0.05*len(i_s))]

u2 = u[i_s]
v2 = v[i_s]

print(len(u))
print(len(u2))



