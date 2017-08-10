#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:49:29 2017

@author: mike
"""

import numpy as np
import beam_ana as ba
import matplotlib.pyplot as plt

#step = 2244
step = 0
#step = 293 # vac beta min
#step = 709

# get normalized coords
#[u,v] = ba.real2norm_coords(ebeam[step]["x"],ebeam[step]["xp"],\
#                            ebeam[step]["beta"],ebeam[step]["alpha"])

u = vbeam[step]["x"]
y = vbeam[step]["xp"]

lips = ba.calc_95ellipse(u,v,step)
rad  = lips["radii"]
cent = lips["center"]
rot  = lips["rot"]

# plot data points
figX, axX1 = plt.subplots(1,1,sharey=True)
axX1.scatter(u,v,c='r',s=0.5)
axX1.set_xlim([-1.1*max(abs(u)),+1.1*max(abs(u))])
axX1.set_ylim([-1.1*max(abs(v)),+1.1*max(abs(v))])

# define ellipse axes
elax1 = [rad[0],0]
elax1 = np.dot(elax1,rot)
elax1mag = np.sqrt(elax1[0]**2+elax1[1]**2)

elax2 = [0,rad[1]]
elax2 = np.dot(elax2,rot)
elax2mag = np.sqrt(elax2[0]**2+elax2[1]**2)

axX1.plot([-elax1[0],+elax1[0]]+cent[0],[-elax1[1],+elax1[1]]+cent[1],c='k')
axX1.plot([-elax2[0],+elax2[0]]+cent[0],[-elax2[1],+elax2[1]]+cent[1],c='k')

# define ellipse foci
focmag = np.sqrt((elax1mag/2)**2+(elax2mag/2)**2)
if np.abs(elax2mag>elax1mag):
    theta = np.arctan2(elax2[1],elax2[0])
else:
    theta = np.arctan2(elax1[1],elax1[0])
foc1 = [focmag*np.cos(theta)+cent[0],focmag*np.sin(theta)+cent[1]]
foc2 = [focmag*np.cos(theta+np.pi)+cent[0],focmag*np.sin(theta+np.pi)+cent[1]]

axX1.scatter(foc1[0],foc1[1],c='b',s=20)
axX1.scatter(foc2[0],foc2[1],c='b',s=20)

# plot ellipse
phi = np.linspace(0.0, 2.0 * np.pi, 1000)
elx = rad[0] * np.cos(phi)
ely = rad[1] * np.sin(phi)

for i in range(len(elx)):
    [elx[i],ely[i]] = np.dot([elx[i],ely[i]],rot) + cent

axX1.plot(elx,ely,c='b')
