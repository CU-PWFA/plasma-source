#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 18:21:12 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
density= 5
power= 7
pw= 70
#%%
PlasmaExpan= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy')
Ext= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_SimExt.npy')
Box= Ext[0]
GridSize= Ext[1]
SimTime= Ext[2]
Frame= Ext[3]

GridNum= PlasmaExpan.shape[0]
FrameNum= PlasmaExpan.shape[1]

t= np.linspace(0, SimTime*1e9, FrameNum)
r= np.linspace(-Box, Box, GridNum)
#%%
PhotonMap= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
x= np.linspace(-350, 350, PhotonMap.shape[0])
Imagrey= np.meshgrid(x, PhotonMap[:, int(PhotonMap.shape[1]-1)])[1]
#%%
plt.pcolormesh(x, r, Imagrey/np.max(Imagrey))
plt.xlabel('x ($\mu$m)')
plt.ylabel('Lineout ($\mu$m)')
plt.colorbar(label= 'Normalized Intensity')

#%%
plt.pcolormesh(t, r, PlasmaExpan, cmap='plasma')
plt.xlabel('Time (ns)')
plt.ylabel('Lineout ($\mu$m)')
plt.colorbar(label= 'Density')
