#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:46:00 2020

@author: valentinalee
"""

#%%
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
#%%
ELECMASS = 9.10938356e-31     # kG
PROTMASS = 1.672621898e-27    # kG
n0=1e23 
kb=1.38064852e-23  # J/K 
#%% Call data
def qReader(FileName):
    hf=h5py.File(FileName,'r')
    Flu= hf['fluids']
    q=np.array(Flu.get('q'))
    Den=q[:, 0]
    Eng=q[:, 4]
    return Den, Eng
#%%
def q_meshReader(FileName):
    hf= h5py.File(FileName,'r')
    Flu= hf['fluids']['domain']['face']
    cells= np.array(Flu.get('cells'))
    vertices= np.array(Flu.get('vertices'))
    return cells, vertices
#%%
BasePath= 'n1_GP1_PW30/'
#BasePath= 'EPExchange_He_Debug1/'
#BasePath= 'EPExchange_Ar_Debug1/'
BaseFileName= 'cap_2D_'

#%%
GridSize=800
Frame=41
Box= 250 #1e-6m

grid_x= np.linspace(-Box*1e-6, Box*1e-6, GridSize)
grid_y= np.linspace(-Box*1e-6, Box*1e-6, GridSize)
grid_X, grid_Y= np.meshgrid(grid_x, grid_y)

den_3d=np.zeros((GridSize, GridSize, Frame))

Gridpath= BasePath+ 'cap_2DGrid.h5'
cells, vertices= q_meshReader(Gridpath)

for n in range (0, Frame):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    print(path)
    den= qReader(path) [0]/(ELECMASS+PROTMASS*18*2)/(100)**3
    grid=np.zeros((len(den), 2))

    for cell in range (0, len(den)):
        corner1x= vertices[cells[cell, 0], 0]
        corner1y= vertices[cells[cell, 0], 1]
        corner2x= vertices[cells[cell, 1], 0]
        corner2y= vertices[cells[cell, 1], 1]
        corner3x= vertices[cells[cell, 2], 0]
        corner3y= vertices[cells[cell, 2], 1]
        corner4x= vertices[cells[cell, 3], 0]
        corner4y= vertices[cells[cell, 3], 1]
        
        grid[cell, 0]= (corner1x+corner2x+corner3x+corner4x)/4
        grid[cell, 1]= (corner1y+corner2y+corner3y+corner4y)/4

    den_3d[:, :, n] = griddata(grid, den, (grid_X, grid_Y), method='linear')     
    #%%
r= np.linspace(-Box, Box, GridSize)
time= np.linspace(0, 60, Frame)

#%%
LineOutArray= np.zeros((GridSize, Frame))
for n in range (0, Frame):
    LineOutArray[:, n]= den_3d[:, :, n][:, int(GridSize/2)]
#%%
plt.figure(1)
plt.title('')
plt.pcolormesh(time, r, LineOutArray)#, norm= colors.LogNorm(vmin=LineOutArray.min(), vmax=LineOutArray.max()))
plt.colorbar()
plt.ylabel('LineOut $\mu$m')
plt.xlabel('Time (ns)')
plt.title('n5_GP9_PW70')
#plt.savefig('')
#plt.title('Helium')
#%%
plt.plot(LineOutArray[75, :])
plt.xlabel('Time (ns)')
plt.ylabel('Density (cm^-3)')
