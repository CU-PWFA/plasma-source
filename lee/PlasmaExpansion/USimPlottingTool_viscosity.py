#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:09:49 2020

@author: valentina_lee
"""

#%%
import numpy as np
import h5py
import matplotlib.pyplot as plt
import animatplot as amp
import matplotlib.colors as colors
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
#%% Call data
def TempReader(FileName):

    hf=h5py.File(FileName,'r')
    Flu= hf['fluids']
    Temp=np.array(Flu.get('Temperature'))
    return Temp
#%% Call data
def ViscReader(FileName):

    hf=h5py.File(FileName,'r')
    Flu= hf['fluids']
    Visc=np.array(Flu.get('Viscosity'))
    return Visc


#%%
#Run in the USim folder
#BasePath= 'Diffusion_NoDF/'
#BasePath= 'NoDF_100ns/'
#BasePath= 'NoDF_20ns/'
BasePath= 'EPExchange_Ar_Debug1/'
#BasePath= 'EPExchange_100ns/'
#BasePath= 'EPExchange_100ns/'
BaseFileName= 'cap_2D_'
#%%
MassDen={}
eng={}
den={}
temp={}
visc={}
GridSize=150
Frame=11
Box= 100

den_3d=np.zeros((GridSize, GridSize, Frame))
en_3d=np.zeros((GridSize, GridSize, Frame))
temp_3d=np.zeros((GridSize, GridSize, Frame))
visc_3d=np.zeros((GridSize, GridSize, Frame))

for n in range (0, Frame):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (GridSize, GridSize))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS*4)/(100)**3
#    den[n]=MassDen[n]/(ELECMASS+PROTMASS*4)/(100)**3
    den_3d[:, :, n]= den[n]
    eng[n]= np.reshape(qReader(path) [1], (GridSize, GridSize))
    en_3d[:, :, n]= eng[n]
#    temp[n]= np.reshape(TempReader(path), (GridSize, GridSize))
#    temp_3d[:, :, n]= temp[n]/11604.45
#    visc[n]= np.reshape(ViscReader(path), (GridSize, GridSize))
#    visc_3d[:, :, n]= visc[n]

#%%
r= np.linspace(-Box, Box, GridSize)
time= np.linspace(0, 20, Frame)
#%%
plt.subplot(3, 3, 1)
f=int(0)
plt.pcolormesh(r, r, den[f], cmap='plasma')
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

#%%
LineOutArray= np.zeros((GridSize, Frame))
for n in range (0, Frame):
    LineOutArray[:, n]= den[n][:, int(GridSize/2)]
#%%
plt.figure(1)
plt.title('')
plt.pcolormesh(time, r, LineOutArray, cmap='plasma')#, norm= colors.LogNorm(vmin=LineOutArray.min(), vmax=LineOutArray.max()))
plt.colorbar(label='Density ($cm^{-3}$)')
plt.ylabel('LineOut ($\mu$m)')
plt.xlabel('Time (ns)')
plt.title('Plasma Expansion by Fluid Simulation')
#plt.title('EPExchange_20ns')

#%%
plt.figure(0)
levs= (1e13, 5e13, 1e14, 5e14, 1e15, 5e15, 1e16, 5e16, 1e17)
plt.contour(time, r, LineOutArray, levs, norm= colors.LogNorm(vmin=LineOutArray.min(), vmax=LineOutArray.max()))
plt.colorbar()
plt.ylabel('LineOut $\mu$m')
plt.xlabel('Time (ns)')

#%%
plt.plot(LineOutArray[:, 30])
#%%    
np.save('LineOutArray_FirstFrame', LineOutArray)
#%%
x= np.linspace(-Box, Box, GridSize)
y= np.linspace(-Box, Box, GridSize)
f= np.linspace(0, 28.899418245468987, 31)

X, Y, F = np.meshgrid(x, y, f)

block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], den_3d, t_axis=2, cmap='plasma')#, vmin=0, vmax=3.5e-4)
cbar= plt.colorbar(block.quad)
cbar.set_label('cm^-3')

plt.title('Density [time]=ns NG_LastFrame')
plt.gca().set_aspect('equal')

timeline = amp.Timeline(f, fps=1.5)
anim = amp.Animation([block], timeline)

anim.controls()

anim.save_gif('pcolormesh')
plt.show()