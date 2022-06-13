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

#%%
#Run in the USim folder
BasePath= 'Diffusion_2D/'
BaseFileName= 'cap_2D_'
#%%
MassDen={}
eng={}
den={}
temp={}
den_3d=np.zeros((100, 100, 21))
en_3d=np.zeros((100, 100, 21))
temp_3d=np.zeros((100, 100, 21))

for n in range (0, 21):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (100, 100))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS)/(100)**3
    den_3d[:, :, n]= den[n]
    eng[n]= np.reshape(qReader(path) [1], (100, 100))
    en_3d[:, :, n]= eng[n]
    temp[n]= np.reshape(TempReader(path), (100, 100))
    temp_3d[:, :, n]= temp[n]/11604.45

#    temp[n]=Eng[n]/n0/kb*2
    

#%%
time= np.linspace(0, 2.8899418245468987E-8, 11)
#%%
plt.pcolormesh(temp[17], cmap='plasma')
plt.title('Plasma Density @ t='+str(time[17]))
plt.colorbar()
plt.axis('scaled')
plt.savefig('den_t='+str(time[10])+".png")


#%%
temp_3d_1= temp_3d
#%%
x= np.linspace(-1500, 1500, 100)
y= np.linspace(-1500, 1500, 100)
f= np.linspace(0, 28.899418245468987, 21)

X, Y, F = np.meshgrid(x, y, f)

block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], en_3d, t_axis=2, cmap='plasma')#, vmin=0, vmax=5)
cbar= plt.colorbar(block.quad)
cbar.set_label('J')

plt.title('Totoal Energy [time]=ns')
plt.gca().set_aspect('equal')
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

timeline = amp.Timeline(f, fps=1.5)
anim = amp.Animation([block], timeline)

anim.controls()

anim.save_gif('pcolormesh')
plt.show()