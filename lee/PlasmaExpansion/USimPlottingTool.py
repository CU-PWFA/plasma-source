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

#%%
#Run in the USim folder
BasePath= 'Diffusion_2D/'
BaseFileName= 'cap_2D_'
#%%
MassDen={}
Eng={}
den={}
temp={}
den_3d=np.zeros((50, 50, 11))
for n in range (0, 11):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (50, 50))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS)/(100)**3
    den_3d[:, :, n]= den[n]
    Eng[n]= np.reshape(qReader(path) [1], (50, 50))
    temp[n]=Eng[n]/n0/kb*2
    

#%%
time= np.linspace(0, 2.8899418245468987E-8, 11)
#%%
plt.pcolormesh(den[10], cmap='plasma')
plt.title('Plasma Density @ t='+str(time[10]))
plt.colorbar()
plt.axis('scaled')
plt.savefig('den_t='+str(time[10])+".png")

#%%
# standard matplotlib stuff
# create the different plotting axes
x= np.linspace(0, 50, 50)
y= np.linspace(0, 50, 50)
f= np.linspace(0, 11, 11)

X, Y= np.meshgrid(x, y)
pcolormesh_data= den_3d

fig, ax=plt.subplots()

ax.set_aspect('equal')
ax.set_xlabel('x(mm)')
ax.set_xlabel('x(mm)')

fig.suptitle('Plasma Density @ t='+str(time[10]))
pcolormesh_block = amp.blocks.Pcolormesh(X, Y, pcolormesh_data, t_axis=2, cmap='plasma')

plt.colorbar(pcolormesh_block.quad)
timeline = amp.Timeline(f)

# now to contruct the animation
anim = amp.Animation(pcolormesh_block, timeline)
#%%
anim.controls()

#anim.save_gif('images/multiblock')
plt.show()

#%%
x= np.linspace(-800, 800, 50)
y= np.linspace(-800, 800, 50)
f= np.linspace(0, 28.899418245468987, 11)

X, Y, F = np.meshgrid(x, y, f)

block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], den_3d, t_axis=2, cmap='plasma')
cbar= plt.colorbar(block.quad)
cbar.set_label('cm^-3')
plt.title('Plasma Density [time]=ns')
plt.gca().set_aspect('equal')
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

timeline = amp.Timeline(f, fps=1.5)
anim = amp.Animation([block], timeline)

anim.controls()

anim.save_gif('pcolormesh')
plt.show()