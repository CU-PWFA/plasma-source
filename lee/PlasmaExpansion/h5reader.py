#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:31:39 2020

@author: valentina_lee
"""
#%%
import numpy as np
import h5py
import matplotlib.pyplot as plt

#%% Call data
def CDReader(FileName):

    hf=h5py.File(FileName,'r')
    CD= hf['fluids']['Cdensity']
    T=np.array(CD.get('timeMesh'))
    density= np.array(CD.get('data'))
    TimeDens=np.zeros((density.shape[0], 2), dtype=complex)
    TimeDens[:, 0]= T[:, 0]
    TimeDens[:, 1]= density[:, 0]
    return TimeDens

#%%
#Run in the USim folder
#BasePath=  'Diffusion_LastFrame/'
BasePath=  'DiffusionWithDF_smallgrid/'
BaseFileName= 'cap_2D_'
#%%
data= np.zeros((5000, 2), dtype= complex) #if see bugs on the data shape, check if the initial data shape is too small
LastShape=0
#%%
for n in range (1, 31):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    TempData= CDReader(path)
    data[int(LastShape):int(LastShape+TempData.shape[0]), 0]=TempData[:, 0]
    data[int(LastShape):int(LastShape+TempData.shape[0]), 1]=TempData[:, 1]
    LastShape= TempData.shape[0]+LastShape
#%%
Data= data[0:LastShape, :]
del data
#%%
d=Data[:, 1]
d_max= np.amax(d)
den= d/d_max*1e17
#%%
plt.plot(LFX, LFden, label='Diffuse to Neutral Gas_LastFrame') 
#plt.plot(abs(Data[:, 0]*1e9),den, label='Diffuse to Vacuum') 
plt.plot(ComplexX,ComplexDen, label='Diffuse to Neutral Gas')
plt.xlabel('t(ns)')
plt.ylabel('density(cm^-3)')
plt.legend()

#plt.title('ComplexModel')
#%%
ComplexDen=den
ComplexX=abs(Data[:, 0]*1e9 )
#%%
LFDen=den
LFX=abs(Data[:, 0]*1e9 )

#%%
np.save('Den.npy', den)
np.save('T.npy', Data[:, 0])
