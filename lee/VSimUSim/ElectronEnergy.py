#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:59:00 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import optimize

#%% Call data
def VSimFileReader(FileName, Species):

    hf=h5py.File(FileName,'r')
    matrix=np.array(hf.get(Species))
    return matrix
#%%
#BasePath= 'CompleteIonization20eV/'
#BaseFileName= 'IonizeFinalResolved_'
#BasePath= 'CompleteIonization15eV/'
#BaseFileName= 'ValentinaNewLaser_'
BasePath= 'LaserProfile_Min/'
BaseFileName= 'ValentinaAugustLow_'
#%%
n=5
Species= 'Electrons'
Filename= BaseFileName+Species+'_'+str(n)+'.h5'
path= BasePath+Filename
matrix_e= VSimFileReader(path, 'Electrons')
xloc_e= matrix_e[:, 0]
yloc_e= matrix_e[:, 1]
Eng_e= matrix_e[:, 5]

#%%
plt.figure(figsize=(10, 10))
plt.scatter(xloc_e, yloc_e, c=Eng_e)
plt.axis("image")
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.colorbar()
#%%
r_e= np.sqrt(xloc_e**2+yloc_e**2)
Tavg= sum(Eng_e)/(Eng_e.shape[0])
print(Tavg)
#%%
plt.plot(r_e, Eng_e, '.')
plt.xlabel('r (m)')
plt.ylabel('Temperature (eV)')

#%%
#BasePath= 'KatsResults/'
#BaseFileName= 'IonizeFinalResolved_'
BasePath= 'CompleteIonization15eV/'
BaseFileName= 'ValentinaNewLaser_'
#%%
n=5
Species= 'HePlus'
Filename= BaseFileName+Species+'_'+str(n)+'.h5'
path= BasePath+Filename
matrix_i= VSimFileReader(path, 'HePlus')
xloc_i= matrix_i[:, 0]
yloc_i= matrix_i[:, 1]
#Eng_i= matrix_i[:, 5]
r_i= np.sqrt(xloc_i**2+yloc_i**2)
#%%
plt.figure(figsize=(10, 10))
plt.scatter(xloc_i, yloc_i, c=Eng_i)
plt.axis("image")
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.colorbar()
#%%
plt.plot(r_i, Eng_i, '.')
plt.xlabel('r (m)')
plt.ylabel('Temperature (eV)')
#%%20eV
Bin=1000
den= np.zeros(Bin)
hist= np.histogram(r_i, Bin)
dr= hist[1][1]-hist[1][0]
r= hist[1]+dr/2
cellv= (140e-6/4000)**2*1
PPC=1
NPIM= 1e23*cellv/PPC
den1= hist[0]*NPIM/(2*np.pi*r[0:int(Bin)]*dr)
#%%15eV
Bin=1000
den= np.zeros(Bin)
hist= np.histogram(r_i, Bin)
dr= hist[1][1]-hist[1][0]
r= hist[1]+dr/2
cellv= (140e-6/2000)**2*1
PPC=1
NPIM= 1e23*cellv/PPC
den1= hist[0]*NPIM/(2*np.pi*r[0:int(Bin)]*dr)
#%%
plt.plot(den1)
#%%
plt.hist(r_i, 30)
plt.xlabel('r (m)')
plt.ylabel('counts')
#%%
plt.plot(hist[1][0:int(Bin)], den1, label= 'VSim')
plt.xlabel('r (m)')
plt.ylabel('Density (m^-3)')
plt.legend()
#%%
#%%
x= hist[1][0:int(Bin)]
y= den1
fitfn= lambda p, x: p[0]*np.exp(-(x**p[2]/(2*p[1]**p[2])))

errfunc = lambda p, x, y: fitfn(p, x) - y;
p0= [1e23, 20e-6, 2];
#p0= [1e17, 20, 1e16, 200, 9];
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)
#%%
#plasma= 1e23*np.exp(-x**8/(2*27.7e-6**8))
plasma15= 1e23*np.exp(-x**2/(2*22.3e-6**2))
#%%
#plt.plot(x, y, '.', label= 'VSim result')
plt.plot(x, fitfn(p1, x), '-', label='VSim')
plt.plot(x, plasma15, label= 'Split-Step result')
plt.legend()

