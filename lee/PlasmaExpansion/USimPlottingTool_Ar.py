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
from scipy import optimize

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
#BasePath= 'EPExchange/'
#BasePath= 'EPExchange_100ns/'
BasePath= 'EPExchange_Ar_20ns/'
BaseFileName= 'cap_2D_'
#%%
MassDen={}
eng={}
den={}
temp={}
visc={}
GridSize=150
Frame=21
Box= 100

den_3d=np.zeros((GridSize, GridSize, Frame))
en_3d=np.zeros((GridSize, GridSize, Frame))
temp_3d=np.zeros((GridSize, GridSize, Frame))
visc_3d=np.zeros((GridSize, GridSize, Frame))

for n in range (0, Frame):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (GridSize, GridSize))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS*36)/(100)**3
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
LineOutArray= np.zeros((GridSize, Frame))
for n in range (0, Frame):
    LineOutArray[:, n]= den[n][:, int(GridSize/2)]
#%%
plt.figure(1)
plt.title('')
plt.pcolormesh(time, r, LineOutArray)#, norm= colors.LogNorm(vmin=LineOutArray.min(), vmax=LineOutArray.max()))
plt.colorbar()
plt.ylabel('LineOut $\mu$m')
plt.xlabel('Time (ns)')
plt.title('Argon')
#%%
plt.plot(LineOutArray[75, :])
plt.xlabel('Time (ns)')
plt.ylabel('Density (cm^-3)')
#%%
Coef= np.zeros((21, 6))
#%%
t=11

x= np.linspace(0, Box, int(GridSize/2))
y= LineOutArray[int(LineOutArray[:, 0].shape[0]/2):int(LineOutArray[:, 0].shape[0]), t]/1e17
#%%
fitfn= lambda p, x: p[0]*np.exp(-((x)**p[2])/(2*p[1]**p[2]))+p[3]*np.exp(-((x-p[5])**2)/(2*p[4]**2))+\
                    p[6]*np.exp(-((x)**p[8])/(2*p[7]**p[8]))

#fitfn= lambda p, x: p[0]*np.exp(-((x)**p[2])/(2*p[1]**p[2]))

errfunc = lambda p, x, y: fitfn(p, x) - y
#p0= [0.7, 40, 1]
#p0= [0.1, 0.1, 0.1, 3, 1, 30, 1, 40, 2]
#p0=[ 0, 63.69525901,  5.76333434,  4.28537451,  1.7140818,  38.09384361, -2.57863762, 66.14103163,  6.75383988]
p0=[ 0.2, 6300.69525901,  5.76333434,  4.28537451,  1.7140818,  38.09384361, 0, 66.14103163,  6.75383988]
#p0= [0.7, 40, 1, 3, 1, 30, 1, 40, 2]
#p0= [0, 0, 2, 5, 2.5, 38, 0, 0, 0]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y), epsfcn= 0.01)
print(p1)
Coef[t, :]= p1
#%%
t=20

x= np.linspace(0, Box, int(GridSize/2))
x= np.append(x[0:15], x[45:75])
y= np.append(LineOutArray[int(LineOutArray[:, 0].shape[0]/2):int(LineOutArray[:, 0].shape[0]/2*1.2), t]/1e17, \
                          LineOutArray[120:int(LineOutArray[:, 0].shape[0]), t]/1e17)

#fitfn= lambda p, x: p[0]*np.exp(-((x)**p[2])/(2*p[1]**p[2]))#+p[3]*np.exp(-((x-p[5])**2)/(2*p[4]**2))+\
 #                   p[6]*np.exp(-((x)**p[8])/(2*p[7]**p[8]))
fitfn= lambda p, x: p[0]*np.exp(-((x)**p[2])/(2*p[1]**p[2]))+p[3]*np.exp(-((x)**p[5])/(2*p[4]**p[5]))

#fitfn= lambda p, x: p[0]*np.exp(-((x)**p[2])/(2*p[1]**p[2]))

errfunc = lambda p, x, y: fitfn(p, x) - y
p0= [0.1, 30, 2, 0.1, 20, 2]
#p0=[2.86811761e-03, 5.51810118e+01, 1.14461051e+02, 1.15694979e-01, 2.70032198e+01, 2.01745158e+00]
#p0= [0.1, 0.1, 0.1, 3, 1, 30, 1, 40, 2]
#p0=[ 0, 63.69525901,  5.76333434,  4.28537451,  1.7140818,  38.09384361, -2.57863762, 66.14103163,  6.75383988]
#p0=[ 0.2, 6300.69525901,  5.76333434,  4.28537451,  1.7140818,  38.09384361, 0, 66.14103163,  6.75383988]
#p0= [0.7, 40, 1, 3, 1, 30, 1, 40, 2]
#p0= [0, 0, 2, 5, 2.5, 38, 0, 0, 0]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y), epsfcn= 0.01)
print(p1)
Coef[t, :]= p1
#%%
#Coef[6, :]= [0, 0, 0, 0, 0, 0]
#%%
xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.plot(x, y, '.')#, label= 'plasma data')
plt.plot(xNew, fitfn(p1, xNew), '-', label='t=10ns')
#plt.plot(xNew, fitfn(p1, xNew), '-', label='fitting')
plt.legend()
plt.xlabel('r ($\mu$m)')
plt.ylabel('Density (1e17 cm^-3)')
#plt.title('t= 0ns')
     #%%
np.save('Coef2', Coef)
#%%
plt.figure(0)
levs= (1e13, 5e13, 1e14, 5e14, 1e15, 5e15, 1e16, 5e16, 1e17)
plt.contour(time, r, LineOutArray, levs, norm= colors.LogNorm(vmin=LineOutArray.min(), vmax=LineOutArray.max()))
plt.colorbar()
plt.ylabel('LineOut $\mu$m')
plt.xlabel('Time (s)')

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