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
#BasePath= 'HOFIwithVis/'
BasePath= 'HOFI/'
#BasePath= 'HOFIStudy2/2.4e24/'
BaseFileName= 'cap_2D_'
#%%
MassDen={}
eng={}
den={}
temp={}
visc={}
GridSize=250
Frame=41
Box= 500

den_3d=np.zeros((GridSize, GridSize, Frame))
en_3d=np.zeros((GridSize, GridSize, Frame))
temp_3d=np.zeros((GridSize, GridSize, Frame))
visc_3d=np.zeros((GridSize, GridSize, Frame))

for n in range (0, Frame):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (GridSize, GridSize))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS)/(100)**3
    den_3d[:, :, n]= den[n]
    eng[n]= np.reshape(qReader(path) [1], (GridSize, GridSize))
    en_3d[:, :, n]= eng[n]
    temp[n]= np.reshape(TempReader(path), (GridSize, GridSize))
    temp_3d[:, :, n]= temp[n]/11604.45
    visc[n]= np.reshape(ViscReader(path), (GridSize, GridSize))
    visc_3d[:, :, n]= visc[n]

#%%
r= np.linspace(-Box, Box, GridSize)
time= np.linspace(0, 4e-9, Frame)
#%%
plt.subplot(3, 3, 1)
f=int(1)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 2)
f=int(4)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 3)
f=int(8)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 4)
f=int(12)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 5)
f=int(16)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 6)
f=int(20)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 7)
f=int(24)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 8)
f=int(28)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

plt.subplot(3, 3, 9)
f=int(31)
plt.pcolormesh(r, r, den[f], cmap='jet', vmax=2.5e18, vmin=0)
plt.title('Plasma Density @ t='+str("{:.2e}".format(time[f])))
plt.colorbar()
plt.axis('scaled')
#plt.savefig('den_t='+str(time[f])+".png")
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')
#%%
plt.plot(r[int(GridSize/2):int(GridSize)], den[0][int(GridSize/2):int(GridSize), int(GridSize/2)],label='0.1ns')
plt.plot(r[int(GridSize/2):int(GridSize)], den[1][int(GridSize/2):int(GridSize), int(GridSize/2)], label='0.2ns')
plt.plot(r[int(GridSize/2):int(GridSize)], den[2][int(GridSize/2):int(GridSize), int(GridSize/2)], label='0.3ns')
plt.plot(r[int(GridSize/2):int(GridSize)], den[3][int(GridSize/2):int(GridSize), int(GridSize/2)], label='0.4ns')
plt.plot(r[int(GridSize/2):int(GridSize)], den[4][int(GridSize/2):int(GridSize), int(GridSize/2)], label='0.5ns')
plt.plot(r[int(GridSize/2):int(GridSize)], den[5][int(GridSize/2):int(GridSize), int(GridSize/2)], label='0.6ns')
plt.plot(r[int(GridSize/2):int(GridSize)], den[6][int(GridSize/2):int(GridSize), int(GridSize/2)], label='0.7ns')
plt.legend()
plt.xlabel('LineOut $\mu$m')
plt.ylabel('Density (cm^-3)')


#%%
#den1e24= den
#den3e23= den
#den4_5e23= den
#den4e23= den
#den5_5e23= den
#den7e23= den
den8e23= den
#%%
#plt.plot(r[int(GridSize/2):int(GridSize)], den1e23[29][int(GridSize/2):int(GridSize), int(GridSize/2)], label='1e17 4.2mbar')
#plt.plot(r[int(GridSize/2):int(GridSize)], den5e23[29][int(GridSize/2):int(GridSize), int(GridSize/2)], label='NoVis')
#plt.plot(r[int(GridSize/2):int(GridSize)], den1e24[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='WithVis')
plt.plot(r[int(GridSize/2):int(GridSize)], den3e23[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='3e17 25mbar')
plt.plot(r[int(GridSize/2):int(GridSize)], den4_5e23[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='4.5e17 35mbar')
plt.plot(r[int(GridSize/2):int(GridSize)], den5_5e23[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='5.5e17 45mbar')
plt.plot(r[int(GridSize/2):int(GridSize)], den7e23[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='7e17 55mbar')
plt.plot(r[int(GridSize/2):int(GridSize)], den8e23[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='8e17 65mbar')
plt.plot(r[int(GridSize/2):int(GridSize)], den1e24[31][int(GridSize/2):int(GridSize), int(GridSize/2)], label='1e18 75mbar')
plt.legend()
plt.xlabel('LineOut $\mu$m')
plt.ylabel('Density (cm^-3)')


#%%
LineOutArray= np.zeros((GridSize, Frame))
#for n in range (0, Frame):
for n in range (0, 40):
    LineOutArray[:, n]= den[n][:, int(GridSize/2)]
#%%
rt= ((5/3)+1)**0.5*(13.7*1.6e-19/PROTMASS)**0.25*(50e-6*(time+0e-10))**0.5*1e6
plt.figure(1)
plt.title('')
plt.pcolormesh(time, r, LineOutArray)#, norm= colors.LogNorm(vmin=LineOutArray.min(), vmax=LineOutArray.max()))
plt.colorbar()
plt.ylabel('LineOut $\mu$m')
plt.xlabel('Time (s)')
#plt.plot(time, rt, label= 'Sedov-Taylor Theory')
#plt.legend()

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