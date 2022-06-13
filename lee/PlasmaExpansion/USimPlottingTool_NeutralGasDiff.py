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
    return Den

#%%
#Run in the USim folder
BasePath= 'NeutralGasDiff30ns/'
BaseFileName= 'diffusion_'
#%%
den={}
den_3d=np.zeros((100, 100, 31))
for n in range (0, 31):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    den[n]= (np.reshape(qReader(path) , (100, 100))*1e-6)#[10:90, 10:90]
    den_3d[:, :, n]= den[n]

#%%
den={}
den_3d=np.zeros((100, 100, 21))
for n in range (0, 21):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    den[n]= (np.reshape(qReader(path) , (100, 100))*1e-6)
    den_3d[:, :, n]= den[n]


#%%
time= np.linspace(0, 30, 21)
#%%
plt.figure(2)
plt.pcolormesh(den[20], cmap='plasma')
plt.title('Plasma Density @ t='+str(time[20]))
plt.colorbar()
plt.axis('scaled')
plt.savefig('den_t='+str(time[20])+".png")


#%%
x= np.linspace(-800, 800, 100)
y= np.linspace(-800, 800, 100)
f= np.linspace(0, 30, 21)

X, Y, F = np.meshgrid(x, y, f)

block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], den_3d, t_axis=2, cmap='plasma')#, vmin=-1.4e15, vmax=1e17)
cbar= plt.colorbar(block.quad)
cbar.set_label('cm^-3')
plt.title('Neutral Density [time]=ns')
plt.gca().set_aspect('equal')
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')    

timeline = amp.Timeline(f, fps=1.5)
anim = amp.Animation([block], timeline)

anim.controls()

anim.save_gif('pcolormesh')
plt.show()

#%%
LastFrameLO= den[20][50, :]

plt.plot(LastFrameLO)
#%%
#adjust the power to find a best fit 
y= LastFrameLO
fitfn= lambda p, x: -p[0]*np.exp(-1/2*(x/p[1])**2)-p[2]*np.exp(-(x**2/(2*p[3]**2))**2)+1e17
errfunc = lambda p, x, y: fitfn(p, x) - y;
p0= [1e16, 200, 1e16, 200];
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y), epsfcn= 2.5e-34);
print(p1)
#%%
plt.plot(x, y, '.', label= 'plasma data')
plt.plot(x, fitfn(p1, x), '-', label='fitting')
plt.legend()
plt.xlabel('$\mu$m')
plt.ylabel('cm^-3')
#plt.text(250, 7e16, 'f(x)=', fontsize=14)
#plt.title('f(x)=-p1[0]*np.exp(-1/2*(x/p1[1])**2)-p1[2]*np.exp(-(x**2/(2*p1[3]**2))**2)+1e17'+str(f'{p1[0]:.1e}'))
