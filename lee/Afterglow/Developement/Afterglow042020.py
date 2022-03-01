vvvv#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 17:41:18 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp2d
from scipy import optimize

#%%
ELECMASS = 9.10938356e-31     # kG
PROTMASS = 1.672621898e-27    # kG
#%% Call data
def qReader(FileName):

    hf=h5py.File(FileName,'r')
    Flu= hf['fluids']
    q=np.array(Flu.get('q'))
    Den=q[:, 0]
    Eng=q[:, 4]
    return Den, Eng
#%%
BasePath= 'DiffusionWithDF_smallgrid/'
BaseFileName= 'cap_2D_'
#%%
x= np.linspace(-1500, 1500, 200)
y= np.linspace(-1500, 1500, 200)
r= x[100:200]
X, Y= np.meshgrid(x, y)
R= np.sqrt(X**2+Y**2)
#%%
MassDen={}
eng={}
den={}

den_3d=np.zeros((200, 200, 31))
en_3d=np.zeros((200, 200, 31))

for n in range (0, 31):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (200, 200))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS*4)
    den_3d[:, :, n]= den[n]
    eng[n]= np.reshape(qReader(path) [1], (200, 200))
    en_3d[:, :, n]= eng[n]
#%%
LineoutTime= np.zeros((100, 31))
for n in range (0, 31):
    LineoutTime[:, n]= den[n][100:200, 100]
#%%
time= np.linspace(0, 100e-9, 31)
f= interp2d(time, r, LineoutTime, kind='cubic')
Time= np.linspace(0, 100e-9, 10000)
LineOutTime= f(Time, r)
#%%
plt.figure(0)
plt.pcolormesh(LineOutTime)
plt.xlabel('time frame')
plt.ylabel('r')
plt.title('Density Lineout - Time')

 #%%create a 1D density decay plot
Ntimestep= int(Time.shape[0])
#alpha= 5e-18 #m^3/s
alpha= 3e-16#5e-18 #m^3/s
PhotonNum= np.zeros((100, 9999))
for f in range (0, 100):
    den= LineOutTime[f, :]
    n0=den[0]
    dndt=-alpha*n0**2
    dt= np.amax(Time)/(Ntimestep-1)
    dn= dndt*dt
    nNew=np.zeros(Ntimestep-1)
    pn= np.zeros(Ntimestep)
    Den=n0
    for t in range (0, Ntimestep-1):
        pn[t]=dn 
        n= Den+(den[t+1]-den[t])+dn
        dndt=-alpha*n**2
        dn= dndt*dt
        if n>0 :
            nNew[t]=n
            Den= n
        else:
            nNew[t]=0
            Den= 0
        
    PhotonNum[f, :]= -pn[1:10000]
       
#%%
plt.pcolormesh(PhotonNum)
plt.xlabel('time frame')
plt.ylabel('r')
plt.title('Radial Photon Number - Time')
#%%Plot the recombination density change
plt.plot(Time, den, (Time+dt)[0:9999], nNew)
plt.xlabel('s')
plt.ylabel('m^-3')
#plt.yscale('log')

#plt.title('Alpha=1e-9cm^3/s_TimeStep~1000')

#%%Run from "Here"
#%%
IntPhoNum= np.zeros(100)
for f in range (0, 100):
    frame= PhotonNum[f, :]
    IntPhoNum[f]= np.sum(frame[0:30])
#%%
plt.plot(r, IntPhoNum)
plt.xlabel('r ($\mu$m)')
plt.ylabel('Photon Number')
#%%
xvar= r
yvar= IntPhoNum
fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)

errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
p0= [10e20, 200, 0];
p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))#, epsfcn= 2.5e-34);
print(p1)
#%%
plt.plot(xvar, yvar, '.', label= 'plasma data')
plt.plot(xvar, fitfn(p1, xvar), '-', label='fitting')
plt.legend()

#%%
AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
#%%
plt.pcolormesh(x, y, AfterglowMap)
plt.axis('scaled')
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

#%%
AfterglowSide= np.sum(AfterglowMap, axis= 0)
plt.plot(AfterglowSide)

#%%"Here"
#%%Afterglow at each time sec
AfterglowTime= np.zeros((200,100))
for sec in range (0, 100):
    IntPhoNum= np.zeros(100)
    for f in range (0, 100):
        frame= PhotonNum[f, :]
        IntPhoNum[f]= np.sum(frame[0+sec*100:100+sec*100])

    xvar= r
    yvar= IntPhoNum
    fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)

    errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
    p0= [10e20, 200, 0];
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
    AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
    AfterglowSide= np.sum(AfterglowMap, axis= 0)
    AfterglowTime[:, sec]= AfterglowSide

#%%
plt.pcolormesh(AfterglowTime)
#%%Afterglow Integral over each time
AfterglowTime= np.zeros((200,100))
for sec in range (0, 100):
    IntPhoNum= np.zeros(100)
    for f in range (0, 100):
        frame= PhotonNum[f, :]
        IntPhoNum[f]= np.sum(frame[0:100+sec*100])

    xvar= r
    yvar= IntPhoNum
    fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)

    errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
    p0= [10e20, 200, 0];
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
    AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
    AfterglowSide= np.sum(AfterglowMap, axis= 0)
    AfterglowTime[:, sec]= AfterglowSide

#%%
Timeplt= np.linspace(0, 100, 100)
plt.pcolormesh(Timeplt, x, AfterglowTime)
plt.xlabel('Integration time (ns)')
plt.ylabel('Lineout ($\mu$m)')
plt.colorbar(label= '# of photons')
