#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 19:12:53 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy import optimize

#%%
ELECMASS = 9.10938356e-31     # kG
PROTMASS = 1.672621898e-27    # kG
#%%
AGdata= np.load('AGdata.npy')
#%%
GridNum= 150
Box= 600#1500
FrameNum= 30
SimTime= 100
NewTimeGrid= 10000
IntResolution= 1e-9#s
FnlTimeStepNum= round(SimTime/IntResolution)

x= np.linspace(-Box, Box, GridNum)
y= np.linspace(-Box, Box, GridNum)
r= x[int(GridNum/2):GridNum]
X, Y= np.meshgrid(x, y)
R= np.sqrt(X**2+Y**2)

#%%
den0= 1e23
TimeR= np.linspace(0, SimTime, r.shape[0])
xvar, yvar= np.meshgrid(TimeR, r)

def fitfn(p, xvar, yvar):
    (p[0]*np.exp(-(xvar*np.cos(-p[5])-yvar*np.sin(-p[5]))**2/(2*p[1]**2))+\
    p[2]*np.exp(-((xvar*np.cos(-p[6])-yvar*np.sin(-p[6]))**3/(2*p[3]**2)**3)))* \
     ((-np.tanh((xvar*np.sin(-p[7])+yvar*np.cos(-p[7]))/p[4])+1)/np.amax(-np.tanh((xvar*np.sin(-p[7])+yvar*np.cos(-p[7]))/p[4])+1))

def us(x):
    return (np.sign(x)+1)/2;

def errfunc(p):
    TimeR= np.linspace(0, SimTime, r.shape[0])
    xvar, yvar= np.meshgrid(TimeR, r)
    
    fit= us(p)
    fit= fitfn(p, xvar, yvar)
    f= interp2d(TimeR, r, fit, kind='cubic')
    Time= np.linspace(0, SimTime, NewTimeGrid)
    LineOutTime= f(Time, r)*den0

    Ntimestep= int(Time.shape[0])
    alpha= 3e-16#5e-18 #m^3/s
    PhotonNum= np.zeros((int(GridNum/2), int(NewTimeGrid-1)))

    for f in range (0, int(GridNum/2)):
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
                    
        PhotonNum[f, :]= -pn[1:int(NewTimeGrid)]
        
        AfterglowTime= np.zeros((int(GridNum), int(FnlTimeStepNum)))

    for sec in range (0, int(FnlTimeStepNum)):
        IntPhoNum= np.zeros(int(GridNum/2))
        for f in range (0, int(GridNum/2)):
            frame= PhotonNum[f, :]
            IntPhoNum[f]= np.sum(frame[0:int(NewTimeGrid/FnlTimeStepNum)+sec*int(NewTimeGrid/FnlTimeStepNum)])

        xvar= r
        yvar= IntPhoNum
        fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)

        errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
        p0= [10e20, 200, 0];
        p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
        AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
        AfterglowSide= np.sum(AfterglowMap, axis= 0)
        AfterglowTime[:, sec]= AfterglowSide

    error=AfterglowTime[: int(FnlTimeStepNum-1)]-AGdata

    return error        
#%%
p0= [0.6, 6, 1.4, 3.3, 10, 18, 0.5, -1.6]
p1, success = optimize.leastsq(errfunc, p0[:], epsfcn= 2.5e-20);
print(p1)

#%%
plt.figure(2)
Timeplt= np.linspace(0, SimTime*1e9, int(FnlTimeStepNum))
plt.pcolormesh(Timeplt, y_, SAfterglow)
plt.xlabel('Integration time (ns)')
plt.ylabel('Lineout ($\mu$m)')
plt.colorbar(label= '# of photons')
plt.title('w= 50 $\mu$m, n=1e17cm^-3')