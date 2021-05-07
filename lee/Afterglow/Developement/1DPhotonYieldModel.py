#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:41:33 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
h=6.63e-34
lambda0=800e-9
c=3e8
f=lambda0/c
#%%create a 1D density decay plot
Ntimestep=10000
t= np.linspace(0, 50, Ntimestep)
den= (-np.tanh(t/10)+1)*1e23
#t100000=t
#%%
alpha= 1e-15#5e-18 #m^3/s
n=den[0]
dndt=-alpha*n**2
dt= 50/(Ntimestep-1)*1e-9
dn= dndt*dt
#%%
nNew=np.zeros(Ntimestep)
#nNew[0]=n
pn= np.zeros(Ntimestep)
for time in range (0, Ntimestep-1):
    pn[time]=dn 
    n= den[time+1]+dn
    dndt=-alpha*n**2
    dt= 50/(Ntimestep-1)*1e-9
    dn= dndt*dt
    nNew[time]=n
#%%
#pn100=pn
#pn1000=pn
#pn10000=pn
pn100000=pn
#%%Photon Yield
Energy=-pn*h*f
#%%
plt.plot(t, den, (t+dt), nNew)
plt.xlabel('ns')
plt.ylabel('m^-3')
plt.title('Alpha=1e-8cm^3/s_TimeStep=10000')

#%%
#plt.plot(t100, pn100, label='100')
#plt.plot(t1000, pn1000,  label= '1000')
#plt.plot(t10000, pn10000, label= '10000')
#plt.plot(t100000, pn100000, label='100000')
plt.plot(t, pn)
plt.xlabel('ns')
plt.ylabel('m^-3')
plt.title('dn_Alpha=1e-8cm^3/s')
#plt.legend()

#%%
plt.plot(t, Energy)
plt.xlabel('ns')
plt.ylabel('J')
plt.title('Energy_Alpha=1e-7cm^3/s')
