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
t= np.linspace(0, 10e-3, Ntimestep)
#t= np.linspace(0, 5e-6, Ntimestep)
den0= 1e23
#%%
alpha= 2.36e-14*1e-6 #m^3/s
#alpha= 3e-16#5e-18 #m^3/s
dndt=-alpha*den0**2
dt= np.amax(t)/(Ntimestep-1)
dn= dndt*dt
#%%
nNew=np.zeros(Ntimestep-1)
pn= np.zeros(Ntimestep)
den= den0
for time in range (0, Ntimestep-1):
    n= den+dn
    dndt=-alpha*n**2
    dt= np.amax(t)/(Ntimestep-1)
    dn= dndt*dt
    pn[time]=dn 
    den= n
    nNew[time]= n

#%%
a518= nNew
#a516= nNew
#%%
plt.plot((t+dt)[0:Ntimestep-1]*1e3, a518, label= '$Recombination Coefficient \\alpha$= 2.36e-14 $(cm^3/s)$')
#plt.plot((t+dt)[0:Ntimestep-1]*1e6, a516, label= 'alpha8=5e-10 (cm^3/s)')

plt.xlabel('Time (ms)')
plt.ylabel('Density ($m^{-3}$)')
plt.legend()
#plt.title('Alpha=1e-8cm^3/s_TimeStep=10000')
